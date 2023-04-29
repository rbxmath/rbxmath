local Tools = require("src.Tools")

local Scalars = {}

local _complex = {}

local _complexFromArray = function (array)
    local result = {array[1], array[2]}

    setmetatable(result, _complex)

    return result
end

local _complexFromTwoInputs = function (a, b)
    local result = {a, b}

    setmetatable(result, _complex)

    return result
end

_complex.__eq = function (left, right)
    return left[1] == right[1] and left[2] == right[2]
end

_complex.__add = function (left, right)
    return _complexFromTwoInputs(left[1] + right[1], left[2] + right[2])
end

_complex.__sub = function (left, right)
    return _complexFromTwoInputs(left[1] - right[1], left[2] - right[2])
end

_complex.__mul = function (left, right)
    return _complexFromTwoInputs(left[1] * right[1] - left[2] * right[2], left[1] * right[2] + left[2] * right[1])
end

_complex.__unm = function(complex)
    return _complexFromTwoInputs(-complex[1], -complex[2])
end

local _complexConjugate = function (complex)
    return _complexFromTwoInputs(complex[1], complex[2])
end

local _complexNorm = function (complex)
    return (complex * _complexConjugate(complex))[1]
end

local _complexTimesScalar = function (scalar, complex)
    return _complexFromTwoInputs(scalar * complex[1], scalar * complex[2])
end

local _complexInverse = function (complex)
    return _complexTimesScalar(1 / _complexNorm(complex), _complexConjugate(complex))
end

_complex.__div = function (left, right)
    return left * _complexInverse(right)
end

local _signString = function (x)
    if x >= 0 then 
        return " + "
    else
        return " - "
    end
end

_complex.__tostring = function (complex)
    local sign = function (x)
        if x >= 0 then 
            return " + "
        else
            return " - "
        end
    end
    return tostring(complex[1]) .. _signString(complex[2]) .. "i " .. tostring(math.abs(complex[2]))
end

local _rational = {}

local _rationalFromArray = function (array)
    if array[2] == 0 then
        error("Division by zero!", -1)
    end

    local result = {math.tointeger(array[1]), math.tointeger(array[2])}

    setmetatable(result, _rational)

    return result
end

local _rationalFromTwoInputs = function (a, b)
    if b == 0 then
        error("Division by zero!", -1)
    end

    local result = {math.tointeger(a), math.tointeger(b)}

    setmetatable(result, _rational)

    return result
end

_rational.__eq = function (left, right)
    return left[1] == right[1] and left[2] == right[2]
end

_rational.__add = function (left, right)
    local denominator = left[2] * right[2]
    local numerator = left[1] * right[2] + right[1] * left[2]

    local divisor = Tools.integers.gcd(denominator, numerator)

    return _rationalFromTwoInputs(numerator / divisor, denominator / divisor)
end

_rational.__sub = function (left, right)
    local denominator = left[2] * right[2]
    local numerator = left[1] * right[2] - right[1] * left[2]

    local divisor = Tools.integers.gcd(denominator, numerator)

    return _rationalFromTwoInputs(numerator / divisor, denominator / divisor)
end

_rational.__mul = function (left, right)
    return _rationalFromTwoInputs(left[1] * right[1], left[2] * right[2])
end

_rational.__div = function (left, right)
    return _rationalFromTwoInputs(left[1] * right[2], left[2] * right[1])
end

_rational.__tostring = function (rational)
    return tostring(rational[1]) .. "/" .. tostring(rational[2])
end

local _rationalFromContinuedFraction = function (continuedFraction)
    local n = #continuedFraction
    local rational = _rationalFromTwoInputs(1, continuedFraction[n])

    for i = n - 1, 1, -1 do
        rational = _rationalFromTwoInputs(continuedFraction[i], 1) + _rationalFromTwoInputs(rational[2], rational[1])
    end

    return rational
end

local _rationalApproximation = function (float, tolerance)
    local n = -math.floor(math.log(tolerance, 10) + 0.5)

    if n <= 0 then
        return math.floor(float)
    end

    local continuedFraction = Tools.continuedFractions.compute(float, n)

    local i = 1;

    while math.abs(float - Tools.continuedFractions.evaluate(continuedFraction)) > tolerance and i < math.maxinteger do
        continuedFraction = Tools.continuedFractions.compute(float, n + i)
        i = i + 1
    end

    return _rationalFromContinuedFraction(continuedFraction)
end

local _ZmodP = {}

local _ZmodPFromArray = function (array)
    local result = {array[1] % array[2], array[2]}

    setmetatable(result, _ZmodP)

    return result
end

local _ZmodPFromTwoInputs = function (a, b)
    local result = {a % b, b}

    setmetatable(result, _ZmodP)

    return result
end

_ZmodP.__eq = function (left, right)
    return left[1] == right[1] and left[2] == right[2]
end

_ZmodP.__add = function (left, right)
    if left[2] ~= right[2] then
        error("Modulus mismatch!", -1)
    end

    return _ZmodPFromTwoInputs(left[1] + right[1], left[2])
end

_ZmodP.__sub = function (left, right)
    if left[2] ~= right[2] then
        error("Modulus mismatch!", -1)
    end

    return _ZmodPFromTwoInputs(left[2] + left[1] - right[1], left[2])
end

_ZmodP.__mul = function (left, right)
    if left[2] ~= right[2] then
        error("Modulus mismatch!", -1)
    end

    return _ZmodPFromTwoInputs(left[1] * right[1], left[2])
end

_ZmodP.__tostring = function (zinteger)
    return tostring(zinteger[1])
end

--- Credit to Wikipedia https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
local _ZmodPInverse = function (zinteger)
    local t = 0
    local r = zinteger[2]
    local newt = 1
    local newr = zinteger[1]

    local quotient
    while newr ~= 0 do
        quotient = r % newr
        t, newt = newt, t - quotient * newt
        r, newr = newr, r - quotient * newr
    end

    if r > 1 then
        error("Input is not invertible!", -1)
    end

    return _ZmodPFromTwoInputs(t + zinteger[2], zinteger[2])
end

_ZmodP.__div = function (left, right)
    return left * _ZmodPInverse(right)
end

Scalars.Complex = {}

Scalars.Complex.new = function (a, b)
    if type(a) == "table" then
        return _complexFromArray(a)
    elseif b == nil then
        return _complexFromTwoInputs(a, 0)
    else
        return _complexFromTwoInputs(a, b)
    end
end

Scalars.Complex.newFromArray = function (array)
    return _complexFromArray(array)
end

Scalars.Complex.norm = function (complex)
    return _complexNorm(complex)
end

Scalars.Complex.conjugate = function (complex)
    return _complexConjugate(complex)
end

Scalars.Complex.inverse = function (complex)
    return _complexInverse(complex)
end

Scalars.Complex.sort = function (array, order)
    local comp = function (left, right)
        if order == "d" then
            if type(left) == "number" and type(right) == "number" then
                return math.abs(left) > math.abs(right)
            elseif getmetatable(left) == getmetatable(right) then
                return _complexNorm(left) > _complexNorm(right)
            elseif getmetatable(left) == _complex then
                return _complexNorm(left) > math.abs(right)
            else
                return math.abs(left) > _complexNorm(right)
            end

        else
            if type(left) == "number" and type(right) == "number" then
                return math.abs(left) < math.abs(right)
            elseif getmetatable(left) == getmetatable(right) then
                return _complexNorm(left) < _complexNorm(right)
            elseif getmetatable(left) == _complex then
                return _complexNorm(left) < math.abs(right)
            else
                return math.abs(left) < _complexNorm(right)
            end
        end
    end
    
    table.sort(array, comp)
end

Scalars.Rational = {}

Scalars.Rational.new = function (a, b)
    return _rationalFromTwoInputs(a, b)
end

Scalars.Rational.newFromArray = function (array)
    return _rationalFromArray(array)
end

Scalars.Rational.toReal = function (rational)
    return rational[1] / rational[2]
end

Scalars.Rational.rationalApproximation = function (float, tolerance)
    return _rationalApproximation(float, tolerance)
end

Scalars.Rational.rationalFromContinuedFraction = function (continuedFraction)
    return _rationalFromContinuedFraction(continuedFraction)
end

Scalars.Rational.inverse = function (rational)
    return _rationalFromTwoInputs(rational[2], rational[1])
end

Scalars.ZmodP = {}

Scalars.ZmodP.new = function (a, b)
    return _ZmodPFromTwoInputs(a, b)
end

Scalars.ZmodP.newFromArray = function (array)
    return _ZmodPFromArray(array)
end

Scalars.ZmodP.inverse = function (zinteger)
    return _ZmodPInverse(zinteger)
end

return Scalars