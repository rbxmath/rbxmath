local Tools = require("Tools")

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

_complex.__tostring = function (complex)
    return tostring(complex[1]) .. " + " .. "sqrt(-1) " .. tostring(complex[2])
end

local _rational = {}

local _rationalFromArray = function (array)
    if array[2] == 0 then
        error("Division by zero!", -1)
    end

    local result = {array[1], array[2]}

    setmetatable(result, _rational)

    return result
end

local _rationalFromTwoInputs = function (a, b)
    if b == 0 then
        error("Division by zero!", -1)
    end

    local result = {a, b}

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

local _ZmodPInverse = function (zinteger)
    
end