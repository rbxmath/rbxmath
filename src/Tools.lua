local Tools = {}

local _primeFactorization = function (a)
    local endPoint = math.floor(math.sqrt(a))
end

local _gcd
_gcd = function (a, b)
    if a % 1 ~= 0 and b % 1 ~= 0 then
        error("Inputs are not integers!", -1)
    end

    if b == 0 then
        return a
    end

    return _gcd(b, a % b)
end

local _integerPower = function (a, b)
    local integer = 1

    for i = 1, b, 1 do
        integer = integer * a
    end

    return integer
end

local _computeContinuedFraction = function (float, n)
    local continuedFractionArray = {}
    local nArray = {float, 1}

    for i = 1, n, 1 do
        continuedFractionArray[i] = math.floor[nArray[i] / nArray[i + 1]]
        nArray[i + 2] = nArray[i] % nArray[i + 1]
        if nArray[i + 1] == 0 then
            break
        end
    end

    return continuedFractionArray
end

local _evaluateContinuedFraction = function (continuedFraction)
    local n = #continuedFraction
    local result = continuedFraction[n]

    for i = n - 1, 1, -1 do
        result = continuedFraction[i] + 1 / result
    end

    return result
end

Tools.integers = {}

Tools.integers.gcd = function (a, b)
    return _gcd(a, b)
end

Tools.integers.power = function (a, b)
    return _integerPower(a, b)
end

Tools.continuedFractions = {}

Tools.continuedFractions.compute = function (float, n)
    return _computeContinuedFraction(float, n)
end

Tools.continuedFractions.evaluate = function (continuedFraction)
    return _evaluateContinuedFraction(continuedFraction)
end

return Tools