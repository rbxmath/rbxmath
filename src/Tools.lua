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

Tools.integers.gcd = function (a, b)
    return _gcd(a, b)
end

return Tools