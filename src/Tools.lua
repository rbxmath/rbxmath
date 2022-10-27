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
        continuedFractionArray[i] = math.floor(nArray[i] / nArray[i + 1])
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

local _subarray = function (array, index)
    local result = {}

    for i = 1, math.min(index, #array) do
        result[i] = array[i]
    end

    return result
end

local _suparray = function (array, index)
    local result = {}

    for i = index, #array do
        result[#result+1] = array[i]
    end

    return result
end

local _binarySearch
_binarySearch = function (element, array)
    if array == nil or #array == 0 then
        return false
    end

    local checkIndex = math.floor(#array / 2)

    if #array == 1 then
        if element == array[1] then
            return 1
        else
            return false
        end
    end

    if element < array[checkIndex] then
        return _binarySearch(element, _subarray(array, checkIndex - 1))
    elseif element > array[checkIndex] then
        local intermediate = _binarySearch(element, _suparray(array, checkIndex + 1))
        if intermediate == false then
            return false
        else
            return checkIndex + intermediate
        end
    else
        return checkIndex
    end
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

Tools.list = {}

Tools.list.binarySearch = function (element, sortedArray)
    return _binarySearch(element, sortedArray)
end

Tools.list.sublist = function (array, index)
    return _subarray(array, index)
end

Tools.list.suplist = function (array, index)
    return _suparray(array, index)
end

Tools.list.tostring = function (array)
    local result = "{"
    for i = 1, #array-1, 1 do
        result = result .. tostring(array[i]) .. ", "
    end
    result = result .. tostring(array[#array]) .. "}"

    return result
end

Tools.list.error = function (left, right)
    if #left ~= #right then
        error("Incomparable lists!", -1)
    end

    local error = 0;

    for i = 1, #left, 1 do
        error = math.max(math.abs(left[i] - right[i]), error)
    end

    return error
end

Tools.list.map = function (f, list)
    local result = {}

    for i, v in ipairs(list) do
        result[i] = f(v)
    end

    return result
end

Tools.list.linspace = function (a, b, n)
    local result = {}

    for i = 0, n-1, 1 do
        result[i + 1] = a + (b - a) * i / (n - 1)
    end
end

return Tools