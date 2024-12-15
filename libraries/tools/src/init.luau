--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = {}

export type Array<T> = { [number] : T }
export type Vector = { [number] : number }
export type Object = { [any] : any }
export type ScalarMap = (number) -> number
export type Tensor = { [number] : { [number] : any } }

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

    for _ = 1, b, 1 do
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

local _binarySearchBetween
_binarySearchBetween = function (element, array)
    if array == nil or #array == 0 then
        return 0
    end

    local checkIndex = math.floor(#array / 2)

    if #array == 1 then
        if element == array[1] then
            return {1}
        elseif element < array[1] then
            return 1
        else
            return 2
        end
    end

    if element < array[checkIndex] then
        return _binarySearchBetween(element, _subarray(array, checkIndex - 1))
    elseif element > array[checkIndex] then
        local intermediate = _binarySearchBetween(element, _suparray(array, checkIndex + 1))
        if type(intermediate) == "table" then
            return {intermediate[1] + checkIndex}
        else
            return checkIndex + intermediate
        end
    else
        return {checkIndex}
    end
end

local function _regulaFalsi (f, t, a, b, tol)
    local tolerance = tol or 10^(-13)
    local leftValue, rightValue, middleValue = f(a) - t, f(b) - t, 0
    local left, right, middle = a, b, (a + b) / 2
    if leftValue * rightValue > 0 then
        return nil
    elseif math.abs(leftValue) < tolerance then
        return left
    elseif math.abs(rightValue) < tolerance then
        return right
    end

    while math.abs(leftValue - rightValue) >= tolerance and math.abs(left - right) >= tolerance do
        middle = (right * leftValue - left * rightValue)/(leftValue - rightValue)
        middleValue = f(middle) - t
        if math.abs(middleValue) < tolerance or math.abs(left - middle) / math.abs(leftValue - rightValue) < tolerance or math.abs(right - middle) / math.abs(leftValue - rightValue) < tolerance then
            return middle
        elseif math.abs(leftValue) < tolerance then
            return left
        elseif math.abs(rightValue) < tolerance then
            return right
        elseif leftValue * middleValue > 0 then
            leftValue = middleValue
            left = middle
        else
            rightValue = middleValue
            right = middle
        end
    end
    
    return middle
end

Tools.integers = {}

function Tools.integers.gcd (a, b)
    return _gcd(a, b)
end

function Tools.integers.power (a, b)
    return _integerPower(a, b)
end

Tools.continuedFractions = {}

function Tools.continuedFractions.compute (float, n)
    return _computeContinuedFraction(float, n)
end

function Tools.continuedFractions.evaluate (continuedFraction)
    return _evaluateContinuedFraction(continuedFraction)
end

Tools.list = {}

function Tools.list.binarySearch (element, sortedArray)
    return _binarySearch(element, sortedArray)
end

Tools.list.binarySearchBetween = _binarySearchBetween

function Tools.list.sublist (array, index)
    return _subarray(array, index)
end

function Tools.list.suplist (array, index)
    return _suparray(array, index)
end

function Tools.list.join (left, right)
    local new = {}
    for _, value in ipairs(left) do
        new[#new + 1] = value
    end
    for _, value in ipairs(right) do
        new[#new + 1] = value
    end
    return new
end

function Tools.list.tostring (array)
    local result = "{"
    for i = 1, #array-1, 1 do
        result = result .. tostring(array[i]) .. ", "
    end
    result = result .. tostring(array[#array]) .. "}"

    return result
end

function Tools.list.deeptostring(array: Object): string
   local result = "{"
   for i = 1, #array - 1, 1 do
      if type(array[i]) == "table" then
	 result = result .. Tools.list.deeptostring(array[i]) .. ", "
      else
	 result = result .. tostring(array[i]) .. ", "
      end
   end
   if type(array[#array]) == "table" then
      result = result .. Tools.list.deeptostring(array[#array]) .. "}"
   else
      result = result .. tostring(array[#array]) .. "}"
   end

   return result
end

function Tools.list.error (left, right)
    if #left ~= #right then
        error("Incomparable lists!", -1)
    end

    local error = 0;

    for i = 1, #left, 1 do
        error = math.max(math.abs(left[i] - right[i]), error)
    end

    return error
end

function Tools.list.map (f, list)
    local result = {}

    for i, v in ipairs(list) do
        result[i] = f(v)
    end

    return result
end

function Tools.list.linspace (a, b, n)
    local result = {}

    for i = 0, n-1, 1 do
        result[i + 1] = a + (b - a) * i / (n - 1)
    end

    return result
end

function Tools.list.copy (array)
    local result = {}

    for i = 1, #array, 1 do
        result[i] = array[i]
    end

    return result
end

function Tools.list.norm (array)
    local sum = 0
    for _, value in ipairs(array) do
        sum = sum + math.pow(value, 2)
    end
    return math.sqrt(sum)
end

function Tools.list.reverse (array)
    local data = {}
    local n = #array
    for i = n, 1, -1 do
        data[n - i + 1] = array[i]
    end
    return data
end

function Tools.list.scale (array, scale)
    local data = {}
    local n = #array
    for i = 1, n, 1 do
        data[i] = scale * array[i]
    end
    return data
end

Tools.solve = {}

function Tools.solve.regulaFalsi (f, t, a, b, tol)
    return _regulaFalsi(f, t, a, b, tol)
end

Tools.combinatorics = {}

function Tools.combinatorics.inversionNumber (permutation)
    local n = #permutation
    local inversionNumber = 0
    for i = 1, n - 1 do
        for j = i + 1, n do
            if permutation[i] > permutation[j] then
                inversionNumber = inversionNumber + 1
            end
        end
    end
    return inversionNumber
end

function _padStringToLength (string, length)
    local result = string
    if (length - #result) % 2 == 1 then
        result = result .. " "
    end
    local width = length - #result
    for _ = 1, width / 2 do
        result = " " .. result .. " "
    end
    return result
end

function _padStringRightToLength (string, length)
    local result = string
    local width = length - #result
    for _ = 1, width do
        result = result .. " "
    end
    return result
end

Tools.admin = {}

function Tools.admin.makeBanners (header, body, width)
    width = width or 50
    local result = "+"
    for _ = 1, width do
        result = result .. "-"
    end
    result = result .. "+\n"
    result = result .. "|"
    result = result .. _padStringToLength(header, width)
    result = result .. "|\n"
    result = result .. "+"
    for _ = 1, width do
        result = result .. "-"
    end
    result = result .. "+\n"
    local iter = 1
    local space, lastSpace = 0, 0
    while iter <= #body do
        lastSpace = space
        space = string.find(body, " ", lastSpace + 1)
        if space == nil then
            if iter + width - 1 >= #body then
                result = result .. "|" .. _padStringRightToLength(string.sub(body, iter, iter + width - 1), width) .. "|\n"
            else
                result = result .. "|" .. _padStringRightToLength(string.sub(body, iter, lastSpace - 1), width) .. "|\n"
                iter = lastSpace + 1
                if iter <= #body then
                    result = result .. "|" .. _padStringRightToLength(string.sub(body, iter, iter + width - 1), width) .. "|\n"
                end
            end
            break
        end
        if space - iter > width then
            result = result .. "|" .. _padStringRightToLength(string.sub(body, iter, lastSpace - 1), width) .. "|\n"
            iter = lastSpace + 1
        end
    end
    result = result .. "+"
    for _ = 1, width do
        result = result .. "-"
    end
    result = result .. "+\n"
    return result
end
        
return Tools
