local shallowCopyArray = function (array)
    local result = {}

    for i = 1, #array do
        result[i] = array[i]
    end

    return result
end

local removeDuplicatesFromSortedList = function (array)
    local result = {array[1]}

    for i = 2, #array do
        if array[i - 1] ~= array[i] then
            result[#result+1] = array[i]
        end
    end

    return result
end

local arrayEquals = function (left, right)
    if #left ~= #right then
        return false
    end

    for i = 1, #left do
        if left[i] ~= right[i] then
            return false
        end
    end

    return true
end

local hadamardProduct = function (left, right)
    if #left ~= #right then
        error("Array length mismatch for Hadamard product!", 2)
    end

    local result = {}

    for i = 1, #left do
        result[i] = left[i] * right[i]
    end

    return result
end

local lexicographicCompare = function (left, right)
    local i = 1

    while i < #left and i < #right do
        if left[i] == right[i] then
            i = i + 1
        else
            return left[i] < right[i]
        end
    end

    if #left >= #right then
        return left[i] < right[i]
    else
        return left[i] <= right[i]
    end
end

local reverseArray = function (array)
    local result = {}

    for i = 1, #array do
        result[#array - i + 1] = array[i]
    end

    return result
end

local concatenateArrays = function (left, right)
    local result = {}

    for i = 1, #left do
        result[i] = left[i]
    end

    for i = 1, #right do
        result[#result+1] = right[i]
    end

    return result
end

local subarray = function (array, index)
    local result = {}

    for i = 1, math.min(index, #array) do
        result[i] = array[i]
    end

    return result
end

local suparray = function (array, index)
    local result = {}

    for i = index, #array do
        result[#result+1] = array[i]
    end
end

local arrayToString = function (array)
    local result = "("

    for i = 1, #array - 1 do
        result = result .. tostring(array[i]) .. " "
    end

    result = result .. tostring(array[#array]) .. ")"

    return result
end

local binarySearch
binarySearch = function (element, array)
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
        return binarySearch(element, subarray(array, checkIndex - 1))
    elseif element > array[checkIndex] then
        return checkIndex + binarySearch(element, suparray(array, checkIndex + 1))
    else
        return checkIndex
    end
end

local integerPower = function (base, exponent)
    local result = base

    for i = 2, exponent do
        result = result * base
    end

    return result
end

local _set = {}

local _setFromArray = function (array)
    local result = shallowCopyArray(array)
    table.sort(result)
    result = removeDuplicatesFromSortedList(result)

    setmetatable(result, _set)

    return result
end

_set.__add = function (left, right)
    local result = _setFromArray({})

    local i = 0
    local j = 0
    local index = 1

    while i < #left or j < #right do
        if i == #left then
            result[index] = right[j + 1]
            j = j + 1
            index = index + 1
        elseif j == #right then
            result[index] = left[i + 1]
            i = i + 1
            index = index + 1
        elseif left[i + 1] < right[j + 1] then
            result[index] = left[i + 1]
            i = i + 1
            index = index + 1
        elseif left[i + 1] == right[j + 1] then
            result[index] = left[i + 1]
            i = i + 1
            j = j + 1
            index = index + 1
        else
            result[index] = right[j + 1]
            j = j + 1
            index = index + 1
        end
    end

    return result
end

_set.__mul = function (left, right)
    local result = _setFromArray({})

    if left == nil or right == nil or #left == 0 or #right == 0 then
        return result
    end

    local i = 1
    local j = 1

    while i <= #left and j <= #right do
        while i <= #left and j <= #right and left[i] < right[j] do
            i = i + 1
        end
        while i <= #left and j <= #right and left[i] > right[j] do
            j = j + 1
        end
        if left[i] == right[j] then
            result[#result+1] = left[i]
            i = i + 1
            j = j + 1
        end
    end

    return result
end

_set.__eq = function (left, right)
    if #left ~= #right then
        return false
    end

    for i = 1, #left do
        if left[i] ~= right[i] then
            return false
        end
    end

    return true
end

_set.__lt = function (left, right)
    if #left >= #right then
        return false
    end

    return left * right == left
end

_set.__le = function (left, right)
    return left * right == left
end

local _duple = {}

local _dupleFromTwoInputs = function (left, right)
    local result = {left, right}

    setmetatable(result, _duple)

    return result
end

_duple.__eq = function (left, right)
    return left[1] == right[1] and left[2] == right[2]
end

_duple.__lt = function (left, right)
    if left[1] == right[1] then
        return left[2] < right[2]
    else
        return left[1] < right[1]
    end
end

_duple.__le = function (left, right)
    return left < right or left == right
end

local directSum = function (left, right)
    if #left ~= #right then
        error("Cardinality mismatch between left and right arrays!", 2)
    end

    local result = {}

    for i = 1, #left do
        result[i] = _dupleFromTwoInputs(left[i], right[i])
    end

    return result
end

local _monomial = {}

local _monomialFromArrayAndPreSet = function (preset, array, coeff)
    if #array ~= #preset then
        error("Cardinality mismatch between power array and symbol set!", 2)
    end

    if preset == nil or #preset == 0 then
        local result = {_setFromArray({}),{}}
        setmetatable(result, _monomial)
        rawset(result, "coeff", coeff)
        return result
    end

    local preresult = directSum(preset, array)
    table.sort(preresult)
    local arrayResult = {preresult[1][2]}
    local presetResult = {preresult[1][1]}

    for i = 2, #preresult do
        if preresult[i - 1][1] ~= preresult[i][1] then
            arrayResult[#arrayResult+1] = preresult[i][2]
            presetResult[#presetResult+1] = preresult[i][1]
        else
            arrayResult[#arrayResult] = arrayResult[#arrayResult] + preresult[i][2]
        end
    end

    local result = {_setFromArray(presetResult), arrayResult}
    setmetatable(result, _monomial)
    rawset(result, "coeff", coeff)

    return result
end

local _monomialScale = function (c, monomial)
    return _monomialFromArrayAndPreSet(shallowCopyArray(monomial[1]), shallowCopyArray(monomial[2]), c * monomial["coeff"])
end

local _monomialDegree = function (monomial)
    local result = 0

    for i = 1, #monomial do
        result = result + monomial[2][i]
    end

    return result
end

local _monomialReRep = function (left, right)
    local symbolSet = left[1] * right[1]

    local zeroArray = {}

    for i = 1, #symbolSet do
        zeroArray[i] = 0
    end

    local leftPreSet = concatenateArrays(left[1], symbolSet)
    local leftArray = concatenateArrays(left[2], zeroArray)
    local rightPreSet = concatenateArrays(right[1], symbolSet)
    local rightArray = concatenateArrays(right[2], zeroArray)

    return _monomialFromArrayAndPreSet(leftPreSet, leftArray, left["coeff"]), _monomialFromArrayAndPreSet(rightPreSet, rightArray, right["coeff"])
end

local _monomialCopy = function (monomial)
    local result = _monomialFromArrayAndPreSet(shallowCopyArray(monomial[1]), shallowCopyArray(monomial[2]), monomial["coeff"])

    return result
end

local _monomialReplace = function (monomial, symbolToReplace, symbol)
    local index = binarySearch(symbolToReplace, monomial[1])

    if index == false then
        return _monomialCopy(monomial)
    else
        local result = shallowCopyArray(monomial[1])
        result[index] = symbol
        return _monomialFromArrayAndPreSet(result, shallowCopyArray(monomial[2]), monomial["coeff"])
    end
end

local _monomialEvaluate = function (monomial, rules)
    local result = _monomialCopy(monomial)
    local constant = 1

    if monomial[1] == nil or monomial[2] == nil then
        return monomial["coeff"]
    end

    for key, value in pairs(rules) do
        local index = binarySearch(key, monomial[1])
        if index ~= false then
            constant = constant * integerPower(value, monomial[2][index])
            result[2][index] = 0
        end
    end

    result["coeff"] = constant * result["coeff"]

    return result
end

local _monomialSimplify = function (monomial)
    local exponentArray = {}
    local symbolPreSet = {}

    if monomial["coeff"] == 0 then
        return _monomialFromArrayAndPreSet({},{},0)
    end

    for i = 1, #monomial[2] do
        if monomial[2][i] ~= 0 then
            exponentArray[#exponentArray+1] = monomial[2][i]
            symbolPreSet[#symbolPreSet+1] = monomial[1][i]
        end
    end

    return _monomialFromArrayAndPreSet(symbolPreSet, exponentArray, monomial["coeff"])
end

_monomial.__add = function (left, right)
    local leftSimple = _monomialSimplify(left)
    local rightSimple = _monomialSimplify(right)

    if leftSimple[1] ~= rightSimple[1] then
        error("Monomial symbol set mismatch!", 2)
    elseif not arrayEquals(leftSimple[2], rightSimple[2]) then
        error("Monomial degree mismatch!", 2)
    end

    return _monomialFromArrayAndPreSet(left[1], left[2], left["coeff"] + right["coeff"])
end

_monomial.__mul = function (left, right)
    return _monomialFromArrayAndPreSet(concatenateArrays(left[1], right[1]), concatenateArrays(left[2], right[2]), left["coeff"] * right["coeff"])
end

_monomial.__unm = function (monomial)
    return _monomialFromArrayAndPreSet(shallowCopyArray(monomial[1]), shallowCopyArray(monomial[2]), -monomial["coeff"])
end

_monomial.__sub = function (left, right)
    return left + -right
end

_monomial.__len = function (monomial)
    return #monomial[1]
end

_monomial.__tostring = function (monomial)
    local result = tostring(monomial["coeff"])

    for i = 1, #monomial do
        result = result .. tostring(monomial[1][i]) .. "^" .. tostring(monomial[2][i])
    end

    return result
end

_monomial.__eq = function (left, right)
    local leftSimple = _monomialSimplify(left)
    local rightSimple = _monomialSimplify(right)

    return leftSimple[1] == rightSimple[1] and arrayEquals(leftSimple[2], rightSimple[2])
end

_monomial.__lt = function (left, right)
    if _monomialDegree(left) ~= _monomialDegree(right) then
        return _monomialDegree(left) < _monomialDegree(right)
    else
        if _monomialDegree(left) == 0 then
            return true
        end
        local leftCopy, rightCopy = _monomialReRep(left, right)
        return lexicographicCompare(reverseArray(leftCopy[2]), reverseArray(rightCopy[2]))
    end
end

_monomial.__le = function (left, right)
    return left == right or left < right
end

local _polynomial = {}

local _polynomialFromArrayOfMonomials = function (array)
    local preresult = {}
    for i = 1, #array do
        preresult[i] = array[i]
    end
    table.sort(preresult)
    local result = {preresult[1]}

    for i = 2, #array do
        if preresult[i - 1] == preresult[i] then
            result[#result] = result[#result] + preresult[i]
        else
            result[#result+1] = preresult[i]
        end
    end

    setmetatable(result, _polynomial)
    return result
end

local _polynomialCopy = function (polynomial)
    local result = {}

    for i = 1, #polynomial do
        result[#result+1] = _monomialCopy(polynomial[i])
    end
end

local _polynomialReplace = function (polynomial, symbolToReplace, symbol)
    local result = {}

    for i = 1, #polynomial do
        result[i] = _monomialReplace(polynomial[i], symbolToReplace, symbol)
    end

    return _polynomialFromArrayOfMonomials(result)
end

local _polynomialScale = function (c, polynomial)
    local result = {}

    for i = 1, #polynomial do
        result[i] = _monomialScale(c, polynomial[i])
    end

    return _polynomialFromArrayOfMonomials(result)
end

local _polynomialSimplify = function (polynomial)
    local result = {}

    for i = 1, #polynomial do
        local simple = _monomialSimplify(polynomial[i])
        if simple["coeff"] ~= 0 then
            result[#result+1] = simple
        end
    end

    return _polynomialFromArrayOfMonomials(result)
end

local _polynomialFromArray = function (array)
    local result = {}

    for i = 1, #array do
        result[i] = _monomialFromArrayAndPreSet({"x"}, {i - 1}, array[i])
    end

    return _polynomialFromArrayOfMonomials(result)
end

local _polynomialEvaluate = function (polynomial, rules)
    local result = {}

    for i = 1, #polynomial do
        local simple = _monomialSimplify(_monomialEvaluate(polynomial[i], rules))
        if simple["coeff"] ~= 0 then
            result[#result+1] = simple
        end
    end

    return _polynomialFromArrayOfMonomials(result)
end

_polynomial.__add = function (left, right)
    local preresult = concatenateArrays(left, right)

    return _polynomialFromArrayOfMonomials(preresult)
end

_polynomial.__unm = function (polynomial)
    local result = {}

    for i = 1, #polynomial do
        result[i] = -polynomial[i]
    end

    return _polynomialFromArrayOfMonomials(result)
end

_polynomial.__sub = function (left, right)
    return left + -right
end

_polynomial.__mul = function (left, right)
    local result = {}

    for i = 1, #left do
        for j = 1, #right do
            result[#result+1] = left[i] * right[j]
        end
    end

    return _polynomialFromArrayOfMonomials(result)
end

_polynomial.__tostring = function (polynomial)
    local result = tostring(polynomial[1])

    for i = 2, #polynomial do
        result = result .. " + " .. tostring(polynomial[i])
    end

    return result
end

local SymbolicAlgebra = {}

SymbolicAlgebra.set = {}

SymbolicAlgebra.set.new = function (array)
    return _setFromArray(array)
end

SymbolicAlgebra.duple = {}

SymbolicAlgebra.duple.new = function (first, second)
    return _dupleFromTwoInputs(first, second)
end

SymbolicAlgebra.monomial = {}

SymbolicAlgebra.monomial.new = function (exponents, symbols, coeff)
    return _monomialFromArrayAndPreSet(exponents, symbols, coeff)
end

SymbolicAlgebra.monomial.copy = function (monomial)
    return _monomialCopy(monomial)
end

SymbolicAlgebra.monomial.replace = function (monomial, symbolToReplace, symbol)
    return _monomialReplace(monomial, symbolToReplace, symbol)
end

SymbolicAlgebra.monomial.scale = function (c, monomial)
    return _monomialScale(c, monomial)
end

SymbolicAlgebra.monomial.simplify = function (monomial)
    return _monomialSimplify(monomial)
end

SymbolicAlgebra.monomial.evaluate = function (monomial, rules)
    return _monomialEvaluate(monomial, rules)
end

SymbolicAlgebra.polynomial = {}

SymbolicAlgebra.polynomial.new = function (array)
    if getmetatable(array[1]) == _monomial then
        return _polynomialSimplify(_polynomialFromArrayOfMonomials(array))
    else
        return _polynomialSimplify(_polynomialFromArray(array))
    end
end

SymbolicAlgebra.polynomial.copy = function (polynomial)
    return _polynomialCopy(polynomial)
end

SymbolicAlgebra.polynomial.replace = function (polynomial, symbolToReplace, symbol)
    return _polynomialReplace(polynomial, symbolToReplace, symbol)
end

SymbolicAlgebra.polynomial.scale = function (c, polynomial)
    return _polynomialScale(c, polynomial)
end

SymbolicAlgebra.polynomial.simplify = function (polynomial)
    return _polynomialSimplify(polynomial)
end

SymbolicAlgebra.polynomial.evaluate = function (polynomial, rules)
    return _polynomialEvaluate(polynomial, rules)
end

return SymbolicAlgebra