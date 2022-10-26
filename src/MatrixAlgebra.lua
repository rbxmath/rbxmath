local MatrixAlgebra = {}

local _linearlyIndexedSparseMatrix = {}

local _linearlyIndexedSparseMatrixFromTableOfTables = function (tableOfTables)
    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _linearlyIndexedSparseMatrix)

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        for j = 1, numberOfColumns do
            if tableOfTables[i][j] ~= 0 then
                sparseForm[rowConstant + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end

local _linearlyIndexedSparseMatrixFromNumericTableOfTables = function (tableOfTables, tol)
    tol = tol or 0

    local numberOfRows = #tableOfTables
    local numberOfColumns = #tableOfTables[1]

    local sparseForm = setmetatable({}, _linearlyIndexedSparseMatrix)

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        for j = 1, numberOfColumns do
            if math.abs(tableOfTables[i][j]) >= tol then
                sparseForm[rowConstant + j] = tableOfTables[i][j]
            end
        end
    end

    rawset(sparseForm, "dimensions", {numberOfRows, numberOfColumns})

    return sparseForm
end

local _linearlyIndexedSparseCopy = function (matrix)
    local copy = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(copy, "dimensions", {matrix.dimensions[1], matrix.dimensions[2]})

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            copy[k] = v
        end
    end

    return copy
end

local _linearlyIndexedSparseIdentity = function (n)
    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(result, "dimensions", {n, n})

    for i = 1, n do
        result[n * (i - 1) + i] = 1
    end

    return result
end

local _linearlyIndexedSparseZero = function (n, m)
    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(result, "dimensions", {n,m})

    return result
end

local _linearlyIndexedSparseRandom = function (n, m, a, b, tol)
    tol = tol or 0

    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(result, "dimensions", {n, m})

    for i = 1, n do
        local rowConstant = m * (i - 1)
        for j = 1, m do
            local val = (b - a) * (math.random()) + a
            if math.abs(val) >= tol then
                result[rowConstant + j] = val
            end
        end
    end

    return result
end

local _linearlyIndexedSparseColumnSwap = function (matrix, n, m)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfRows do
        local rowConstant = numberOfColumns * (i - 1)
        matrix[rowConstant + n],  matrix[rowConstant + m] =  matrix[rowConstant + m],  matrix[rowConstant + n]
    end

    return matrix
end

local _linearlyIndexedSparseRightPermutationMatrix = function (n, permutations)
    local result = _linearlyIndexedSparseIdentity(n)

    for i = 1, n do
        local destination = permutations[i]
        if destination ~= nil then
            result = _linearlyIndexedSparseColumnSwap(result, i, destination)
        end
    end

    return result
end

local _linearlyIndexedSparseRowSwap = function (matrix, n, m)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfColumns do
        local rowConstant = numberOfColumns * (n - 1)
        local rowConstant2 = numberOfColumns * (m - 1)
        matrix[rowConstant + i],  matrix[rowConstant2 + i] =  matrix[rowConstant2 + i],  matrix[rowConstant + i]
    end

    return matrix
end

local _linearlyIndexedSparsePermutationMatrix = function (n, permutations)
    local result = _linearlyIndexedSparseIdentity(n)

    for i = 1, n do
        local destination = permutations[i]
        if destination ~= nil then
            result = _linearlyIndexedSparseRowSwap(result, i, destination)
        end
    end

    return result
end

local _linearlyIndexedSparseShear = function (matrix, n, m, c, tol)
    tol = tol or 0

    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfColumns do
        local val1 = matrix[numberOfColumns * (n - 1) + i]
        local val2 = matrix[numberOfColumns * (m - 1) + i]
        if val1 ~= nil and val2 ~= nil then
            local val3 = val2 + c * val1
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (m - 1) + i] = nil
            else
                matrix[numberOfColumns * (m - 1) + i] = val3
            end
        elseif val1 ~= nil then
            local val3 = c * val1
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (m - 1) + i] = nil
            else
                matrix[numberOfColumns * (m - 1) + i] = val3
            end
        end
    end

    return matrix
end

local _linearlyIndexedSparseRowAdd = function (matrix, row, n, c, tol)
    tol = tol or 0

    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfColumns do
        local val1 = matrix[numberOfColumns * (n - 1) + i]
        local val2 = row[i]
        if val1 ~= nil and val2 ~= nil then
            local val3 = val1 + c * val2
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (n - 1) + i] = nil
            else
                matrix[numberOfColumns * (n - 1) + i] = val3
            end
        elseif val2 ~= nil then
            local val3 = c * val2
            if math.abs(val3) < tol then
                matrix[numberOfColumns * (n - 1) + i] = nil
            else
                matrix[numberOfColumns * (n - 1) + i] = val3
            end
        end
    end

    return matrix
end

local _linearlyIndexedSparseTranspose = function (matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(result, "dimensions", {numberOfColumns, numberOfRows})

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            local columnNumber = k % numberOfColumns
            if columnNumber == 0 then
                columnNumber = numberOfColumns
            end
            local rowNumber = (k - columnNumber) / numberOfColumns + 1
            result[numberOfRows * (columnNumber - 1) + rowNumber] = v
        end
    end

    return result
end

local _linearlyIndexedSparseLU = function (matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute LU of sparse rectangular matrix.")
    end

    local permuations = {}

    local l = _linearlyIndexedSparseIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[numberOfColumns * (i - 1) + i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[numberOfColumns * (j - 1) + i] or 0
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate ~= nil and maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            _linearlyIndexedSparseRowSwap(matrix, i, maxRow)
            for k = 1, i - 1 do
                l[numberOfColumns * (i - 1) + k], l[numberOfColumns * (maxRow - 1) + k] = l[numberOfColumns * (maxRow - 1) + k], l[numberOfColumns * (i - 1) + k]
            end
            permuations[i] = maxRow
        end

        max = matrix[numberOfColumns * (i - 1) + i]

        for j = i + 1, numberOfRows do
            local val = matrix[numberOfColumns * (j - 1) + i]
            if val ~= nil then
                local valOverMax = val / max
                l[numberOfColumns * (j - 1) + i] = valOverMax
                matrix[numberOfColumns * (j - 1) + i] = nil
                for k = i + 1, numberOfColumns do
                    local val1 = matrix[numberOfColumns * (j - 1) + k]
                    local val2 = matrix[numberOfColumns * (i - 1) + k]
                    if val1 ~= nil and val2 ~= nil then
                        matrix[numberOfColumns * (j - 1) + k] = val1 - val2 * valOverMax
                    elseif val2 ~= nil then
                        matrix[numberOfColumns * (j - 1) + k] = -val2 * valOverMax
                    end
                end
            end
        end
    end

    local permuationMatrix = _linearlyIndexedSparsePermutationMatrix(numberOfRows, permuations)

    return {l, matrix, permuationMatrix}
end

local _linearlyIndexedSparseInverse = function (matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute inverse of sparse rectangular matrix.")
    end

    local result = _linearlyIndexedSparseIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[numberOfColumns * (i - 1) + i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[numberOfColumns * (j - 1) + i] or 0
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate ~= nil and maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            _linearlyIndexedSparseRowSwap(matrix, i, maxRow)
            _linearlyIndexedSparseRowSwap(result, i, maxRow)
        end

        max = matrix[numberOfColumns * (i - 1) + i]

        for j = i + 1, numberOfRows do
            local val = matrix[numberOfColumns * (j - 1) + i]
            if val ~= nil then
                local valOverMax = val / max
                matrix[numberOfColumns * (j - 1) + i] = nil
                _linearlyIndexedSparseShear(result, i, j, -valOverMax)
                for k = i + 1, numberOfColumns do
                    local val1 = matrix[numberOfColumns * (j - 1) + k]
                    local val2 = matrix[numberOfColumns * (i - 1) + k]
                    if val1 ~= nil and val2 ~= nil then
                        matrix[numberOfColumns * (j - 1) + k] = val1 - val2 * valOverMax
                    elseif val2 ~= nil then
                        matrix[numberOfColumns * (j - 1) + k] = -val2 * valOverMax
                    end
                end
            end
        end
    end

    for i = numberOfRows, 1, -1 do
        local rowConstant = numberOfColumns * (i - 1)
        local val = matrix[rowConstant + i]
        for j = 1, i - 1 do
            local val1 = matrix[numberOfColumns * (i - j - 1) + i]
            if val1 ~= nil then
                _linearlyIndexedSparseShear(result, i, i - j, -val1 / val)
            end
        end
        for j = 1, numberOfColumns do
            result[rowConstant + j] = result[rowConstant + j] / val
        end
    end

    return result
end

local _linearlyIndexedSparseMatrixSparsify = function (matrix, tol)
    tol = tol or 0

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            if v < tol then
                matrix[k] = nil
            end
        end
    end

    return matrix
end

local _sparseVector = {}

local _sparseVectorFromArray = function (array)
    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {})

    rawset(result, "length", 0)

    for i = 1, #array do
        if array[i] ~= 0 then
            result[i] = array[i]
            result.keys[result.keys.length + 1] = i
        end
    end

    return result
end

local _sparseVectorFromNumericArray = function (array, tol)
    tol = tol or 0

    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {})

    rawset(result, "length", #array)

    for i = 1, #array do
        if math.abs(array[i]) >= tol then
            result[i] = array[i]
            result.keys[#result.keys + 1] = i
        end
    end

    return result
end

local _sparseVectorCopy = function (vector)
    local copy = setmetatable({}, _sparseVector)

    rawset(copy, "keys", {})

    rawset(copy, "length", vector.length)

    for k in vector.keys do
        copy.keys[#copy.keys + 1] = k
    end

    return copy
end

local _sparseStandardBasisVector = function (n, m)
    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {m})

    rawset(result, "length", n)

    result[m] = 1

    return result
end

local _sparseZeroVector = function (n)
    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {})

    rawset(result, "length", n)

    return result
end

local _sparseRandomVector = function (n, a, b, tol)
    tol = tol or 0

    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {})

    rawset(result, "length", n)

    for i = 1, n do
        local val = (b - a) * (math.random()) + a
        if val >= tol then
            result[i] = val
            result.keys[#result.keys + 1] = i
        end
    end

    return result
end

local _sparseOnesVector = function (n)
    local result = setmetatable({}, _sparseVector)

    rawset(result, "keys", {})

    rawset(result, "length", n)

    for i = 1, n do
        result[i] = 1
        result.keys[#result.keys + 1] = i
    end

    return result
end

local _keyUnion = function (left, right)
    local j = 1

    if #left < #right then
        left, right = right, left
    end

    local keys = {}

    local index = 1

    for i = 1, #left do
        local val1 = left[i]
        local val2 = right[j]

        while val2 and val1 > val2 do
            keys[index] = val2
            index = index + 1
            j = j + 1
            val2 = right[j]
        end
        if val1 == val2 then
            keys[index] = val1
            index = index + 1
            j = j + 1
        else
            keys[index] = val1
            index = index + 1
        end
    end

    return keys
end

local _sparseVectorSwap = function (vector, n, m)
    vector[n], vector[m] = vector[m], vector[n]

    return vector
end

local _sparseMatrix = {}

local _sparseMatrixFromArrayOfArrays = function (ArrayOfArrays)
    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "dimensions", {#ArrayOfArrays, #ArrayOfArrays[1]})

    rawset(result, "keys", {})

    for i = 1, #ArrayOfArrays do
        local row = _sparseVectorFromArray(ArrayOfArrays[i])
        if #row.keys ~= 0 then
            result[i] = row
            result.keys[#result.keys+1] = i
        end
    end

    return result
end

local _sparseMatrixIdentity = function (n)
    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "dimensions", {n, n})

    rawset(result, "keys", {})

    for i = 1, n do
        result[i] = _sparseStandardBasisVector(n, i)
        result.keys[#result.keys+1] = i
    end

    return result
end

local _sparseMatrixRowPermute = function (matrix, permuations)
    for n, m in ipairs(permuations) do
        matrix[n], matrix[m] = matrix[m], matrix[n]
    end

    local keys = {}

    for k, v in ipairs(matrix) do
        keys[#keys+1] = k
    end

    matrix.keys = keys

    return matrix
end

local _sparseLeftPermutationMatrix = function (n, permuations)
    local result = _sparseMatrixIdentity(n)

    return _sparseMatrixRowPermute(result, permuations)
end

local _sparseMatrixLU = function (matrix)
    
end

local _sparseMatrixColumn = function (matrix, n)
    local result = setmetatable({}, _sparseVector)

    local keys = {}

    rawset(result, "length", matrix.dimensions[1])

    for k, i in ipairs(matrix.keys) do
        result[i] = matrix[i][n]
        keys[#keys+1] = i
    end

    rawset(result, "keys", keys)

    return result
end

local _sparseMatrixColumnKeys = function (matrix)
    local keys = {}

    for k, i in ipairs(matrix.keys) do
        keys = _keyUnion(keys, matrix[i].keys)
    end

    return keys
end

local _sparseMatrixFlatten = function (matrix)
    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    rawset(result, "dimensions", {numberOfRows, numberOfColumns})

    for k, i in ipairs(matrix.keys) do
        local row = matrix[i]
        local rowConstant = numberOfColumns * (i - 1)
        for kk, ii in ipairs(row.keys) do
            result[rowConstant + ii] = row[ii]
        end
    end

    return result
end

local _sparseMatrixUnFlatten = function (matrix)
    local result = setmetatable({}, _sparseMatrix)

    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    rawset(result, "dimensions", {numberOfRows, numberOfColumns})

    local keys = {}

    local rowNumber = 0

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            local columnNumber = k % numberOfColumns
            if columnNumber == 0 then
                columnNumber = numberOfColumns
            end
            local tempRowNumber = math.floor((k - columnNumber) / numberOfColumns + 1.1)
            if tempRowNumber > rowNumber then
                local row = _sparseZeroVector(numberOfColumns)
                row[columnNumber] = v
                row.keys[#row.keys+1] = columnNumber
                result[tempRowNumber] = row
                keys[#keys+1] = tempRowNumber
                rowNumber = tempRowNumber
            else
                local row = result[rowNumber]
                row[columnNumber] = v
                row.keys[#row.keys+1] = columnNumber
                keys[#keys+1] = rowNumber
            end
        end
    end

    rawset(result, "keys", keys)

    return result
end

local _matrix = {}

local _matrixFromTableOfTables = function (tableOfTables)
    local result = setmetatable(tableOfTables, _matrix)

    rawset(result, "dimensions", {#tableOfTables,#tableOfTables[1]})

    return result
end

local _nthStandardBasisVector = function (n, m)
    local e = {}

    for i = 1, m do
        if i == n then
            e[i] = 1
        else
            e[i] = 0
        end
    end

    return e
end

local _zeroVector = function (n)
    local z = {}

    for i = 1, n do
        z[i] = 0
    end

    return z
end

local _matrixIdentity = function (n)
    local result = {}

    for i = 1, n do
        result[i] = _nthStandardBasisVector(i, n)
    end

    return _matrixFromTableOfTables(result)
end

local _matrixZero = function (n, m)
    local result = {}

    for i = 1, n do
        result[i] = _zeroVector(m)
    end

    return _matrixFromTableOfTables(result)
end

_sparseVector.__add = function (left, right)
    if left.length ~= right.length then
        error("Cannot add sparse vectors of unequal length.", 2)
    end

    local keys = _keyUnion(left.keys, right.keys)

    local result = setmetatable({}, _sparseVector)

    rawset(result, "length", left.length)

    local newKeys = {}

    for k, i in ipairs(keys) do
        local val1 = left[i]
        local val2 = right[i]
        if val1 and val2 then
            local val3 = val1 + val2
            if val3 == 0 then
                result[i] = nil
            else
                result[i] = val3
                newKeys[#newKeys + 1] = i
            end
        elseif val2 then
            result[i] = val2
            newKeys[#newKeys + 1] = i
        else
            result[i] = val1
            newKeys[#newKeys + 1] = i
        end
    end

    result.keys = newKeys

    return result
end

_sparseVector.__sub = function (left, right)
    if left.length ~= right.length then
        error("Cannot add sparse vectors of unequal length.", 2)
    end

    local keys = _keyUnion(left.keys, right.keys)

    local result = setmetatable({}, _sparseVector)

    rawset(result, "length", left.length)

    local newKeys = {}

    for k, i in ipairs(keys) do
        local val1 = left[i]
        local val2 = right[i]
        if val1 and val2 then
            local val3 = val1 - val2
            if val3 == 0 then
                result[i] = nil
            else
                result[i] = val3
                newKeys[#newKeys + 1] = i
            end
        elseif val2 then
            result[i] = -val2
            newKeys[#newKeys + 1] = i
        else
            result[i] = val1
            newKeys[#newKeys + 1] = i
        end
    end

    result.keys = newKeys

    return result
end

_sparseVector.__unm = function (vector)
    local copy = _sparseVectorCopy(vector)

    for k, i in ipairs(copy.keys) do
        copy[i] = -copy[i]
    end

    return copy
end

_sparseVector.__mul = function (left, right)
    if left.length ~= right.length then
        error("Cannot take inner product of sparse vectors of unequal length.", 2)
    end

    local keys = _keyUnion(left.keys, right.keys)

    local result = 0

    for k, i in ipairs(keys) do
        local val1 = left[i]
        local val2 = right[i]
        if val1 and val2 then
            result = result + val1 * val2
        end
    end

    return result
end

_sparseVector.__tostring = function (vector)
    local result = "{"

    local length = vector.length

    for i = 1, length - 1 do
        local val = vector[i]
        if val then
            result = result .. val .. ","
        else
            result = result .. 0 .. ","
        end
    end

    local val = vector[length]
    if val then
        result = result .. val
    else
        result = result .. 0
    end

    result = result .. "}"

    return result
end

_sparseMatrix.__add = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] and left.dimensions[2] ~= right.dimensions[2] then
        error("Cannot add sparse matrices of unequal dimensions.", 2)
    end

    local keys = _keyUnion(left.keys, right.keys)

    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "dimensions", {left.dimensions[1], left.dimensions[2]})

    local newKeys = {}

    for k, i in ipairs(keys) do
        local val1 = left[i]
        local val2 = right[i]
        if val1 and val2 then
            local val3 = val1 + val2
            if val3 == 0 then
                result[i] = nil
            else
                result[i] = val3
                newKeys[#newKeys + 1] = i
            end
        elseif val2 then
            result[i] = val2
            newKeys[#newKeys + 1] = i
        else
            result[i] = val1
            newKeys[#newKeys + 1] = i
        end
    end

    result.keys = newKeys

    return result
end

_sparseMatrix.__sub = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] and left.dimensions[2] ~= right.dimensions[2] then
        error("Cannot add sparse matrices of unequal dimensions.", 2)
    end

    local keys = _keyUnion(left.keys, right.keys)

    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "dimensions", {left.dimensions[1], left.dimensions[2]})

    local newKeys = {}

    for k, i in ipairs(keys) do
        local val1 = left[i]
        local val2 = right[i]
        if val1 and val2 then
            local val3 = val1 - val2
            if val3 == 0 then
                result[i] = nil
            else
                result[i] = val3
                newKeys[#newKeys + 1] = i
            end
        elseif val2 then
            result[i] = -val2
            newKeys[#newKeys + 1] = i
        else
            result[i] = val1
            newKeys[#newKeys + 1] = i
        end
    end

    result.keys = newKeys

    return result
end

_sparseMatrix.__mul = function (left, right)
    local liLeft = _sparseMatrixFlatten(left)
    local liRight = _sparseMatrixFlatten(right)

    local result = _sparseMatrixUnFlatten(liLeft * liRight)

    return result
end

_sparseMatrix.__tostring = function (matrix)
    local result = "{"

    local length = matrix.dimensions[1]
    local width = matrix.dimensions[2]

    for i = 1, length - 1 do
        local val = matrix[i]
        if val then
            result = result .. tostring(val) .. ","
        else
            result = result .. tostring(_sparseZeroVector(width)) .. ","
        end
    end

    local val = matrix[length]
    if val then
        result = result .. tostring(val)
    else
        result = result .. tostring(_sparseZeroVector(width))
    end

    result = result .. "}"

    return result
end

_linearlyIndexedSparseMatrix.__add = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] or left.dimensions[2] ~= right.dimensions[2] then
        error("Attempting to add sparse matrices of different sizes.")
    end

    local inflatedLeft = _sparseMatrixUnFlatten(left)
    local inflatedRight = _sparseMatrixUnFlatten(right)

    return _sparseMatrixFlatten(inflatedLeft + inflatedRight)
end

_linearlyIndexedSparseMatrix.__sub = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] or left.dimensions[2] ~= right.dimensions[2] then
        error("Attempting to add sparse matrices of different sizes.")
    end

    local inflatedLeft = _sparseMatrixUnFlatten(left)
    local inflatedRight = _sparseMatrixUnFlatten(right)

    return _sparseMatrixFlatten(inflatedLeft - inflatedRight)
end

_linearlyIndexedSparseMatrix.__mul = function (left, right)
    if left.dimensions[2] ~= right.dimensions[1] then
        error("Attempting to multiply incompatible sparse matrices.")
    end

    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    rawset(result, "dimensions", {left.dimensions[1], right.dimensions[2]})

    for k, v in pairs(left) do
        if type(k) == "number" then
            local columnNumber = k % left.dimensions[2]
            if columnNumber == 0 then
                columnNumber = left.dimensions[2]
            end
            local rowNumber = (k - columnNumber) / left.dimensions[2]
            local rowConstant = right.dimensions[2] * (columnNumber - 1)
            local rowConstant2 = right.dimensions[2] * rowNumber
            for i = 1, right.dimensions[2] do
                local val1 = result[rowConstant2 + i]
                local val2 = right[rowConstant + i]
                if type(val1) == "number" and val2 ~= nil then
                    val1 = val1 + v * right[rowConstant + i]
                    result[rowConstant2 + i] = val1
                elseif val2 ~= nil then
                    result[rowConstant2 + i] = v * right[rowConstant + i]
                end
            end
        end
    end

    return result
end

_linearlyIndexedSparseMatrix.__tostring = function (matrix)
    local result = "{"

    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfRows - 1 do
        result = result .. "{"
        for j = 1, numberOfColumns - 1 do
            local val = matrix[numberOfColumns * (i - 1) + j] or 0
            result = result .. tostring(val) .. ","
        end
        local val = matrix[numberOfColumns * (i - 1) + numberOfColumns] or 0
        result = result .. tostring(val) .. "},"
    end

    result = result .. "{"
    for j = 1, numberOfColumns - 1 do
        local val = matrix[numberOfColumns * (numberOfRows - 1) + j] or 0
        result = result .. tostring(val) .. ","
    end
    local val = matrix[numberOfColumns * (numberOfRows - 1) + numberOfColumns] or 0
    result = result .. tostring(val) .. "}"

    result = result .. "}"

    return result
end

_matrix.__add = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] or left.dimensions[2] ~= right.dimensions[2] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for k, v in ipairs(left) do
        result[k] = {}
        for kk, vv in ipairs(v) do
            result[k] = vv + right[k][kk]
        end
    end

    return _matrixFromTableOfTables(result)
end

_matrix.__sub = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] or left.dimensions[2] ~= right.dimensions[2] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for k, v in ipairs(left) do
        result[k] = {}
        for kk, vv in ipairs(v) do
            result[k] = vv - right[k][kk]
        end
    end

    return _matrixFromTableOfTables(result)
end

MatrixAlgebra.sparseVector = {}

MatrixAlgebra.sparseVector.new = function (array)
    return _sparseVectorFromArray(array)
end

MatrixAlgebra.sparseVector.newNumeric = function (array, tol)
    return _sparseVectorFromNumericArray(array, tol)
end

MatrixAlgebra.sparseVector.ones = function (n)
    return _sparseOnesVector(n)
end

MatrixAlgebra.sparseVector.zero = function (n)
    return _sparseZeroVector(n)
end

MatrixAlgebra.sparseVector.random = function (n, a, b, tol)
    return _sparseRandomVector(n, a, b, tol)
end

MatrixAlgebra.sparseMatrix = {}

MatrixAlgebra.sparseMatrix.new = function (ArrayOfArrays)
    return _sparseMatrixFromArrayOfArrays(ArrayOfArrays)
end

MatrixAlgebra.sparseMatrix.identity = function (n)
    return _sparseMatrixIdentity(n)
end

MatrixAlgebra.sparseMatrix.flatten = function (matrix)
    return _sparseMatrixFlatten(matrix)
end

MatrixAlgebra.liSparseMatrix = {}

MatrixAlgebra.liSparseMatrix.new = function (tableOfTables)
    return _linearlyIndexedSparseMatrixFromTableOfTables(tableOfTables)
end

MatrixAlgebra.liSparseMatrix.newNumeric = function (tableOfTables, tol)
    return _linearlyIndexedSparseMatrixFromNumericTableOfTables(tableOfTables, tol)
end

MatrixAlgebra.liSparseMatrix.identity = function (n)
    return _linearlyIndexedSparseIdentity(n)
end

MatrixAlgebra.liSparseMatrix.zero = function (n, m)
    return _linearlyIndexedSparseZero(n, m)
end

MatrixAlgebra.liSparseMatrix.random = function (n, m, a, b, tol)
    return _linearlyIndexedSparseRandom(n, m, a, b, tol)
end

MatrixAlgebra.liSparseMatrix.lu = function (input)
    local matrix = _linearlyIndexedSparseCopy(input)
    return _linearlyIndexedSparseLU(matrix)
end

MatrixAlgebra.liSparseMatrix.inverse = function (input)
    local matrix = _linearlyIndexedSparseCopy(input)
    return _linearlyIndexedSparseInverse(matrix)
end

MatrixAlgebra.liSparseMatrix.transpose = function (input)
    local matrix = _linearlyIndexedSparseCopy(input)
    return _linearlyIndexedSparseTranspose(matrix)
end

MatrixAlgebra.liSparseMatrix.sparsify = function (input, tol)
    local matrix = _linearlyIndexedSparseCopy(input)
    return _linearlyIndexedSparseMatrixSparsify(matrix, tol)
end

MatrixAlgebra.liSparseMatrix.unflatten = function (matrix)
    return _sparseMatrixUnFlatten(matrix)
end

MatrixAlgebra.matrix = {}

MatrixAlgebra.matrix.new = function (tableOfTables)
    return _matrixFromTableOfTables(tableOfTables)
end

MatrixAlgebra.matrix.identity = function (n)
    return _matrixIdentity(n)
end

MatrixAlgebra.matrix.zero = function (n, m)
    return _matrixZero(n, m)
end

return MatrixAlgebra
