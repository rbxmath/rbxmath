local Tools = require("src.Tools")

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

local _linearlyIndexedSparseMatrixSquareSolve
_linearlyIndexedSparseMatrixSquareSolve = function (matrix, vector)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot solve sparse rectangular system with this function.")
    end

    local columnVector = _linearlyIndexedSparseTranspose(_linearlyIndexedSparseMatrixFromTableOfTables({vector}))

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
            error("Sparse matrix system is not solvable")
        end

        if maxRow ~= i then
            _linearlyIndexedSparseRowSwap(matrix, i, maxRow)
            _linearlyIndexedSparseRowSwap(columnVector, i, maxRow)
        end

        max = matrix[numberOfColumns * (i - 1) + i]

        for j = i + 1, numberOfRows do
            local val = matrix[numberOfColumns * (j - 1) + i]
            if val ~= nil then
                local valOverMax = val / max
                local columnVal1, columnVal2 = columnVector[j], columnVector[i]
                if columnVal1 ~= nil and columnVal2 ~= nil then
                    columnVector[j] = columnVal1 - valOverMax * columnVal2
                elseif columnVal2 ~= nil then
                    columnVector[j] = -valOverMax * columnVal2
                end
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

    for i = numberOfRows, 1, -1 do
        local temp = 0
        for j = i+1, numberOfColumns, 1 do
            temp = temp + matrix[numberOfColumns * (i - 1) + j] * columnVector[j]
        end
        columnVector[i] = columnVector[i] or 0
        columnVector[i] = (columnVector[i] - temp) / matrix[numberOfColumns * (i - 1) + i]
    end

    return columnVector
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

    rawset(result, "length", #array)

    for i = 1, #array do
        if array[i] ~= 0 then
            result[i] = array[i]
            result.keys[#result.keys + 1] = i
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

    for _, k in ipairs(vector.keys) do
        if vector[k] ~= nil then
            copy.keys[#copy.keys + 1] = k
            copy[k] = vector[k]
        end
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
    table.sort(left)
    table.sort(right)

    local result = {}

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

local _sparseVectorSwap = function (vector, n, m)
    vector[n], vector[m] = vector[m], vector[n]

    return vector
end

local _sparseVectorScale = function (vector, scale)
    if type(scale) == "table" then
        for i, v in ipairs(vector.keys) do
            vector[v] = scale[v] * vector[v]
        end
    elseif type(scale) == "number" then
        for i, v in ipairs(vector.keys) do
            vector[v] = scale * vector[v]
        end
    end

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

local _sparseMatrixFromNumericArrayOfArrays = function (ArrayOfArrays, tol)
    tol = tol or 0
    local result = setmetatable({}, _sparseMatrix)

    rawset(result, "dimensions", {#ArrayOfArrays, #ArrayOfArrays[1]})

    rawset(result, "keys", {})

    for i = 1, #ArrayOfArrays do
        local row = _sparseVectorFromNumericArray(ArrayOfArrays[i], tol)
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

local _sparseMatrixTranspose = function (matrix)
    local result = setmetatable({}, _sparseMatrix)

    result.dimensions = {}
    result.dimensions[1], result.dimensions[2] = matrix.dimensions[2], matrix.dimensions[1]

    result.keys = {}

    for i = 1, matrix.dimensions[2] do
        local array = {}
        for j = 1, matrix.dimensions[1] do
            array[j] = matrix[j][i] or 0
        end
        local row = _sparseVectorFromArray(array)
        if #row.keys ~= 0 then
            result[i] = row
            result.keys[#result.keys+1] = i
        end
    end

    return result
end

local _sparseMatrixLU = function (matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute LU of sparse rectangular matrix.")
    end

    local permuations = {}

    local l = _sparseMatrixIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i] or 0
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
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            for k = 1, i - 1 do
                l[i][k], l[maxRow][k] = l[maxRow][k], l[i][k]
            end
            permuations[i] = maxRow
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            if val ~= nil then
                local valOverMax = val / max
                l[j][i] = valOverMax
                matrix[j][i] = nil
                for k = i + 1, numberOfColumns do
                    local val1 = matrix[j][k]
                    local val2 = matrix[i][k]
                    if val1 ~= nil and val2 ~= nil then
                        matrix[j][k] = val1 - val2 * valOverMax
                    elseif val2 ~= nil then
                        matrix[j][k] = -val2 * valOverMax
                    end
                end
            end
        end
    end

    local permuationMatrix = _sparseLeftPermutationMatrix(numberOfRows, permuations)

    return {l, matrix, permuationMatrix}
end

local _sparseMatrixInverse = function (matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute inverse of sparse rectangular matrix.")
    end

    local result = _sparseMatrixIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i] or 0
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
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            result[i], result[maxRow] = result[maxRow], result[i]
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            if val ~= nil then
                local valOverMax = val / max
                matrix[j][i] = nil
                result[j] = result[j] + _sparseVectorScale(_sparseVectorCopy(result[i]), -valOverMax)
                for k = i + 1, numberOfColumns do
                    local val1 = matrix[j][k]
                    local val2 = matrix[i][k]
                    if val1 ~= nil and val2 ~= nil then
                        matrix[j][k] = val1 - val2 * valOverMax
                    elseif val2 ~= nil then
                        matrix[j][k] = -val2 * valOverMax
                    end
                end
            end
        end
    end

    for i = numberOfRows, 1, -1 do
        local val = matrix[i][i]
        for j = 1, i - 1 do
            local val1 = matrix[i - j][i]
            if val1 ~= nil then
                result[i - j] = result[i - j] + _sparseVectorScale(_sparseVectorCopy(result[i]), -val1 / val)
            end
        end
        for j = 1, numberOfColumns do
            if result[i][j] then
                result[i][j] = result[i][j] / val
            end
        end
    end

    return result
end

local _sparseMatrixSquareSolve
_sparseMatrixSquareSolve = function (matrix, vector)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    if numberOfRows ~= numberOfColumns then
        error("Cannot solve sparse rectangular system with this function.")
    end

    if #matrix.keys ~= numberOfRows then
        error("Cannot solve degenerate system with this function.")
    end

    local columnVector = _sparseMatrixTranspose(_sparseMatrixFromArrayOfArrays({vector}))

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i] or 0
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate ~= nil and maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix system is not solvable")
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            columnVector[i], columnVector[maxRow] = columnVector[maxRow], columnVector[i]
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            if val ~= nil then
                local valOverMax = val / max
                local columnVal1, columnVal2 = columnVector[j][1], columnVector[i][1]
                if columnVal1 ~= nil and columnVal2 ~= nil then
                    columnVector[j][1] = columnVal1 - valOverMax * columnVal2
                elseif columnVal2 ~= nil then
                    columnVector[j][1] = -valOverMax * columnVal2
                end
                matrix[j][i] = nil
                for k = i + 1, numberOfColumns do
                    local val1 = matrix[j][k]
                    local val2 = matrix[i][k]
                    if val1 ~= nil and val2 ~= nil then
                        matrix[j][k] = val1 - val2 * valOverMax
                    elseif val2 ~= nil then
                        matrix[j][k] = -val2 * valOverMax
                    end
                end
            end
        end
    end

    for i = numberOfRows, 1, -1 do
        local temp = 0
        for j = i+1, numberOfColumns, 1 do
            temp = temp + matrix[i][j] * columnVector[j][1]
        end
        columnVector[i][1] = columnVector[i][1] or 0
        if matrix[i][i] == nil then
            error("Sparse matrix system is not solvable")
        end
        columnVector[i][1] = (columnVector[i][1] - temp) / matrix[i][i]
    end

    return columnVector
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

local _sparseMatrixCopy = function (matrix)
    local result = {}
    setmetatable(result, _sparseMatrix)
    result.keys = Tools.list.copy(matrix.keys)
    result.dimensions = Tools.list.copy(matrix.dimensions)

    for _, v in ipairs(matrix.keys) do
        result[v] = _sparseVectorCopy(matrix[v])
    end

    return result
end

local _matrix = {}

local _matrixFromTableOfTables = function (tableOfTables)
    local result = {}

    for i = 1, #tableOfTables do
        result[i] = {}
        for j = 1, #tableOfTables[1] do
            result[i][j] = tableOfTables[i][j]
        end
    end

    result = setmetatable(result, _matrix)

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

local _matrixRowPermute = function (matrix, permuations)
    for n, m in ipairs(permuations) do
        matrix[n], matrix[m] = matrix[m], matrix[n]
    end
    return matrix
end

local _leftPermutationMatrix = function (n, permuations)
    local result = _matrixIdentity(n)

    return _matrixRowPermute(result, permuations)
end

local _matrixTranspose = function (matrix)
    local arrayOfArrays = _matrixZero(#matrix[1], #matrix)

    for i = 1, #matrix do
        for j = 1, #matrix[1] do
            arrayOfArrays[j][i] = matrix[i][j]
        end
    end

    return _matrixFromTableOfTables(arrayOfArrays)
end

local _matrixLU = function (matrix)
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute LU of rectangular matrix.")
    end

    local permuations = {}

    local l = _matrixIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i]

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i]
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            for k = 1, i - 1 do
                l[i][k], l[maxRow][k] = l[maxRow][k], l[i][k]
            end
            permuations[i] = maxRow
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            local valOverMax = val / max
            l[j][i] = valOverMax
            matrix[j][i] = 0
            for k = i + 1, numberOfColumns do
                matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
            end
        end
    end

    local permuationMatrix = _leftPermutationMatrix(numberOfRows, permuations)

    return {l, matrix, permuationMatrix}
end

local _matrixInverse = function (matrix)
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute inverse of sparse rectangular matrix.")
    end

    local result = _matrixIdentity(numberOfRows)

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i]
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix is not invertible")
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            result[i], result[maxRow] = result[maxRow], result[i]
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            local valOverMax = val / max
            matrix[j][i] = 0
            for k = 1, numberOfColumns do
                if k > i then
                    matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
                end
                result[j][k] = result[j][k] - result[i][k] * valOverMax
            end
        end
    end

    for i = numberOfRows, 1, -1 do
        local val = matrix[i][i]
        for j = 1, i - 1 do
            local val1 = matrix[i - j][i]
            for k = 1, numberOfColumns do
                result[i - j][k] = result[i - j][k] - val1 * result[i][k] / val
            end
        end
        for j = 1, numberOfColumns do
            result[i][j] = result[i][j] / val
        end
    end

    return result
end

local _matrixSquareSolve
_matrixSquareSolve = function (matrix, vector)
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    if numberOfRows ~= numberOfColumns then
        error("Cannot solve sparse rectangular system with this function.")
    end

    local columnVector = _matrixTranspose(_matrixFromTableOfTables({vector}))

    for i = 1, numberOfColumns - 1 do
        local maxRow = i
        local max = matrix[i][i]

        for j = i, numberOfRows do
            local maxCandidate = matrix[j][i]
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Sparse matrix system is not solvable")
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            columnVector[i], columnVector[maxRow] = columnVector[maxRow], columnVector[i]
        end

        max = matrix[i][i]

        for j = i + 1, numberOfRows do
            local val = matrix[j][i]
            local valOverMax = val / max
            local columnVal1, columnVal2 = columnVector[j][1], columnVector[i][1]
            columnVector[j][1] = columnVal1 - valOverMax * columnVal2
            matrix[j][i] = 0
            for k = i + 1, numberOfColumns do
                matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
            end
        end
    end

    local result = {}

    for i = numberOfRows, 1, -1 do
        local temp = 0
        for j = i+1, numberOfColumns, 1 do
            temp = temp + matrix[i][j] * columnVector[j][1]
        end
        columnVector[i][1] = columnVector[i][1]
        if matrix[i][i] == 0 then
            error("Sparse matrix system is not solvable")
        end
        columnVector[i][1] = (columnVector[i][1] - temp) / matrix[i][i]
    end

    for i = 1, numberOfRows do
        result[i] = columnVector[i][1]
    end

    return result
end

local _matrixFlatten = function (matrix)
    local result = setmetatable({}, _linearlyIndexedSparseMatrix)

    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    rawset(result, "dimensions", {numberOfRows, numberOfColumns})

    for i = 1, numberOfRows do
        local row = matrix[i]
        local rowConstant = numberOfColumns * (i - 1)
        for ii = 1, numberOfColumns do
            result[rowConstant + ii] = row[ii]
        end
    end

    return result
end

local _matrixUnFlatten = function (matrix)
    local result = setmetatable({}, _matrix)

    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    local rowNumber = 0

    for k, v in pairs(matrix) do
        if type(k) == "number" then
            local columnNumber = k % numberOfColumns
            if columnNumber == 0 then
                columnNumber = numberOfColumns
            end
            local tempRowNumber = math.floor((k - columnNumber) / numberOfColumns + 1.1)
            if tempRowNumber > rowNumber then
                local row = _zeroVector(numberOfColumns)
                row[columnNumber] = v
                result[tempRowNumber] = row
                rowNumber = tempRowNumber
            else
                local row = result[rowNumber]
                row[columnNumber] = v
            end
        end
    end

    return result
end

local _matrixCopy = function (matrix)
    local result = {}

    for i = 1, #matrix do
        result[i] = {}
        for j = 1, #matrix[1] do
            result[i][j] = matrix[i][j]
        end
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
    if left.dimensions[2] ~= right.dimensions[1] then
        error("Attempting to multiply incompatible sparse matrices.")
    end

    local rightTranspose = _sparseMatrixTranspose(right)

    local result = setmetatable({}, _sparseMatrix)
    result.keys, result.dimensions = {}, {}
    result.dimensions[1], result.dimensions[2] = left.dimensions[1], right.dimensions[2]

    for _, i in ipairs(left.keys) do
        local row = _sparseZeroVector(result.dimensions[2])
        for _, j in ipairs(rightTranspose.keys) do
            local val = left[i] * rightTranspose[j]
            if val ~= nil then
                row[j] = val
                row.keys[#row.keys+1] = j
            end
        end
        if #row.keys ~= 0 then
            result[i] = row
            result.keys[#result.keys+1] = i
        end
    end

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

    local size = left.dimensions[1] * left.dimensions[2]
    local copy = _linearlyIndexedSparseCopy(left)

    for i = 1, size, 1 do
        local leftVal = left[i]
        local rightVal = right[i]
        if leftVal then
            if rightVal then
                copy[i] = leftVal + rightVal
            else
                copy[i] = leftVal
            end
        elseif rightVal then
            copy[i] = rightVal
        end
    end

    return copy
end

_linearlyIndexedSparseMatrix.__sub = function (left, right)
    if left.dimensions[1] ~= right.dimensions[1] or left.dimensions[2] ~= right.dimensions[2] then
        error("Attempting to add sparse matrices of different sizes.")
    end

    local size = left.dimensions[1] * left.dimensions[2]
    local copy = _linearlyIndexedSparseCopy(left)

    for i = 1, size, 1 do
        local leftVal = left[i]
        local rightVal = right[i]
        if leftVal then
            if rightVal then
                copy[i] = leftVal - rightVal
            else
                copy[i] = leftVal
            end
        elseif rightVal then
            copy[i] = -rightVal
        end
    end

    return copy
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
    if #left ~= #right or #left[1] ~= #right[1] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for i = 1, #left do
        result[i] = {}
        for j = 1, #left[1] do
            result[i][j] = left[i][j] + right[i][j]
        end
    end

    return _matrixFromTableOfTables(result)
end

_matrix.__sub = function (left, right)
    if #left ~= #right or #left[1] ~= #right[1] then
        error("Attempting to add matrices of different sizes.")
    end

    local result = {}

    for i = 1, #left do
        result[i] = {}
        for j = 1, #left[1] do
            result[i][j] = left[i][j] - right[i][j]
        end
    end

    return _matrixFromTableOfTables(result)
end

_matrix.__mul = function (left, right)
    if #left[1] ~= #right then
        error("Attempting to multiply incompatible sparse matrices.")
    end

    local result = setmetatable({}, _matrix)

    for i = 1, #left do
        local row = {}
        for j = 1, #right[1] do
            local val = 0
            for k = 1, #left[1] do
                val = val + left[i][k] * right[k][j]
            end
            row[j] = val
        end
        result[i] = row
    end

    return result
end

_matrix.__tostring = function (matrix)
    local result = "{"

    local length = #matrix

    for i = 1, length - 1 do
        local val = matrix[i]
        result = result .. Tools.list.tostring(val) .. ","
    end

    local val = matrix[length]
    result = result .. Tools.list.tostring(val)

    result = result .. "}"

    return result
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

MatrixAlgebra.sparseMatrix.newNumeric = function (ArrayOfArrays, tol)
    return _sparseMatrixFromNumericArrayOfArrays(ArrayOfArrays, tol)
end

MatrixAlgebra.sparseMatrix.zero = function (n, m)
    local ArrayOfArrays = {}
    for i = 1, n, 1 do
        ArrayOfArrays[i] = {}
        for j = 1, m, 1 do
            ArrayOfArrays[i][j] = 0
        end
    end

    return _sparseMatrixFromArrayOfArrays(ArrayOfArrays)
end

MatrixAlgebra.sparseMatrix.identity = function (n)
    return _sparseMatrixIdentity(n)
end

MatrixAlgebra.sparseMatrix.flatten = function (matrix)
    return _sparseMatrixFlatten(matrix)
end

MatrixAlgebra.sparseMatrix.lu = function (matrix)
    local temp = _sparseMatrixCopy(matrix)
    return _sparseMatrixLU(temp)
end

MatrixAlgebra.sparseMatrix.inverse = function (matrix)
    local temp = _sparseMatrixCopy(matrix)
    return _sparseMatrixInverse(temp)
end

MatrixAlgebra.sparseMatrix.solve = function (matrix, vector)
    local temp = _sparseMatrixCopy(matrix)
    return _sparseMatrixSquareSolve(temp, vector)
end

MatrixAlgebra.sparseMatrix.random = function (n, m, a, b, tol)
    local ArrayOfArrays = {}

    for i = 1, n do
        ArrayOfArrays[i] = {}
        for j = 1, m do
            ArrayOfArrays[i][j] = (b - a) * math.random() + a
        end
    end

    return _sparseMatrixFromNumericArrayOfArrays(ArrayOfArrays, tol)
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

MatrixAlgebra.liSparseMatrix.solve = function (matrix, vector)
    local mat = _linearlyIndexedSparseCopy(matrix)
    return _linearlyIndexedSparseMatrixSquareSolve(mat, vector)
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

MatrixAlgebra.liSparseMatrix.apply = function (matrix, vector)
    local columnVector = _linearlyIndexedSparseTranspose(_linearlyIndexedSparseMatrixFromTableOfTables({vector}))
    columnVector = matrix * columnVector;
    return Tools.list.copy(columnVector)
end

MatrixAlgebra.liSparseMatrix.scale = function (matrix, c)
    local copy = _linearlyIndexedSparseCopy(matrix)
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]
    local size = numberOfRows * numberOfColumns
    if type(c) == "number" then
        for i = 1, size, 1 do
            local v = copy[i]
            if v ~= nil then
                copy[i] = c * v
            end
        end
    elseif type(c) == "table" then
        for i = 1, numberOfRows, 1 do
            for j = 1, numberOfColumns, 1 do
                local v = copy[numberOfRows * (i - 1) + j]
                if v ~= nil then
                    copy[numberOfRows * (i - 1) + j] = c[i] * v
                end
            end
        end
    end
    return copy
end

MatrixAlgebra.liSparseMatrix.toArrayOfArrays = function (matrix)
    local result = {}
    local numberOfRows = matrix.dimensions[1]
    local numberOfColumns = matrix.dimensions[2]

    for i = 1, numberOfRows, 1 do
        result[i] = {}
        for j = 1, numberOfColumns, 1 do
            result[i][j] = matrix[numberOfRows * (i - 1) + j]
        end
    end

    return result
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

MatrixAlgebra.matrix.flatten = function (matrix)
    return _matrixFlatten(matrix)
end

MatrixAlgebra.matrix.lu = function (matrix)
    local temp = _matrixCopy(matrix)
    return _matrixLU(temp)
end

MatrixAlgebra.matrix.inverse = function (matrix)
    local temp = _matrixCopy(matrix)
    return _matrixInverse(temp)
end

MatrixAlgebra.matrix.solve = function (matrix, vector)
    local temp = _matrixCopy(matrix)
    return _matrixSquareSolve(temp, vector)
end

MatrixAlgebra.matrix.random = function (n, m, a, b)
    local ArrayOfArrays = {}

    for i = 1, n do
        ArrayOfArrays[i] = {}
        for j = 1, m do
            ArrayOfArrays[i][j] = (b - a) * math.random() + a
        end
    end

    return _matrixFromTableOfTables(ArrayOfArrays)
end

MatrixAlgebra.matrix.apply = function (matrix, vector)
    local columnVector = _matrixTranspose(_matrixFromTableOfTables({vector}))
    columnVector = matrix * columnVector;

    local result = {}
    for i = 1, #columnVector do
        result[i] = columnVector[i][1]
    end

    return result
end

MatrixAlgebra.matrix.scale = function (matrix, c)
    local copy = _matrixCopy(matrix)
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]
    if type(c) == "number" then
        for i = 1, numberOfRows, 1 do
            for j = 1, numberOfColumns, 1 do
                copy[i][j] = c * copy[i][j]
            end
        end
    elseif type(c) == "table" then
        for i = 1, numberOfRows, 1 do
            for j = 1, numberOfColumns, 1 do
                copy[i][j] = c[i] * copy[i][j]
            end
        end
    end
    return copy
end

MatrixAlgebra.matrix.map = function (matrix, f)
    local copy = _matrixCopy(matrix)
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]
    for i = 1, numberOfRows, 1 do
        for j = 1, numberOfColumns, 1 do
            copy[i][j] = f(copy[i][j])
        end
    end
    return copy
end

MatrixAlgebra.matrix.toArrayOfArrays = function (matrix)
    local result = {}
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    for i = 1, numberOfRows, 1 do
        result[i] = {}
        for j = 1, numberOfColumns, 1 do
            result[i][j] = matrix[i][j]
        end
    end

    return result
end

return MatrixAlgebra
