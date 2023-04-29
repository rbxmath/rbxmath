local Tools = require("src.Tools")
local Scalars = require("src.Scalars")

local Matrices = {}

Matrices.constants = {}
Matrices.constants.STRASSENLIMIT = 256

local Matrix = {
    length = 0,
    width = 0
}

function Matrix:new (list)
    if type(list[1]) == "table" then
        setmetatable(list, self)
        self.__index = self
        list.length = #list
        list.width = #list[1]
        return list
    else
        local matrix = {}
        for key, value in ipairs(list) do
            matrix[key] = {value}
        end
        setmetatable(matrix, self)
        self.__index = self
        matrix.length = #list
        matrix.width = 1
        return matrix
    end
end

--[[

----------------
Matrix Utilities
----------------

]]

function Matrix:copy ()
    local data = {}
    for i = 1, self.length do
        data[i] = {}
        for j = 1, self.width do
            data[i][j] = self[i][j]
        end
    end
    return Matrix:new(data)
end

function Matrix:submatrix (rowStart, rowEnd, columnStart, columnEnd)
    local data = {}
    local row, column = 1, 1
    for i = rowStart, rowEnd, 1 do
        data[row] = {}
        column = 1
        for j = columnStart, columnEnd, 1 do
            data[row][column] = self[i][j]
            column = column + 1
        end
        row = row + 1
    end
    return Matrix:new(data)
end

function Matrix:setSubmatrix (rowStart, rowEnd, columnStart, columnEnd, matrix)
    local row, column = 1, 1
    for i = rowStart, rowEnd, 1 do
        column = 1
        for j = columnStart, columnEnd, 1 do
            self[i][j] = matrix[row][column]
            column = column + 1
        end
        row = row + 1
    end
    return self
end

function Matrix:toScaled (lambda)
    for i = 1, self.length do
        for j = 1, self.width do
            self[i][j] = lambda * self[i][j]
        end
    end
    return self
end

function Matrix:scaled (lambda)
    local copy = self:copy()
    for i = 1, copy.length do
        for j = 1, copy.width do
            copy[i][j] = lambda * copy[i][j]
        end
    end
    return copy
end

function Matrix:padTo (length, width)
    local data = {}
    for i = 1, length, 1 do
        data[i] = {}
        for j = 1, width, 1 do
            data[i][j] = self[i][j] or 0
        end
    end
    return Matrix:new(data)
end

function Matrix:strassenSubdivide ()
    local size = math.max(self.length, self.width)
    if size % 2 == 1 then
        size = size + 1
    end
    local data1 = {}
    local data2 = {}
    local data3 = {}
    local data4 = {}
    size = size / 2
    for i = 1, size, 1 do
        local iPlus = i + size
        data1[i], data2[i], data3[i], data4[i] = {}, {}, {}, {}
        for j = 1, size, 1 do
            local jPlus = j + size
            data1[i][j], data2[i][j], data3[i][j], data4[i][j] = self[i][j] or 0, self[i][jPlus] or 0, self[iPlus][j] or 0, self[iPlus][jPlus] or 0
        end
    end
    return {
        Matrix:new(data1),
        Matrix:new(data2),
        Matrix:new(data3),
        Matrix:new(data4)
    }
end

Matrix.__add = function (left, right)
    if left.length ~= right.length or left.width ~= right.width then
        error("Attempting to add matrices of incompatible dimension!", -1)
    end

    local data = {}
    for i = 1, left.length, 1 do
        data[i] = {}
        for j = 1, left.width, 1 do
            data[i][j] = left[i][j] + right[i][j]
        end
    end

    return Matrix:new(data)
end

Matrix.__sub = function (left, right)
    if left.length ~= right.length or left.width ~= right.width then
        error("Attempting to add matrices of incompatible dimension!", -1)
    end

    local data = {}
    for i = 1, left.length, 1 do
        data[i] = {}
        for j = 1, left.width, 1 do
            data[i][j] = left[i][j] - right[i][j]
        end
    end

    return Matrix:new(data)
end

Matrix.__mul = function (left, right)
    if left.width ~= right.length then
        error("Attempting to multiply matrices of incompatible dimension!", -1)
    end
    local data
    if math.max(left.length, left.width, right.length, right.width) >= Matrices.constants.STRASSENLIMIT then
        --local tic = os.clock()
        local strassenLeft = left:strassenSubdivide()
        local strassenRight = right:strassenSubdivide()
        --print(os.clock() - tic, math.max(left.length, left.width, right.length, right.width))
        local M1 = (strassenLeft[1] + strassenLeft[4]) * (strassenRight[1] + strassenRight[4])
        local M2 = (strassenLeft[3] + strassenLeft[4]) * strassenRight[1]
        local M3 = strassenLeft[1] * (strassenRight[2] - strassenRight[4])
        local M4 = strassenLeft[4] * (strassenRight[3] - strassenRight[1])
        local M5 = (strassenLeft[1] + strassenLeft[2]) * strassenRight[4]
        local M6 = (strassenLeft[3] - strassenLeft[1]) * (strassenRight[1] + strassenRight[2])
        local M7 = (strassenLeft[2] - strassenLeft[4]) * (strassenRight[3] + strassenRight[4])
        local strassenResult = {M1 + M4 - M5 + M7, M3 + M5, M2 + M4, M1 - M2 + M3 + M6}
        data = {}
        for i = 1, M1.length, 1 do
            local iPlus = i + M1.length
            data[i] = {}
            data[iPlus] = {}
            for j = 1, M1.width, 1 do
                local jPlus = j + M1.length
                data[i][j] = strassenResult[1][i][j]
                data[i][jPlus] = strassenResult[2][i][j]
                data[iPlus][j] = strassenResult[3][i][j]
                data[iPlus][jPlus] = strassenResult[4][i][j]
            end
        end
        return Matrix:new(data):submatrix(1, left.length, 1, right.width)
    else
        data = {}
        for i = 1, left.length do
            local row = {}
            for j = 1, right.width do
                local val = 0
                for k = 1, left.width do
                    val = val + left[i][k] * right[k][j]
                end
                row[j] = val
            end
            data[i] = row
        end
        return Matrix:new(data)
    end
end

Matrix.__tostring = function (matrix)
    local result = "{"

    local length = matrix.length

    for i = 1, length - 1 do
        local val = matrix[i]
        result = result .. Tools.list.tostring(val) .. ","
    end

    local val = matrix[length]
    result = result .. Tools.list.tostring(val)

    result = result .. "}"

    return result
end

function _padStringToLength (string, length)
    local result = string
    if (length - #result) % 2 == 1 then
        result = result .. " "
    end
    local width = length - #result
    for i = 1, width / 2 do
        result = " " .. result .. " "
    end
    return result
end

function Matrix:pretty ()
    local length = self.length
    local width = self.width

    local result = "";

    if length == 1 then
        result = "("
        for i = 1, width - 1 do
            result = result .. _padStringToLength(tostring(self[1][i]), 20) .. " "
        end
        result = result .. _padStringToLength(tostring(self[1][width]), 20) .. ")"
        return result
    end

    for i = 1, length do
        if i == 1 then
            result = result .. "/"
        elseif i == length then
            result = result .. "\\"
        else
            result = result .. "|"
        end
        for j = 1, width - 1 do
            result = result .. _padStringToLength(tostring(self[i][j]), 20) .. " "
        end
        result = result .. _padStringToLength(tostring(self[i][width]), 20)
        if i == 1 then
            result = result .. "\\\n"
        elseif i == length then
            result = result .. "/"
        else
            result = result .. "|\n"
        end
    end
    return result
end

--[[

---------------
Common Matrices
---------------

]]

function Matrix:zero (n, m)
    m = m or n
    local data = {}
    for i = 1, n do
        data[i] = {}
        for j = 1, m do
            data[i][j] = 0
        end
    end
    return Matrix:new(data)
end

function Matrix:random (n, m, a, b)
    m = m or n
    local data = {}
    for i = 1, n do
        data[i] = {}
        for j = 1, m do
            data[i][j] = (b - a) * math.random() + a
        end
    end
    return Matrix:new(data)
end

function Matrix:identity (n)
    local data = {}
    for i = 1, n do
        data[i] = {}
        for j = 1, n do
            if i == j then
                data[i][j] = 1
            else
                data[i][j] = 0
            end
        end
    end
    return Matrix:new(data)
end

function Matrix:permutation (permutation)
    local matrix = Matrix:identity(#permutation)
    local data = {}
    for key, value in ipairs(permutation) do
        data[key] = matrix[value]
    end
    return Matrix:new(data)
end

function Matrix:permuted (permutation)
    local data = {}
    for key, value in ipairs(permutation) do
        data[key] = self[value]
    end
    return Matrix:new(data)
end

function Matrix:toPermuted (permutation)
    local matrix = self:copy()
    for key, value in ipairs(permutation) do
        self[key] = matrix[value]
    end
    return self
end

--[[

-----------
Matrix Maps
-----------

]]

function Matrix:transpose ()
    local data = {}
    for j = 1, self.width do
        data[j] = {}
        for i = 1, self.length do
            data[j][i] = self[i][j]
        end
    end
    return Matrix:new(data)
end

function Matrix:toTranspose ()
    local data = {}
    for j = 1, self.width do
        data[j] = {}
        for i = 1, self.length do
            data[j][i] = self[i][j]
        end
    end
    self = data
    return Matrix:new(data)
end

function Matrix:inverse ()
    local length, width = self.length, self.width
    local matrix = self:copy()

    if length ~= width then
        error("Cannot compute inverse of rectangular matrix.", -1)
    end

    local result = matrix:identity(length)

    for i = 1, width - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, length do
            local maxCandidate = matrix[j][i]
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Matrix is not invertible.", -1)
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            result[i], result[maxRow] = result[maxRow], result[i]
        end

        max = matrix[i][i]

        for j = i + 1, length do
            local val = matrix[j][i]
            local valOverMax = val / max
            matrix[j][i] = 0
            for k = 1, width do
                if k > i then
                    matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
                end
                result[j][k] = result[j][k] - result[i][k] * valOverMax
            end
        end
    end

    for i = length, 1, -1 do
        local val = matrix[i][i]
        for j = 1, i - 1 do
            local val1 = matrix[i - j][i]
            for k = 1, width do
                result[i - j][k] = result[i - j][k] - val1 * result[i][k] / val
            end
        end
        for j = 1, width do
            result[i][j] = result[i][j] / val
        end
    end

    return result
end

function Matrix:toInverse ()
    local length, width = self.length, self.width
    local matrix = self:copy()

    if length ~= width then
        error("Cannot compute inverse of rectangular matrix.", -1)
    end

    local result = matrix:identity(length)

    for i = 1, width - 1 do
        local maxRow = i
        local max = matrix[i][i] or 0

        for j = i, length do
            local maxCandidate = matrix[j][i]
            maxCandidate = math.abs(maxCandidate)
            if maxCandidate > max then
                max = maxCandidate
                maxRow = j
            end
        end

        if max == 0 then
            error("Matrix is not invertible.", -1)
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            result[i], result[maxRow] = result[maxRow], result[i]
        end

        max = matrix[i][i]

        for j = i + 1, length do
            local val = matrix[j][i]
            local valOverMax = val / max
            matrix[j][i] = 0
            for k = 1, width do
                if k > i then
                    matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
                end
                result[j][k] = result[j][k] - result[i][k] * valOverMax
            end
        end
    end

    for i = length, 1, -1 do
        local val = matrix[i][i]
        for j = 1, i - 1 do
            local val1 = matrix[i - j][i]
            for k = 1, width do
                result[i - j][k] = result[i - j][k] - val1 * result[i][k] / val
            end
        end
        for j = 1, width do
            result[i][j] = result[i][j] / val
        end
    end
    self = result
    return self
end

--[[

--------------
Linear Systems
--------------

]]

function Matrix:solve (vector)
    local matrix = self:copy()
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    if numberOfRows ~= numberOfColumns then
        error("Cannot solve rectangular system with this function.")
    end

    local columnVector = Matrix:new(vector)

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
            error("Matrix system is not solvable")
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
            error("Matrix system is not solvable")
        end
        columnVector[i][1] = (columnVector[i][1] - temp) / matrix[i][i]
    end

    for i = 1, numberOfRows do
        result[i] = columnVector[i][1]
    end

    return result
end

function Matrix:LUDecomposition () 
    local matrix = self:copy()
    local numberOfRows = #matrix
    local numberOfColumns = #matrix[1]

    if numberOfRows ~= numberOfColumns then
        error("Cannot compute LU of rectangular matrix.")
    end

    local permutation = {}
    for i = 1, numberOfRows do
        permutation[i] = i
    end

    local l = Matrix:identity(numberOfRows)

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
            error("Matrix is not invertible")
        end

        if maxRow ~= i then
            matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
            for k = 1, i - 1 do
                l[i][k], l[maxRow][k] = l[maxRow][k], l[i][k]
            end
            permutation[i], permutation[maxRow] = permutation[maxRow], permutation[i]
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

    return {l, matrix, permutation}
end

--[[

------------------
Matrix Scalar Maps
------------------

]]

function Matrix:determinant ()
    if self.length == 2 and self.width == 2 then
        return self[1][1] * self[2][2] - self[1][2] * self[2][1]
    end
    local LUDecomposition
    local test = pcall(function () LUDecomposition = self:LUDecomposition() end)
    if test == false then
        return 0
    end
    local determinant = 1
    for i = 1, self.length do
        determinant = determinant * LUDecomposition[2][i][i]
    end
    return determinant * math.pow(-1, Tools.combinatorics.inversionNumber(LUDecomposition[3]))
end

function Matrix:trace ()
    local n = math.min(self.length, self.width)
    local sum = 0
    for i = 1, n do
        sum = sum + matrix[i][i]
    end
    return sum
end

function Matrix:dot (matrix)
    return (self:transpose() * matrix):trace()
end

function Matrix:column (i)
    if i > self.width then
        error("Matrix doesn't have " .. tostring(i) .. " columns.")
    end
    local column = {}
    for j = 1, self.width do
        column[j] = self[i][j]
    end
    return column
end

--[[

----------------------
Eigenvalue Computation
----------------------

]]

-- The following is taken from an excellent and thorough book called Fundamentals of Matrix Computation by Watkins
function Matrix:hessenbergForm ()
    local a
    if self.length ~= self.width then
        a = self:padTo(math.max(self.length, self.width))
    else
        a = self:copy()
    end
    local n = a.length
    local b = Matrix:zero(n, 1)
    for k = 1, n - 2 do
        local maxList = {}
        for i = k + 1, n do
            maxList[#maxList + 1] = math.abs(a[i][k])
        end
        local beta = math.max(table.unpack(maxList))
        local gamma = 0
        if beta ~= 0 then
            local sum = 0
            for i = k + 1, n do
                a[i][k] = a[i][k] / beta
                sum = sum + math.pow(a[i][k], 2)
            end
            local tau = math.sqrt(sum)
            if a[k + 1][k] < 0 then tau = -tau end
            local eta = a[k + 1][k] + tau
            a[k + 1][k] = 1
            for i = k + 2, n do
                a[i][k] = a[i][k] / eta
            end
            gamma = eta / tau
            tau = tau * beta

            local temp = a:submatrix(k + 1, n, k + 1, n)
            local temp2 = a:submatrix(k + 1, n, k, k)
            b:setSubmatrix(k + 1, n, 1, 1, (temp2:transpose() * temp):toScaled(-gamma):transpose())
            a:setSubmatrix(k + 1, n, k + 1, n, temp + temp2 * b:submatrix(k + 1, n, 1, 1):transpose())

            temp = a:submatrix(1, n, k + 1, n)
            temp2 = a:submatrix(k + 1, n, k, k)
            b:setSubmatrix(1, n, 1, 1, (temp * temp2):toScaled(-gamma))
            a:setSubmatrix(1, n, k + 1, n, temp + b:submatrix(1, n, 1, 1) * temp2:transpose())
            a:setSubmatrix(k + 1, n, k, k, temp2:toScaled(0))

            a[k + 1][k] = -tau
        end
    end
    return a
end
function Matrix:toHessenbergForm ()
    local a
    if self.length ~= self.width then
        a = self:padTo(math.max(self.length, self.width))
        self = a
    else
        a = self
    end
    local n = a.length
    local b = Matrix:zero(n, 1)
    for k = 1, n - 2 do
        local maxList = {}
        for i = k + 1, n do
            maxList[#maxList + 1] = math.abs(a[i][k])
        end
        local beta = math.max(table.unpack(maxList))
        local gamma = 0
        if beta ~= 0 then
            local sum = 0
            for i = k + 1, n do
                a[i][k] = a[i][k] / beta
                sum = sum + math.pow(a[i][k], 2)
            end
            local tau = math.sqrt(sum)
            if a[k + 1][k] < 0 then tau = -tau end
            local eta = a[k + 1][k] + tau
            a[k + 1][k] = 1
            for i = k + 2, n do
                a[i][k] = a[i][k] / eta
            end
            gamma = eta / tau
            tau = tau * beta

            local temp = a:submatrix(k + 1, n, k + 1, n)
            local temp2 = a:submatrix(k + 1, n, k, k)
            b:setSubmatrix(k + 1, n, 1, 1, (temp2:transpose() * temp):toScaled(-gamma):transpose())
            a:setSubmatrix(k + 1, n, k + 1, n, temp + temp2 * b:submatrix(k + 1, n, 1, 1):transpose())

            temp = a:submatrix(1, n, k + 1, n)
            temp2 = a:submatrix(k + 1, n, k, k)
            b:setSubmatrix(1, n, 1, 1, (temp * temp2):toScaled(-gamma))
            a:setSubmatrix(1, n, k + 1, n, temp + b:submatrix(1, n, 1, 1) * temp2:transpose())
            a:setSubmatrix(k + 1, n, k, k, temp2:toScaled(0))

            a[k + 1][k] = -tau
        end
    end
    return a
end

function Matrix:francisOne (tol, maxIters)
    tol = tol or 10^(-7)
    maxIters = maxIters or 1000
    local matrix
    if self.length == 1 and self.width == 1 then
        return {self[1][1]}
    elseif self.length ~= self.width then
        matrix = self:padTo(math.max(self.length, self.width))
        self = matrix
    else
        matrix = self
    end
    local n = matrix.length
    if n == 2 then
        local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
        if temp1 < 0 then
            return {Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2), Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2)}
        else
            return {(matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2, (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2}
        end
    end
    local iterations = 0
    local temp = 0
    local tempMin = 0
    while iterations < maxIters do
        local min = 10^(12)
        local minimizer = 0
        for i = 1, n - 1 do
            if math.abs(matrix[i + 1][i]) < min then
                min = math.abs(matrix[i + 1][i])
                minimizer = i
            end
            if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
                return Tools.list.join(matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters), matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters))
            end
        end
        local shift = matrix[n][n]
        local u = {matrix[1][1] - shift + math.sqrt(math.pow(matrix[1][1] - shift, 2) + math.pow(matrix[2][1], 2)), matrix[2][1]}
        local gamma = 2 / (math.pow(u[1], 2) + math.pow(u[2], 2))
        u = Matrix:new(u)
        local Q = Matrix:identity(2) - (u * u:transpose()):toScaled(gamma)
        for i = 1, n do
            matrix:setSubmatrix(1, 2, i, i, Q * matrix:submatrix(1, 2, i, i))
        end
        for i = 1, n do
            matrix:setSubmatrix(i, i, 1, 2, matrix:submatrix(i, i, 1, 2) * Q:transpose())
        end
        --[[for k = 1, n - 2 do
            u = matrix:submatrix(k + 1, k + 2, k, k)
            u[1][1] = u[1][1] + math.sqrt(math.pow(u[1][1], 2) + math.pow(u[2][1], 2))
            gamma = 2 / (math.pow(u[1][1], 2) + math.pow(u[2][1], 2))
            Q = Matrix:identity(2) - (u * u:transpose()):toScaled(gamma)
            matrix:setSubmatrix(k + 1, k + 2, k, n, Q * matrix:submatrix(k + 1, k + 2, k, n))
            matrix:setSubmatrix(k, n, k + 1, k + 2, matrix:submatrix(k, n, k + 1, k + 2) * Q:transpose())
        end]]
        matrix:toHessenbergForm()
        iterations = iterations + 1
        temp = minimizer
        tempMin = min
        if iterations == maxIters then
            --print("Francis One failed to converge in " .. tostring(maxIters) .. " iterations! Breaking on " .. tostring(min) .. ".")
            return Tools.list.join(matrix:submatrix(1, minimizer, 1, minimizer):francisOne(tol), matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisOne(tol))
        end
    end
    return matrix
end

function Matrix:francisTwo (tol, maxIters)
    tol = tol or 10^(-7)
    maxIters = maxIters or 10000
    local matrix
    if self.length == 1 and self.width == 1 then
        return {self[1][1]}
    elseif self.length ~= self.width then
        matrix = self:padTo(math.max(self.length, self.width))
        self = matrix
    else
        matrix = self
    end
    local n = matrix.length
    if n == 2 then
        local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
        if temp1 < 0 then
            return {Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2), Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2)}
        else
            return {(matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2, (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2}
        end
    end
    local iterations = 0
    local temp = 0
    local tempMin = 0
    local p = 2
    while iterations < maxIters do
        local min = 10^(12)
        local minimizer = 0
        -- Check if the lower diagonal has any zero elements and split
        for i = 1, n - 1 do
            if math.abs(matrix[i + 1][i]) < min then
                min = math.abs(matrix[i + 1][i])
                minimizer = i
            end
            if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
                return Tools.list.join(matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters), matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters))
            end
        end
        local shift = matrix:submatrix(n - 1, n, n - 1, n):francisOne()
        local u
        local gamma
        -- Compute the first column of shifted matrix
        if type(shift[1]) == "table" then
            local a11 = Scalars.Complex.new(matrix[1][1], 0)
            local a12 = Scalars.Complex.new(matrix[1][2], 0)
            local a21 = Scalars.Complex.new(matrix[2][1], 0)
            local a22 = Scalars.Complex.new(matrix[2][2], 0)
            local rho1 = shift[1]
            local rho2 = shift[2]
            u = {( (a11 - rho1) * (a11 - rho2) + a12 * a21 )[1], ( a21 * ((a11 + a22) - (rho1 + rho2)) )[1], matrix[3][2] * matrix[2][1]}
            tau = Tools.list.norm(u)
            u[1] = u[1] + tau
            tau = Tools.list.norm(u)
            gamma = 2 / math.pow(tau, 2)
        else
            local a11 = matrix[1][1]
            local a12 = matrix[1][2]
            local a21 = matrix[2][1]
            local a22 = matrix[2][2]
            local rho1 = shift[1]
            local rho2 = shift[2]
            u = {( (a11 - rho1) * (a11 - rho2) + a12 * a21 ), ( a21 * ((a11 + a22) - (rho1 + rho2)) ), matrix[3][2] * matrix[2][1]}
            tau = Tools.list.norm(u)
            u[1] = u[1] + tau
            tau = Tools.list.norm(u)
            gamma = 2 / math.pow(tau, 2)
        end
        -- Compute reflector
        u = Matrix:new(u)
        local Q = Matrix:identity(u.length) - (u * u:transpose()):toScaled(gamma)
        for i = 1, n do
            matrix:setSubmatrix(1, p + 1, i, i, Q * matrix:submatrix(1, p + 1, i, i))
        end
        for i = 1, n do
            matrix:setSubmatrix(i, i, 1, p + 1, matrix:submatrix(i, i, 1, p + 1) * Q:transpose())
        end
        --[[for k = 1, n - p - 1 do
            u = matrix:submatrix(k + 1, k + p + 1, k, k)
            local tau = math.sqrt((u:transpose() * u)[1][1])
            u[1][1] = u[1][1] + tau
            tau = math.sqrt((u:transpose() * u)[1][1])
            gamma = 2 / math.pow(tau, 2)
            Q = Matrix:identity(p + 1) - (u * u:transpose()):toScaled(gamma)
            matrix:setSubmatrix(k + 1, k + p + 1, k, n, Q * matrix:submatrix(k + 1, k + p + 1, k, n))
            matrix:setSubmatrix(k, n, k + 1, k + p + 1, matrix:submatrix(k, n, k + 1, k + p + 1) * Q:transpose())
        end]]
        matrix:toHessenbergForm()
        iterations = iterations + 1
        temp = minimizer
        tempMin = min
        if iterations == maxIters then
            --print("Francis Two failed to converge in " .. tostring(maxIters) .. " iterations! Breaking on " .. tostring(min) .. ".")
            return Tools.list.join(matrix:submatrix(1, minimizer, 1, minimizer):francisTwo(tol), matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisTwo(tol))
        end
    end
    return matrix
end

function Matrix:francisThree (tol, maxIters)
    tol = tol or 10^(-10)
    maxIters = maxIters or 10000
    local matrix
    if self.length == 1 and self.width == 1 then
        return {self[1][1]}
    elseif self.length ~= self.width then
        matrix = self:padTo(math.max(self.length, self.width))
        self = matrix
    else
        matrix = self
    end
    local n = matrix.length
    if n == 2 then
        local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
        if temp1 < 0 then
            return {Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2), Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2)}
        else
            return {(matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2, (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2}
        end
    elseif n == 3 then
        return self:francisTwo(tol, maxIters)
    end
    local iterations = 0
    local temp = 0
    local tempMin = 0
    local p = 3
    while iterations < maxIters do
        local min = 10^(12)
        local minimizer = 0
        -- Check if the lower diagonal has any zero elements and split
        for i = 1, n - 1 do
            if math.abs(matrix[i + 1][i]) < min then
                min = math.abs(matrix[i + 1][i])
                minimizer = i
            end
            if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
                return Tools.list.join(matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters), matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters))
            end
        end
        local shift = matrix:submatrix(n - 2, n, n - 2, n):francisOne(tol * 10, maxIters)
        Scalars.Complex.sort(shift)
        local u
        local gamma
        -- Compute the first column of shifted matrix
        if type(shift[1]) == "table" or type(shift[2]) == "table" then
            local a11 = Scalars.Complex.new(matrix[1][1], 0)
            local a12 = Scalars.Complex.new(matrix[1][2], 0)
            local a13 = Scalars.Complex.new(matrix[1][3], 0)
            local a21 = Scalars.Complex.new(matrix[2][1], 0)
            local a22 = Scalars.Complex.new(matrix[2][2], 0)
            local a23 = Scalars.Complex.new(matrix[2][3], 0)
            local a32 = Scalars.Complex.new(matrix[3][2], 0)
            local a33 = Scalars.Complex.new(matrix[3][3], 0)
            local rho1, rho2, rho3
            if type(shift[1]) == "table" then
                rho1 = shift[1]
            else
                rho1 = Scalars.Complex.new(shift[1], 0)
            end
            if type(shift[2]) == "table" then
                rho2 = shift[2]
            else
                rho2 = Scalars.Complex.new(shift[2], 0)
            end
            if type(shift[3]) == "table" then
                rho3 = shift[3]
            else
                rho3 = Scalars.Complex.new(shift[3], 0)
            end
            u = {(a21 * (a13 * a32 + a12 * (a11 + a22 - rho2 - rho3)) + (a11 - rho1) * (a11 * a11 + a12 * a21 + rho2 * rho3 - a11 * (rho2 + rho3)))[1], 
                (a21 * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho1) * (a11 + a22 - rho2 - rho3) + rho2 * rho3 - a22 * (rho2 + rho3)))[1], 
                (a21 * a32 * (a11 + a22 + a33 - rho1 - rho2 - rho3))[1],
                matrix[2][1] * matrix[3][2] * matrix[4][3]
            }
            tau = Tools.list.norm(u)
            u[1] = u[1] + tau
            tau = Tools.list.norm(u)
            gamma = 2 / math.pow(tau, 2)
        else
            local a11 = matrix[1][1]
            local a12 = matrix[1][2]
            local a13 = matrix[1][3]
            local a21 = matrix[2][1]
            local a22 = matrix[2][2]
            local a23 = matrix[2][3]
            local a32 = matrix[3][2]
            local a33 = matrix[3][3]
            local rho1 = shift[1]
            local rho2 = shift[2]
            local rho3 = shift[3]
            u = {
                (a21 * (a13 * a32 + a12 * (a11 + a22 - rho2 - rho3)) + (a11 - rho1) * (a11^2 + 
                a12 * a21 + rho2 * rho3 - a11 * (rho2 + rho3))), 
                (a21 * (a12 * a21 + a22^2 + 
                a23 * a32 + (a11 - rho1) * (a11 + a22 - rho2 - rho3) + rho2 * rho3 - 
                a22 * (rho2 + rho3))), 
                (a21 * a32 * (a11 + a22 + a33 - rho1 - rho2 - rho3)),
                matrix[2][1] * matrix[3][2] * matrix[4][3]
            }
            tau = Tools.list.norm(u)
            u[1] = u[1] + tau
            tau = Tools.list.norm(u)
            gamma = 2 / math.pow(tau, 2)
        end
        -- Compute reflector
        u = Matrix:new(u)
        local Q = Matrix:identity(u.length) - (u * u:transpose()):toScaled(gamma)
        for i = 1, n do
            matrix:setSubmatrix(1, p + 1, i, i, Q * matrix:submatrix(1, p + 1, i, i))
        end
        for i = 1, n do
            matrix:setSubmatrix(i, i, 1, p + 1, matrix:submatrix(i, i, 1, p + 1) * Q:transpose())
        end
        --[[for k = 1, n - p - 1 do
            u = matrix:submatrix(k + 1, k + p + 1, k, k)
            local tau = math.sqrt((u:transpose() * u)[1][1])
            u[1][1] = u[1][1] + tau
            tau = math.sqrt((u:transpose() * u)[1][1])
            gamma = 2 / math.pow(tau, 2)
            Q = Matrix:identity(p + 1) - (u * u:transpose()):toScaled(gamma)
            matrix:setSubmatrix(k + 1, k + p + 1, k, n, Q * matrix:submatrix(k + 1, k + p + 1, k, n))
            matrix:setSubmatrix(k, n, k + 1, k + p + 1, matrix:submatrix(k, n, k + 1, k + p + 1) * Q:transpose())
        end]]
        matrix:toHessenbergForm()
        iterations = iterations + 1
        temp = minimizer
        tempMin = min
        if iterations == maxIters then
            print("Francis Three failed to converge in " .. tostring(maxIters) .. " iterations! Breaking on " .. tostring(min) .. ".")
            return Tools.list.join(matrix:submatrix(1, minimizer, 1, minimizer):francisThree(tol), matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisThree(tol))
        end
    end
    return matrix
end

function Matrix:francisSix (tol, maxIters)
    -- Set Defaults
    tol = tol or 10^(-12)
    maxIters = maxIters or 10000
    -- Create Matrix
    local matrix
    if self.length == 1 and self.width == 1 then
        return {self[1][1]}
    elseif self.length ~= self.width then
        matrix = self:padTo(math.max(self.length, self.width))
        self = matrix
    else
        matrix = self
    end
    -- Handle the case where n is too small
    local n = matrix.length
    if n == 2 then
        local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
        if temp1 < 0 then
            return {Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2), Scalars.Complex.new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2)}
        else
            return {(matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2, (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2}
        end
    elseif n == 3 then
        return self:francisTwo(tol, maxIters)
    elseif n < 7 then
        return self:francisThree(tol, maxIters)
    end
    -- Setup variables
    local iterations = 0
    local temp = 0
    local tempMin = 0
    local p = 6
    while iterations < maxIters do
        local min = 10^(12)
        local minimizer = 0
        -- Check if the lower diagonal has any zero elements and split
        for i = 1, n - 1 do
            if math.abs(matrix[i + 1][i]) < min then
                min = math.abs(matrix[i + 1][i])
                minimizer = i
            end
            if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
                return Tools.list.join(matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters), matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters))
            end
        end
        local shift = matrix:submatrix(n - p + 1, n, n - p + 1, n):francisThree(tol * 10, maxIters)
        Scalars.Complex.sort(shift)
        local u
        local gamma
        -- Compute the first column of shifted matrix
        local a11, a12, a13, a14, a15, a16, a17 = Scalars.Complex.new(matrix[1][1], 0), Scalars.Complex.new(matrix[1][2], 0), Scalars.Complex.new(matrix[1][3], 0), Scalars.Complex.new(matrix[1][4], 0), Scalars.Complex.new(matrix[1][5], 0), Scalars.Complex.new(matrix[1][6], 0), Scalars.Complex.new(matrix[1][7], 0)
        local a21, a22, a23, a24, a25, a26, a27 = Scalars.Complex.new(matrix[2][1], 0), Scalars.Complex.new(matrix[2][2], 0), Scalars.Complex.new(matrix[2][3], 0), Scalars.Complex.new(matrix[2][4], 0), Scalars.Complex.new(matrix[2][5], 0), Scalars.Complex.new(matrix[2][6], 0), Scalars.Complex.new(matrix[2][7], 0)
        local a32, a33, a34, a35, a36, a37 = Scalars.Complex.new(matrix[3][2], 0), Scalars.Complex.new(matrix[3][3], 0), Scalars.Complex.new(matrix[3][4], 0), Scalars.Complex.new(matrix[3][5], 0), Scalars.Complex.new(matrix[3][6], 0), Scalars.Complex.new(matrix[3][7], 0)
        local a43, a44, a45, a46, a47 = Scalars.Complex.new(matrix[4][3], 0), Scalars.Complex.new(matrix[4][4], 0), Scalars.Complex.new(matrix[4][5], 0), Scalars.Complex.new(matrix[4][6], 0), Scalars.Complex.new(matrix[4][7], 0)
        local a54, a55, a56, a57 = Scalars.Complex.new(matrix[5][4], 0), Scalars.Complex.new(matrix[5][5], 0), Scalars.Complex.new(matrix[5][6], 0), Scalars.Complex.new(matrix[5][7], 0)
        local a65, a66, a67 = Scalars.Complex.new(matrix[6][5], 0), Scalars.Complex.new(matrix[6][6], 0), Scalars.Complex.new(matrix[6][7], 0)
        local a76, a77 = Scalars.Complex.new(matrix[7][6], 0), Scalars.Complex.new(matrix[7][7], 0)
        local rho1, rho2, rho3, rho4, rho5, rho6 = Scalars.Complex.new(shift[1]), Scalars.Complex.new(shift[2]), Scalars.Complex.new(shift[3]), Scalars.Complex.new(shift[4]), Scalars.Complex.new(shift[5]), Scalars.Complex.new(shift[6])
        u = {
            (a21*((a22 - rho2)*((a22 - rho3)*((a22 - rho4)*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                a32*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a12*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a32*(a23*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + (a33 - rho4)*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a43*(a12*a24 + a13*a34 + a15*a54 + a14*(a11 + a44 - rho5 - rho6)) + 
                a13*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a12*(a21*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a11 - rho4)*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6)))) + 
                a32*(a23*((a22 - rho4)*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                a32*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a12*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                (a33 - rho3)*(a23*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a33 - rho4)*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a43*(a12*a24 + a13*a34 + a15*a54 + a14*(a11 + a44 - rho5 - rho6)) + 
                a13*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a43*(a24*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + a34*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                (a44 - rho4)*(a12*a24 + a13*a34 + a15*a54 + a14*(a11 + a44 - rho5 - rho6)) + 
                a54*(a12*a25 + a13*a35 + a14*a45 + a16*a65 + a15*(a11 + a55 - rho5 - rho6)) + 
                a14*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a13*(a21*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a11 - rho4)*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6)))) + 
                a12*(a21*((a22 - rho4)*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                a32*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a12*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                (a11 - rho3)*(a21*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a11 - rho4)*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))))) + 
                (a11 - rho1)*(a21*((a22 - rho3)*((a22 - rho4)*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                a32*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a12*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a32*(a23*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + (a33 - rho4)*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a43*(a12*a24 + a13*a34 + a15*a54 + a14*(a11 + a44 - rho5 - rho6)) + 
                a13*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                a12*(a21*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a11 - rho4)*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6)))) + 
                (a11 - rho2)*(a21*((a22 - rho4)*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                a32*(a12*a23 + a14*a43 + a13*(a11 + a33 - rho5 - rho6)) + 
                a12*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))) + 
                (a11 - rho3)*(a21*(a13*a32 + a12*(a11 + a22 - rho5 - rho6)) + 
                (a11 - rho4)*(a11 * a11 + a12*a21 + rho5*rho6 - a11*(rho5 + rho6))))))[1]
            ,
            (a21*(a12*a21*(a32*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + a12*a21*(a11 + a22 - rho5 - rho6) + 
                (a22 - rho4)*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6)) + 
                (a11 - rho3)*(a12*a21 + a22 * a22 + a23*a32 + (a11 - rho4)*(a11 + a22 - rho5 - rho6) + rho5*rho6 - a22*(rho5 + rho6)))
                + a32*(a13*a21*(a12*a21 + a22 * a22 + a23*a32 + (a11 - rho4)*(a11 + a22 - rho5 - rho6) + rho5*rho6 - 
                a22*(rho5 + rho6)) + (a33 - rho3)*((a33 - rho4)*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + 
                a43*(a14*a21 + a23*a34 + a25*a54 + a24*(a22 + a44 - rho5 - rho6)) + a13*a21*(a11 + a22 - rho5 - rho6) + 
                a23*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6))) + 
                a43*(a34*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + 
                (a44 - rho4)*(a14*a21 + a23*a34 + a25*a54 + a24*(a22 + a44 - rho5 - rho6)) + 
                a54*(a15*a21 + a23*a35 + a24*a45 + a26*a65 + a25*(a22 + a55 - rho5 - rho6)) + a14*a21*(a11 + a22 - rho5 - rho6) + 
                a24*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6))) + 
                a23*(a32*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + a12*a21*(a11 + a22 - rho5 - rho6) + 
                (a22 - rho4)*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6)))) + 
                (a22 - rho2)*(a12*a21*(a12*a21 + a22 * a22 + a23*a32 + (a11 - rho4)*(a11 + a22 - rho5 - rho6) + rho5*rho6 - 
                a22*(rho5 + rho6)) + a32*((a33 - rho4)*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + 
                a43*(a14*a21 + a23*a34 + a25*a54 + a24*(a22 + a44 - rho5 - rho6)) + a13*a21*(a11 + a22 - rho5 - rho6) + 
                a23*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6))) + 
                (a22 - rho3)*(a32*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + a12*a21*(a11 + a22 - rho5 - rho6) + 
                (a22 - rho4)*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6)))) + 
                (a11 - rho1)*(a12*a21*(a12*a21 + a22 * a22 + a23*a32 + (a11 - rho4)*(a11 + a22 - rho5 - rho6) + rho5*rho6 - 
                a22*(rho5 + rho6)) + a32*((a33 - rho4)*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + 
                a43*(a14*a21 + a23*a34 + a25*a54 + a24*(a22 + a44 - rho5 - rho6)) + a13*a21*(a11 + a22 - rho5 - rho6) + 
                a23*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6))) + 
                (a22 - rho3)*(a32*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + a12*a21*(a11 + a22 - rho5 - rho6) + 
                (a22 - rho4)*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6))) + 
                (a11 - rho2)*(a32*(a13*a21 + a24*a43 + a23*(a22 + a33 - rho5 - rho6)) + a12*a21*(a11 + a22 - rho5 - rho6) + 
                (a22 - rho4)*(a12*a21 + a22 * a22 + a23*a32 + rho5*rho6 - a22*(rho5 + rho6)) + 
                (a11 - rho3)*(a12*a21 + a22 * a22 + a23*a32 + (a11 - rho4)*(a11 + a22 - rho5 - rho6) + rho5*rho6 - 
                a22*(rho5 + rho6))))))[1]
            ,
            (a21*a32*(a13*a21*a32*(a11 + a22 + a33 - rho4 - rho5 - rho6) + 
                a23*a32*(a12*a21 + a23*a32 + a33 * a33 + a34*a43 + (a22 - rho4)*(a22 + a33 - rho5 - rho6) + rho5*rho6 - 
                a33*(rho5 + rho6)) + a12*a21*(a12*a21 + a23*a32 + a33 * a33 + a34*a43 + (a22 - rho4)*(a22 + a33 - rho5 - rho6) + 
                (a11 - rho3)*(a11 + a22 + a33 - rho4 - rho5 - rho6) + rho5*rho6 - a33*(rho5 + rho6)) + 
                a43*(a14*a21*a32 + (a44 - rho4)*(a24*a32 + a35*a54 + a34*(a33 + a44 - rho5 - rho6)) + 
                a54*(a25*a32 + a34*a45 + a36*a65 + a35*(a33 + a55 - rho5 - rho6)) + a24*a32*(a22 + a33 - rho5 - rho6) + 
                a34*(a23*a32 + a33 * a33 + a34*a43 + rho5*rho6 - a33*(rho5 + rho6))) + 
                (a33 - rho3)*(a13*a21*a32 + a43*(a24*a32 + a35*a54 + a34*(a33 + a44 - rho5 - rho6)) + a23*a32*(a22 + a33 - rho5 - rho6) + 
                (a33 - rho4)*(a23*a32 + a33 * a33 + a34*a43 + rho5*rho6 - a33*(rho5 + rho6))) + 
                (a22 - rho2)*(a13*a21*a32 + a43*(a24*a32 + a35*a54 + a34*(a33 + a44 - rho5 - rho6)) + a23*a32*(a22 + a33 - rho5 - rho6) + 
                a12*a21*(a11 + a22 + a33 - rho4 - rho5 - rho6) + 
                (a33 - rho4)*(a23*a32 + a33 * a33 + a34*a43 + rho5*rho6 - a33*(rho5 + rho6)) + 
                (a22 - rho3)*(a12*a21 + a23*a32 + a33 * a33 + a34*a43 + (a22 - rho4)*(a22 + a33 - rho5 - rho6) + rho5*rho6 - 
                a33*(rho5 + rho6))) + (a11 - rho1)*(a13*a21*a32 + a43*(a24*a32 + a35*a54 + a34*(a33 + a44 - rho5 - rho6)) + 
                a23*a32*(a22 + a33 - rho5 - rho6) + a12*a21*(a11 + a22 + a33 - rho4 - rho5 - rho6) + 
                (a33 - rho4)*(a23*a32 + a33 * a33 + a34*a43 + rho5*rho6 - a33*(rho5 + rho6)) + 
                (a22 - rho3)*(a12*a21 + a23*a32 + a33 * a33 + a34*a43 + (a22 - rho4)*(a22 + a33 - rho5 - rho6) + rho5*rho6 - 
                a33*(rho5 + rho6)) + (a11 - rho2)*(a12*a21 + a23*a32 + a33 * a33 + a34*a43 + 
                (a22 - rho4)*(a22 + a33 - rho5 - rho6) + (a11 - rho3)*(a11 + a22 + a33 - rho4 - rho5 - rho6) + rho5*rho6 - 
                a33*(rho5 + rho6)))))[1]
            ,
            (a21*a32*a43*(a13*a21*a32 + a24*a32*a43 + a54*(a35*a43 + a46*a65 + a45*(a44 + a55 - rho5 - rho6)) + 
                a34*a43*(a33 + a44 - rho5 - rho6) + a23*a32*(a22 + a33 + a44 - rho4 - rho5 - rho6) + 
                a12*a21*(a11 + a22 + a33 + a44 - rho3 - rho4 - rho5 - rho6) + 
                (a44 - rho4)*(a34*a43 + a44 * a44 + a45*a54 + rho5*rho6 - a44*(rho5 + rho6)) + 
                (a33 - rho3)*(a23*a32 + a34*a43 + a44 * a44 + a45*a54 + (a33 - rho4)*(a33 + a44 - rho5 - rho6) + rho5*rho6 - 
                a44*(rho5 + rho6)) + (a22 - rho2)*(a12*a21 + a23*a32 + a34*a43 + a44 * a44 + a45*a54 + 
                (a33 - rho4)*(a33 + a44 - rho5 - rho6) + (a22 - rho3)*(a22 + a33 + a44 - rho4 - rho5 - rho6) + rho5*rho6 - 
                a44*(rho5 + rho6)) + (a11 - rho1)*(a12*a21 + a23*a32 + a34*a43 + a44 * a44 + a45*a54 + 
                (a33 - rho4)*(a33 + a44 - rho5 - rho6) + (a22 - rho3)*(a22 + a33 + a44 - rho4 - rho5 - rho6) + 
                (a11 - rho2)*(a11 + a22 + a33 + a44 - rho3 - rho4 - rho5 - rho6) + rho5*rho6 - a44*(rho5 + rho6))))[1]
            ,
            (a21*a32*a43*a54*(a12*a21 + a23*a32 + a34*a43 + a45*a54 + a55 * a55 + a56*a65 + (a44 - rho4)*(a44 + a55 - rho5 - rho6) + 
                (a33 - rho3)*(a33 + a44 + a55 - rho4 - rho5 - rho6) + (a22 - rho2)*(a22 + a33 + a44 + a55 - rho3 - rho4 - rho5 - rho6) + 
                (a11 - rho1)*(a11 + a22 + a33 + a44 + a55 - rho2 - rho3 - rho4 - rho5 - rho6) + rho5*rho6 - a55*(rho5 + rho6)))[1]
            ,
            (a21*a32*a43*a54*a65*(a11 + a22 + a33 + a44 + a55 + a66 - rho1 - rho2 - rho3 - rho4 - rho5 - rho6))[1]
            ,
            (a21*a32*a43*a54*a65*a76)[1]
        }
        tau = Tools.list.norm(u)
        u[1] = u[1] + tau
        tau = Tools.list.norm(u)
        gamma = 2 / math.pow(tau, 2)
        -- Compute reflector
        u = Matrix:new(u)
        local Q = Matrix:identity(u.length) - (u * u:transpose()):toScaled(gamma)
        for i = 1, n do
            matrix:setSubmatrix(1, p + 1, i, i, Q * matrix:submatrix(1, p + 1, i, i))
        end
        for i = 1, n do
            matrix:setSubmatrix(i, i, 1, p + 1, matrix:submatrix(i, i, 1, p + 1) * Q:transpose())
        end
        --[[for k = 1, n - p - 1 do
            u = matrix:submatrix(k + 1, k + p + 1, k, k)
            local tau = math.sqrt((u:transpose() * u)[1][1])
            u[1][1] = u[1][1] + tau
            tau = math.sqrt((u:transpose() * u)[1][1])
            gamma = 2 / math.pow(tau, 2)
            Q = Matrix:identity(p + 1) - (u * u:transpose()):toScaled(gamma)
            matrix:setSubmatrix(k + 1, k + p + 1, k, n, Q * matrix:submatrix(k + 1, k + p + 1, k, n))
            matrix:setSubmatrix(k, n, k + 1, k + p + 1, matrix:submatrix(k, n, k + 1, k + p + 1) * Q:transpose())
        end]]
        matrix:toHessenbergForm()
        iterations = iterations + 1
        temp = minimizer
        tempMin = min
        if iterations == maxIters then
            print("Francis Six failed to converge in " .. tostring(maxIters) .. " iterations! Breaking on " .. tostring(min) .. ".")
            return Tools.list.join(matrix:submatrix(1, minimizer, 1, minimizer):francisSix(tol), matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisSix(tol))
        end
    end
    return matrix
end

function Matrix:eigenvalues ()
    local matrix = self:hessenbergForm()
    local eigList = matrix:francisSix()
    Scalars.Complex.sort(eigList, "d")
    return eigList
end

Matrices.constants.ComputeStrassenLimit = function ()
    local strassenLimit = Matrices.constants.STRASSENLIMIT
    for i = 1, 100, 1 do
        local data = {}
        for j = 1, math.pow(2, i), 1 do
            data[j] = {}
            for k = 1, math.pow(2, i) do
                data[j][k] = 1
            end
        end
        local matrix = Matrix:new(data)
        Matrices.constants.STRASSENLIMIT = math.floor(math.pow(2, i))
        local tic = os.clock()
        for j = 1, 10, 1 do
            local store = matrix * matrix
        end
        local time = (os.clock() - tic)
        Matrices.constants.STRASSENLIMIT = math.floor(math.pow(2, i + 1))
        tic = os.clock()
        for j = 1, 10, 1 do
            local store = matrix * matrix
        end
        local time2 = (os.clock() - tic)
        if time < time2 then
            print(time, time2)
            Matrices.constants.STRASSENLIMIT = math.floor(math.pow(2, i))
            return nil
        end
    end
end

Matrices.Matrix = Matrix

return Matrices