--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

--[[
   +--------------------------------------------------+
   |                Lua Math: Matrices                |
   +--------------------------------------------------+
   |This portion of the library focuses on the        |
   |creation of many crucial objects in linear        |
   |algebra. For most users, the Matrix object will be|
   |the most useful. Matrices can be multiplied,      |
   |added, and subtracted with the standard symbols,  |
   |and their contents can be accessed as you would   |
   |access information from a table of tables. For the|
   |adventorous, there are the ComplexMatrix and      |
   |SparseMatrix objects. Each of these has their own |
   |quirks, but they are all compatible with the      |
   |standard operations (+, -, *).                    |
   +--------------------------------------------------+
]]

local Tools = require(script.Parent.Tools)
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
   type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object
local Scalars = require(script.Parent.Scalars)
local Complex = Scalars.Complex
local Vectors = require(script.Parent.Vectors)
local Numerics = require(script.Parent.Numerics)

local Matrices = {}

Matrices.constants = {}
Matrices.constants.STRASSENLIMIT = 256
Matrices.constants.SPARSESTRASSENLIMIT = 256
Matrices.constants.COMPLEXSTRASSENLIMIT = 64

--[[
Throughout the library we will use the following conventions:

 - functions that start with "to" or "set" will mutate the object 
     all other functions will leave the matrix unchanged
]]--

local Matrix = {
   length = 0,
   width = 0,
}

function Matrix:new(list: Tensor | Vector): Tensor
   if type(list[1]) == "table" then
      setmetatable(list, self)
      self.__index = self
      list.length = #list
      list.width = #list[1]
      return list
   else
      local matrix = {}
      for key, value in ipairs(list) do
	 matrix[key] = { value }
      end
      setmetatable(matrix, self)
      self.__index = self
      matrix.length = #list
      matrix.width = 1
      return matrix
   end
end

--[[
   +--------------------------------------------------+
   |                 Matrix Utilities                 |
   +--------------------------------------------------+
   |This section contains many useful functions       |
   |including those to copy and manipulate matrices.  |
   |Functions of the form "to_____" or "set_____" will|
   |change the underlying matrix, while others will   |
   |return a shallow copy.                            |
   +--------------------------------------------------+
]]

function Matrix:copy(): Tensor
   local data = {}
   for i = 1, self.length do
      data[i] = {}
      for j = 1, self.width do
	 data[i][j] = self[i][j]
      end
   end
   return Matrix:new(data)
end

function Matrix:submatrix(rowStart: number, rowEnd: number, columnStart: number, columnEnd: number): Tensor
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

function Matrix:setSubmatrix(rowStart: number, rowEnd: number, columnStart: number, columnEnd: number, matrix: Tensor)
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

function Matrix:addBand(value: number, position: number)
   local cp = self:copy()
   if position == 0 then
      for i = 1, math.min(self.length, self.width) do
	 cp[i][i] = value
      end
   elseif position > 0 then
      for i = 1, math.min(self.length, self.width - position) do
	 cp.data[i][i + position] = value
      end
   else
      for i = 1, math.min(self.length + position, self.width) do
	 cp.data[i - position][i] = value
      end
   end
   return cp
end

function Matrix:toAddedBand(value: number, position: number)
   if position == 0 then
      for i = 1, math.min(self.length, self.width) do
	 self[i][i] = value
      end
   elseif position > 0 then
      for i = 1, math.min(self.length, self.width - position) do
	 self[i][i + position] = value
      end
   else
      for i = 1, math.min(self.length + position, self.width) do
	 self[i - position][i] = value
      end
   end
   return self
end

function Matrix:toScaled(lambda: number): Tensor
   for i = 1, self.length do
      for j = 1, self.width do
	 self[i][j] = lambda * self[i][j]
      end
   end
   return self
end

function Matrix:scaled(lambda: number): Tensor
   local copy = self:copy()
   for i = 1, copy.length do
      for j = 1, copy.width do
	 copy[i][j] = lambda * copy[i][j]
      end
   end
   return copy
end

function Matrix:padTo(length: number, width: number): Tensor
   local data = {}
   for i = 1, length, 1 do
      data[i] = {}
      for j = 1, width, 1 do
	 data[i][j] = self[i][j] or 0
      end
   end
   return Matrix:new(data)
end

function Matrix:strassenSubdivide(): Array<Tensor>
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
	 data1[i][j], data2[i][j], data3[i][j], data4[i][j] =
	    self[i][j] or 0, self[i][jPlus] or 0, self[iPlus][j] or 0, self[iPlus][jPlus] or 0
      end
   end
   return {
      Matrix:new(data1),
      Matrix:new(data2),
      Matrix:new(data3),
      Matrix:new(data4),
   }
end

function Matrix:column(i: number): Vector
   if i > self.width then
      error("Matrix doesn't have " .. tostring(i) .. " columns.")
   end
   local column = {}
   for j = 1, self.length do
      column[j] = self[j][i]
   end
   return column
end

function Matrix:getColumn(i: number): Vector
   return self:column(i)
end

function Matrix:row(i: number): Vector
   if i > self.length then
      error("Matrix doesn't have " .. tostring(i) .. " rows.")
   end
   local row = {}
   for j = 1, self.width do
      row[j] = self[i][j]
   end
   return row
end

function Matrix:getRow(i: number): Vector
   return self:row(i)
end

local function _padStringToLength(string: string, length: number): string
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

function Matrix:pretty(n: number, m: number): string
   local length = self.length
   local width = self.width

   n = n or 20
   m = m or 20

   local result = ""

   if length == 1 then
      result = "("
      for i = 1, width - 1 do
	 result = result .. _padStringToLength(string.sub(tostring(self[1][i]), 1, n), m) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self[1][width]), m), 1, n) .. ")"
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
	 result = result .. _padStringToLength(string.sub(tostring(self[i][j]), 1, n), m) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self[i][width]), 1, n), m)
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
   +--------------------------------------------------+
   |                Matrix Metamethods                |
   +--------------------------------------------------+
   |This section contains all of the metamethods for  |
   |matrices. Addition and subtraction are relatively |
   |standard, but multiplication is an implementation |
   |of Strassen's method. The size of matrix at which |
   |Strassen multiplication will be used is set in    |
   |Matrices.constant.STRASSENLIMIT.                  |
   +--------------------------------------------------+
]]

function Matrix.__add(left: Tensor, right: Tensor): Tensor
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

function Matrix.__sub(left: Tensor, right: Tensor): Tensor
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

function Matrix.__mul(left: Tensor, right: Tensor): Tensor
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
      local strassenResult = { M1 + M4 - M5 + M7, M3 + M5, M2 + M4, M1 - M2 + M3 + M6 }
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

function Matrix.__tostring(matrix: Tensor): string
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

function Matrix.__eq(left: Tensor, right: Tensor): boolean
   if left.length ~= right.length or left.width ~= right.width then
      return false
   else
      for i = 1, left.length do
	 for j = 1, left.width do
	    if left[i][j] ~= right[i][j] then
	       return false
	    end
	 end
      end
   end
   return true
end

--[[
   +--------------------------------------------------+
   |                 Common Matrices                  |
   +--------------------------------------------------+
   |This section contains methods for constructing    |
   |many common matrices such as the identity and zero|
   |matrices.                                         |
   +--------------------------------------------------+
]]

function Matrix:zero(n: number, m: number): Tensor
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

function Matrix:random(n: number, m: number, a: number, b: number): Tensor
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

function Matrix:identity(n: number): Tensor
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

function Matrix:permutation(permutation: Vector): Tensor
   local matrix = Matrix:identity(#permutation)
   local data = {}
   for key, value in ipairs(permutation) do
      data[key] = matrix[value]
   end
   return Matrix:new(data)
end

function Matrix:permuted(permutation: Vector): Tensor
   local data = {}
   for key, value in ipairs(permutation) do
      data[key] = self[value]
   end
   return Matrix:new(data)
end

function Matrix:toPermuted(permutation: Vector): Tensor
   local matrix = self:copy()
   for key, value in ipairs(permutation) do
      self[key] = matrix[value]
   end
   return self
end

function Matrix:getColumnVector(i: number): Vector
   local data = {}
   for j = 1, self.length do
      data[j] = self[j][i]
   end
   return data
end

--[[
   +--------------------------------------------------+
   |                   Matrix Maps                    |
   +--------------------------------------------------+
   |This section contains the common maps that send   |
   |matrices to matrices. This includes functions like|
   |the transpose and inverse.                        |
   +--------------------------------------------------+
]]

function Matrix:transpose()
   local data = {}
   for j = 1, self.width do
      data[j] = {}
      for i = 1, self.length do
	 data[j][i] = self[i][j]
      end
   end
   return Matrix:new(data)
end

function Matrix:toTranspose()
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

function Matrix:inverse()
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

function Matrix:toInverse()
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
   +--------------------------------------------------+
   |                  Linear Systems                  |
   +--------------------------------------------------+
   |This section contains methods pertaining to       |
   |solving systems of linear equations. This includes|
   |linear solve and the LU factorization.            |
   +--------------------------------------------------+
]]

function Matrix:solve(vector)
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
      for j = i + 1, numberOfColumns, 1 do
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

function Matrix:LUDecomposition()
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

   return l, matrix, permutation
end

function Matrix:CholeskyDecomposition()
   local matrix = self:copy()
   local numberOfRows = #matrix
   local numberOfColumns = #matrix[1]

   if numberOfRows ~= numberOfColumns then
      error("Cannot compute Cholesky of rectangular matrix.")
   end

   for i = 1, numberOfColumns do
      for j = 1, i - 1 do
	 matrix[i][i] -= matrix[j][i] ^ 2
      end
      if matrix[i][i] <= 0 then
	 error("Matrix is not positive definite.")
      end
      matrix[i][i] = math.sqrt(matrix[i][i])
      for j = i + 1, numberOfColumns do
	 for k = 1, i - 1 do
	    matrix[i][j] -= matrix[k][i] * matrix[k][j]
	 end
	 matrix[i][j] /= matrix[i][i]
      end
   end

   for i = 1, numberOfColumns do
      for j = 1, i - 1 do
	 matrix[i][j] = 0
      end
   end

   return matrix
end

function Matrix:GramSchmidtQR(tolerance)
   tolerance = tolerance or 10 ^ -13
   local numberOfRows = #self
   local numberOfColumns = #self[1]
   local Q = Matrix:new(Vectors.gramSchmidt(self:transpose(), tolerance))
   local R = Q * self
   for i = 1, #R do
      for j = 1, math.min(i - 1, #R[1]) do
	 R[i][j] = 0
      end
   end
   Q = Q:transpose()
   return Q, R
end

function Matrix:HouseholderQR(tolerance: number)
   tolerance = tolerance or 10 ^ -13
   local a = self:copy()
   local g = {}
   for i = 1, a.width - 1 do
      local beta = 0
      for j = i, a.length do
	 beta = math.max(beta, math.abs(a[j][i]))
      end
      local sum = 0
      for j = i, a.length do
	 a[j][i] /= beta
	 sum = sum + math.pow(a[j][i], 2)
      end
      local tau = math.sqrt(sum)
      if a[i][i] < 0 then
	 tau = -tau
      end
      local eta = a[i][i] + tau
      a[i][i] = 1
      for j = i + 1, a.length do
	 a[j][i] /= eta
      end
      local gamma = eta / tau
      g[#g + 1] = gamma
      tau *= beta
      
      local everythingButIthCol = a:submatrix(i, a.length, i + 1, a.width)
      local ithCol = a:submatrix(i, a.length, i, i)
      
      a:setSubmatrix(i, a.length, i + 1, a.width,
		     everythingButIthCol -
		     ithCol:scaled(gamma) * ithCol:transpose() * everythingButIthCol)

      a[i][i] = -tau
   end

   return a, g
end

function Matrix.ExpandHouseholderQR(householderQR, gammaVector)
   local R = householderQR:copy()
   for i = 1, #R do
      for j = 1, math.min(i - 1, #R[1]) do
	 R[i][j] = 0
      end
   end
   local Q = Matrix:identity(#R, #R)
   for i = #R - 1, 1, -1 do
      local uPartial = householderQR:submatrix(i + 1, #R, i, i):transpose()[1]
      local u = {1}
      for i = 1, #uPartial do
	 u[i + 1] = uPartial[i]
      end
      u = Matrix:new(u)
      qSub = Q:submatrix(i, #R, i, #R)
      Q:setSubmatrix(i, #R, i, #R, qSub -
		     u:scaled(gammaVector[i]) * u:transpose() * qSub)
   end
   return Q, R
end

function Matrix:FullHouseholderQR()
   local a, g = self:HouseholderQR()
   return Matrix.ExpandHouseholderQR(a, g)
end

function Matrix:GivensQR(tol)
   tol = tol or 10 ^ -13
   local QT = Matrix:identity(#self)
   local R = self:copy()
   for i = 1, self.width do
      for j = self.length, i + 1, -1 do
	 local a = R[i][i]
	 local b = R[j][i]
	 if math.abs(b) < tol then
	    continue
	 end
	 local r = Numerics.hypot(a, b)
	 R[i][i], R[j][i] = r, 0
	 local c = a / r
	 local s = -b / r
	 for k = i + 1, self.width do
	    local x = R[i][k]
	    local y = R[j][k]
	    R[i][k] = c * x - s * y
	    R[j][k] = s * x + c * y
	 end
	 for k = 1, self.width do
	    local x = QT[i][k]
	    local y = QT[j][k]
	    QT[i][k] = c * x - s * y
	    QT[j][k] = s * x + c * y
	 end
      end
   end
   return QT:transpose(), R
end	 

function Matrix:newReflector(u: Vector)
   local data = {}
   local length = #u
   local norm = Tools.list.norm(u)
   local gamma = 2 / norm ^ 2
   for i = 1, length do
      data[i] = {}
      for j = 1, length do
	 if i == j then
	    data[i][j] = 1 - gamma * u[i] * u[j]
	 else
	    data[i][j] = -gamma * u[i] * u[j]
	 end
      end
   end
   return Matrix:new(data)
end

function Matrix.applyReflector(u: Vector, v: Vector)
   local dot = Vectors.dot(u, v)
   local norm = Vectors.norm(u)
   local gamma = 2 / norm ^ 2
   local newVector = Vectors.scale(gamma * dot, u)
   return Vectors.sub(v, newVector)
end

function Matrix.applyUnitReflector(u: Vector, v: Vector)
   local dot = Vectors.dot(u, v)
   local newVector = Vectors.scale(2 * dot, u)
   return Vectors.sub(v, newVector)
end

function Matrix:LinearLeastSquare(vector: Vector, tolerance: number)
   tolerance = tolerance or 10 ^ -13
   local a = self:copy()
   for i = 1, a.width do
      local beta = 0
      for j = i, a.length do
	 beta = math.max(beta, math.abs(a[j][i]))
      end
      local sum = 0
      for j = i, a.length do
	 a[j][i] /= beta
	 sum = sum + math.pow(a[j][i], 2)
      end
      local tau = math.sqrt(sum)
      if a[i][i] < 0 then
	 tau = -tau
      end
      local eta = a[i][i] + tau
      a[i][i] = 1
      for j = i + 1, a.length do
	 a[j][i] /= eta
      end
      local gamma = eta / tau
      tau *= beta

      local b = Matrix:zero(a.length, a.width)

      local temp = a:submatrix(i, a.length, i + 1, a.width)
      local temp2 = a:submatrix(i, a.length, i, i)
      b:setSubmatrix(i, a.length, 1, 1, (temp2:transpose() * temp):toScaled(-gamma):transpose())
      a:setSubmatrix(i, a.length, i + 1, a.width, temp + temp2 * b:submatrix(i, a.length, 1, 1):transpose())
      vector = Matrix:applyReflector(a:getColumnVector(i), vector)

      a[i][i] = -tau
   end
   return a:submatrix(1, math.min(self.length, self.width), 1, math.min(self.length, self.width)):solve(vector)
end

--[[
   +--------------------------------------------------+
   |                   Scalar Maps                    |
   +--------------------------------------------------+
   |This section contains many common scalar maps     |
   |including the trace and determinant.              |
   +--------------------------------------------------+
]]

function Matrix:determinant()
   if self.length == 2 and self.width == 2 then
      return self[1][1] * self[2][2] - self[1][2] * self[2][1]
   end
   local LUDecomposition
   local test = pcall(function()
	 LUDecomposition = self:LUDecomposition()
   end)
   if test == false then
      return 0
   end
   local determinant = 1
   for i = 1, self.length do
      determinant = determinant * LUDecomposition[2][i][i]
   end
   return determinant * math.pow(-1, Tools.combinatorics.inversionNumber(LUDecomposition[3]))
end

function Matrix:trace()
   local n = math.min(self.length, self.width)
   local sum = 0
   for i = 1, n do
      sum = sum + self[i][i]
   end
   return sum
end

function Matrix:dot(matrix)
   return (self:transpose() * matrix):trace()
end

function Matrix:frobenius()
   local sum = 0
   for i = 1, #self do
      for j = 1, #self[1] do
	 sum += self[i][j] ^ 2
      end
   end
   return math.sqrt(sum)
end

--[[
   +--------------------------------------------------+
   |             Eigenvalue Computations              |
   +--------------------------------------------------+
   |This section contains the math needed to compute  |
   |the eigenvalues of a matrix. The primary tool for |
   |this is the Francis algorithm. This has been      |
   |implemented for the Rayleigh shifts of degree 1,  |
   |2, 3, and 6. Much of the algorithms in this       |
   |section come from the book Fundamentals of Matrix |
   |Computation by Watkins.                           |
   +--------------------------------------------------+
]]

function Matrix:hessenbergForm()
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
	 if a[k + 1][k] < 0 then
	    tau = -tau
	 end
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

function Matrix:toHessenbergForm()
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
      local beta = 0
      for i = k + 1, n do
	 beta = math.max(beta, math.abs(a[i][k]))
      end
      local gamma = 0
      if beta ~= 0 then
	 local sum = 0
	 for i = k + 1, n do
	    a[i][k] = a[i][k] / beta
	    sum = sum + math.pow(a[i][k], 2)
	 end
	 local tau = math.sqrt(sum)
	 if a[k + 1][k] < 0 then
	    tau = -tau
	 end
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

function Matrix:francisOne(tol, maxIters)
   tol = tol or 10 ^ -7
   maxIters = maxIters or 1000
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }
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
	 return {
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2),
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2),
	 }
      else
	 return {
	    (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2,
	    (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2,
	 }
      end
   end
   local iterations = 0
   local temp = 0
   local tempMin = 0
   while iterations < maxIters do
      local min = 10 ^ 12
      local minimizer = 0
      for i = 1, n - 1 do
	 if math.abs(matrix[i + 1][i]) < min then
	    min = math.abs(matrix[i + 1][i])
	    minimizer = i
	 end
	 if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
	    return Tools.list.join(
	       matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters),
	       matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters)
	    )
	 end
      end
      local shift = matrix[n][n]
      local u = {
	 matrix[1][1] - shift + math.sqrt(math.pow(matrix[1][1] - shift, 2) + math.pow(matrix[2][1], 2)),
	 matrix[2][1],
      }
      local gamma = 2 / (math.pow(u[1], 2) + math.pow(u[2], 2))
      u = Matrix:new(u)
      local Q = Matrix:identity(2) - (u * u:transpose()):toScaled(gamma)
      for i = 1, n do
	 matrix:setSubmatrix(1, 2, i, i, Q * matrix:submatrix(1, 2, i, i))
      end
      for i = 1, n do
	 matrix:setSubmatrix(i, i, 1, 2, matrix:submatrix(i, i, 1, 2) * Q:transpose())
      end
      matrix:toHessenbergForm()
      iterations = iterations + 1
      temp = minimizer
      tempMin = min
      if iterations == maxIters then
	 --print("Francis One failed to converge in " .. tostring(maxIters) .. " iterations! Breaking on " .. tostring(min) .. ".")
	 return Tools.list.join(
	    matrix:submatrix(1, minimizer, 1, minimizer):francisOne(tol),
	    matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisOne(tol)
	 )
      end
   end
   return matrix
end

function Matrix:francisTwo(tol, maxIters)
   tol = tol or 10 ^ -7
   maxIters = maxIters or 10000
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }
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
	 return {
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2),
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2),
	 }
      else
	 return {
	    (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2,
	    (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2,
	 }
      end
   end
   local iterations = 0
   local temp = 0
   local tempMin = 0
   local p = 2
   while iterations < maxIters do
      local min = 10 ^ 12
      local minimizer = 0
      -- Check if the lower diagonal has any zero elements and split
      for i = 1, n - 1 do
	 if math.abs(matrix[i + 1][i]) < min then
	    min = math.abs(matrix[i + 1][i])
	    minimizer = i
	 end
	 if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
	    return Tools.list.join(
	       matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters),
	       matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters)
	    )
	 end
      end
      local shift = matrix:submatrix(n - 1, n, n - 1, n):francisOne()
      local u
      local gamma
      -- Compute the first column of shifted matrix
      if type(shift[1]) == "table" then
	 local a11 = Complex:new(matrix[1][1], 0)
	 local a12 = Complex:new(matrix[1][2], 0)
	 local a21 = Complex:new(matrix[2][1], 0)
	 local a22 = Complex:new(matrix[2][2], 0)
	 local rho1 = shift[1]
	 local rho2 = shift[2]
	 u = {
	    ((a11 - rho1) * (a11 - rho2) + a12 * a21)[1],
	    (a21 * ((a11 + a22) - (rho1 + rho2)))[1],
	    matrix[3][2] * matrix[2][1],
	 }
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
	 u = {
	    ((a11 - rho1) * (a11 - rho2) + a12 * a21),
	    (a21 * ((a11 + a22) - (rho1 + rho2))),
	    matrix[3][2] * matrix[2][1],
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
      matrix:toHessenbergForm()
      iterations = iterations + 1
      temp = minimizer
      tempMin = min
      if iterations == maxIters then
	 return Tools.list.join(
	    matrix:submatrix(1, minimizer, 1, minimizer):francisTwo(tol),
	    matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisTwo(tol)
	 )
      end
   end
   return matrix
end

function Matrix:francisThree(tol, maxIters)
   tol = tol or 10 ^ -10
   maxIters = maxIters or 10000
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }
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
	 return {
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2),
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2),
	 }
      else
	 return {
	    (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2,
	    (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2,
	 }
      end
   elseif n == 3 then
      return self:francisTwo(tol, maxIters)
   end
   local iterations = 0
   local temp = 0
   local tempMin = 0
   local p = 3
   while iterations < maxIters do
      local min = 10 ^ 12
      local minimizer = 0
      -- Check if the lower diagonal has any zero elements and split
      for i = 1, n - 1 do
	 if math.abs(matrix[i + 1][i]) < min then
	    min = math.abs(matrix[i + 1][i])
	    minimizer = i
	 end
	 if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
	    return Tools.list.join(
	       matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters),
	       matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters)
	    )
	 end
      end
      local shift = matrix:submatrix(n - 2, n, n - 2, n):francisOne(tol * 10, maxIters)
      Complex:sort(shift)
      local u
      local gamma
      -- Compute the first column of shifted matrix
      if type(shift[1]) == "table" or type(shift[2]) == "table" then
	 local a11 = Complex:new(matrix[1][1], 0)
	 local a12 = Complex:new(matrix[1][2], 0)
	 local a13 = Complex:new(matrix[1][3], 0)
	 local a21 = Complex:new(matrix[2][1], 0)
	 local a22 = Complex:new(matrix[2][2], 0)
	 local a23 = Complex:new(matrix[2][3], 0)
	 local a32 = Complex:new(matrix[3][2], 0)
	 local a33 = Complex:new(matrix[3][3], 0)
	 local rho1, rho2, rho3
	 if type(shift[1]) == "table" then
	    rho1 = shift[1]
	 else
	    rho1 = Complex:new(shift[1], 0)
	 end
	 if type(shift[2]) == "table" then
	    rho2 = shift[2]
	 else
	    rho2 = Complex:new(shift[2], 0)
	 end
	 if type(shift[3]) == "table" then
	    rho3 = shift[3]
	 else
	    rho3 = Complex:new(shift[3], 0)
	 end
	 u = {
	    (
	       a21 * (a13 * a32 + a12 * (a11 + a22 - rho2 - rho3))
	       + (a11 - rho1) * (a11 * a11 + a12 * a21 + rho2 * rho3 - a11 * (rho2 + rho3))
	    )[1],
	    (
	       a21
	       * (
		  a12 * a21
		  + a22 * a22
		  + a23 * a32
		  + (a11 - rho1) * (a11 + a22 - rho2 - rho3)
		  + rho2 * rho3
		  - a22 * (rho2 + rho3)
		 )
	    )[1],
	    (a21 * a32 * (a11 + a22 + a33 - rho1 - rho2 - rho3))[1],
	    matrix[2][1] * matrix[3][2] * matrix[4][3],
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
	    (
	       a21 * (a13 * a32 + a12 * (a11 + a22 - rho2 - rho3))
	       + (a11 - rho1) * (a11 ^ 2 + a12 * a21 + rho2 * rho3 - a11 * (rho2 + rho3))
	    ),
	    (
	       a21
	       * (
		  a12 * a21
		  + a22 ^ 2
		  + a23 * a32
		  + (a11 - rho1) * (a11 + a22 - rho2 - rho3)
		  + rho2 * rho3
		  - a22 * (rho2 + rho3)
		 )
	    ),
	    (a21 * a32 * (a11 + a22 + a33 - rho1 - rho2 - rho3)),
	    matrix[2][1] * matrix[3][2] * matrix[4][3],
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
      matrix:toHessenbergForm()
      iterations = iterations + 1
      temp = minimizer
      tempMin = min
      if iterations == maxIters then
	 print(
	    "Francis Three failed to converge in "
	    .. tostring(maxIters)
	    .. " iterations! Breaking on "
	    .. tostring(min)
	    .. "."
	 )
	 return Tools.list.join(
	    matrix:submatrix(1, minimizer, 1, minimizer):francisThree(tol),
	    matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisThree(tol)
	 )
      end
   end
   return matrix
end

function Matrix:francisSix(tol, maxIters)
   -- Set Defaults
   tol = tol or 10 ^ -12
   maxIters = maxIters or 10000
   -- Create Matrix
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }
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
	 return {
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, math.sqrt(-temp1) / 2),
	    Complex:new((matrix[1][1] + matrix[2][2]) / 2, -math.sqrt(-temp1) / 2),
	 }
      else
	 return {
	    (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2,
	    (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2,
	 }
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
      local min = 10 ^ 12
      local minimizer = 0
      -- Check if the lower diagonal has any zero elements and split
      for i = 1, n - 1 do
	 if math.abs(matrix[i + 1][i]) < min then
	    min = math.abs(matrix[i + 1][i])
	    minimizer = i
	 end
	 if math.abs(matrix[i + 1][i]) < tol or (temp == minimizer and math.abs(tempMin - min) < tol) then
	    return Tools.list.join(
	       matrix:submatrix(1, i, 1, i):francisSix(tol, maxIters),
	       matrix:submatrix(i + 1, n, i + 1, n):francisSix(tol, maxIters)
	    )
	 end
      end
      local shift = matrix:submatrix(n - p + 1, n, n - p + 1, n):francisThree(tol * 10, maxIters)
      Complex:sort(shift)
      local u
      local gamma
      -- Compute the first column of shifted matrix
      local a11, a12, a13, a14, a15, a16, a17 =
	 Complex:new(matrix[1][1], 0),
	 Complex:new(matrix[1][2], 0),
	 Complex:new(matrix[1][3], 0),
	 Complex:new(matrix[1][4], 0),
	 Complex:new(matrix[1][5], 0),
	 Complex:new(matrix[1][6], 0),
	 Complex:new(matrix[1][7], 0)
      local a21, a22, a23, a24, a25, a26, a27 =
	 Complex:new(matrix[2][1], 0),
	 Complex:new(matrix[2][2], 0),
	 Complex:new(matrix[2][3], 0),
	 Complex:new(matrix[2][4], 0),
	 Complex:new(matrix[2][5], 0),
	 Complex:new(matrix[2][6], 0),
	 Complex:new(matrix[2][7], 0)
      local a32, a33, a34, a35, a36, a37 =
	 Complex:new(matrix[3][2], 0),
	 Complex:new(matrix[3][3], 0),
	 Complex:new(matrix[3][4], 0),
	 Complex:new(matrix[3][5], 0),
	 Complex:new(matrix[3][6], 0),
	 Complex:new(matrix[3][7], 0)
      local a43, a44, a45, a46, a47 =
	 Complex:new(matrix[4][3], 0),
	 Complex:new(matrix[4][4], 0),
	 Complex:new(matrix[4][5], 0),
	 Complex:new(matrix[4][6], 0),
	 Complex:new(matrix[4][7], 0)
      local a54, a55, a56, a57 =
	 Complex:new(matrix[5][4], 0),
	 Complex:new(matrix[5][5], 0),
	 Complex:new(matrix[5][6], 0),
	 Complex:new(matrix[5][7], 0)
      local a65, a66, a67 = Complex:new(matrix[6][5], 0), Complex:new(matrix[6][6], 0), Complex:new(matrix[6][7], 0)
      local a76, a77 = Complex:new(matrix[7][6], 0), Complex:new(matrix[7][7], 0)
      local rho1, rho2, rho3, rho4, rho5, rho6 =
	 Complex:new(shift[1]),
	 Complex:new(shift[2]),
	 Complex:new(shift[3]),
	 Complex:new(shift[4]),
	 Complex:new(shift[5]),
	 Complex:new(shift[6])
      u = {
	 (
	    a21
	    * ((a22 - rho2) * ((a22 - rho3) * ((a22 - rho4) * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a32 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a12 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a32 * (a23 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a33 - rho4) * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a43 * (a12 * a24 + a13 * a34 + a15 * a54 + a14 * (a11 + a44 - rho5 - rho6)) + a13 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a12 * (a21 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a11 - rho4) * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6)))) + a32 * (a23 * ((a22 - rho4) * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a32 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a12 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + (a33 - rho3) * (a23 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a33 - rho4) * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a43 * (a12 * a24 + a13 * a34 + a15 * a54 + a14 * (a11 + a44 - rho5 - rho6)) + a13 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a43 * (a24 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a34 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + (a44 - rho4) * (a12 * a24 + a13 * a34 + a15 * a54 + a14 * (a11 + a44 - rho5 - rho6)) + a54 * (a12 * a25 + a13 * a35 + a14 * a45 + a16 * a65 + a15 * (a11 + a55 - rho5 - rho6)) + a14 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a13 * (a21 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a11 - rho4) * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6)))) + a12 * (a21 * ((a22 - rho4) * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a32 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a12 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + (a11 - rho3) * (a21 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a11 - rho4) * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6)))))
	    + (a11 - rho1)
	    * (a21 * ((a22 - rho3) * ((a22 - rho4) * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a32 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a12 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a32 * (a23 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a33 - rho4) * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a43 * (a12 * a24 + a13 * a34 + a15 * a54 + a14 * (a11 + a44 - rho5 - rho6)) + a13 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + a12 * (a21 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a11 - rho4) * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6)))) + (a11 - rho2) * (a21 * ((a22 - rho4) * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + a32 * (a12 * a23 + a14 * a43 + a13 * (a11 + a33 - rho5 - rho6)) + a12 * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6))) + (a11 - rho3) * (a21 * (a13 * a32 + a12 * (a11 + a22 - rho5 - rho6)) + (a11 - rho4) * (a11 * a11 + a12 * a21 + rho5 * rho6 - a11 * (rho5 + rho6)))))
	 )[1],
	 (
	    a21
	    * (
	       a12
	       * a21
	       * (a32 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a12 * a21 * (a11 + a22 - rho5 - rho6) + (a22 - rho4) * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6)) + (a11 - rho3) * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho4) * (a11 + a22 - rho5 - rho6) + rho5 * rho6 - a22 * (rho5 + rho6)))
	       + a32 * (a13 * a21 * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho4) * (a11 + a22 - rho5 - rho6) + rho5 * rho6 - a22 * (rho5 + rho6)) + (a33 - rho3) * ((a33 - rho4) * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a43 * (a14 * a21 + a23 * a34 + a25 * a54 + a24 * (a22 + a44 - rho5 - rho6)) + a13 * a21 * (a11 + a22 - rho5 - rho6) + a23 * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))) + a43 * (a34 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + (a44 - rho4) * (a14 * a21 + a23 * a34 + a25 * a54 + a24 * (a22 + a44 - rho5 - rho6)) + a54 * (a15 * a21 + a23 * a35 + a24 * a45 + a26 * a65 + a25 * (a22 + a55 - rho5 - rho6)) + a14 * a21 * (a11 + a22 - rho5 - rho6) + a24 * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))) + a23 * (a32 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a12 * a21 * (a11 + a22 - rho5 - rho6) + (a22 - rho4) * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))))
	       + (a22 - rho2) * (a12 * a21 * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho4) * (a11 + a22 - rho5 - rho6) + rho5 * rho6 - a22 * (rho5 + rho6)) + a32 * ((a33 - rho4) * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a43 * (a14 * a21 + a23 * a34 + a25 * a54 + a24 * (a22 + a44 - rho5 - rho6)) + a13 * a21 * (a11 + a22 - rho5 - rho6) + a23 * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))) + (a22 - rho3) * (a32 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a12 * a21 * (a11 + a22 - rho5 - rho6) + (a22 - rho4) * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))))
	       + (a11 - rho1)
	       * (a12 * a21 * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho4) * (a11 + a22 - rho5 - rho6) + rho5 * rho6 - a22 * (rho5 + rho6)) + a32 * ((a33 - rho4) * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a43 * (a14 * a21 + a23 * a34 + a25 * a54 + a24 * (a22 + a44 - rho5 - rho6)) + a13 * a21 * (a11 + a22 - rho5 - rho6) + a23 * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))) + (a22 - rho3) * (a32 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a12 * a21 * (a11 + a22 - rho5 - rho6) + (a22 - rho4) * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6))) + (a11 - rho2) * (a32 * (a13 * a21 + a24 * a43 + a23 * (a22 + a33 - rho5 - rho6)) + a12 * a21 * (a11 + a22 - rho5 - rho6) + (a22 - rho4) * (a12 * a21 + a22 * a22 + a23 * a32 + rho5 * rho6 - a22 * (rho5 + rho6)) + (a11 - rho3) * (a12 * a21 + a22 * a22 + a23 * a32 + (a11 - rho4) * (a11 + a22 - rho5 - rho6) + rho5 * rho6 - a22 * (rho5 + rho6))))
	      )
	 )[1],
	 (
	    a21
	    * a32
	    * (
	       a13 * a21 * a32 * (a11 + a22 + a33 - rho4 - rho5 - rho6)
	       + a23 * a32 * (a12 * a21 + a23 * a32 + a33 * a33 + a34 * a43 + (a22 - rho4) * (a22 + a33 - rho5 - rho6) + rho5 * rho6 - a33 * (rho5 + rho6))
	       + a12 * a21 * (a12 * a21 + a23 * a32 + a33 * a33 + a34 * a43 + (a22 - rho4) * (a22 + a33 - rho5 - rho6) + (a11 - rho3) * (a11 + a22 + a33 - rho4 - rho5 - rho6) + rho5 * rho6 - a33 * (rho5 + rho6))
	       + a43 * (a14 * a21 * a32 + (a44 - rho4) * (a24 * a32 + a35 * a54 + a34 * (a33 + a44 - rho5 - rho6)) + a54 * (a25 * a32 + a34 * a45 + a36 * a65 + a35 * (a33 + a55 - rho5 - rho6)) + a24 * a32 * (a22 + a33 - rho5 - rho6) + a34 * (a23 * a32 + a33 * a33 + a34 * a43 + rho5 * rho6 - a33 * (rho5 + rho6)))
	       + (a33 - rho3) * (a13 * a21 * a32 + a43 * (a24 * a32 + a35 * a54 + a34 * (a33 + a44 - rho5 - rho6)) + a23 * a32 * (a22 + a33 - rho5 - rho6) + (a33 - rho4) * (a23 * a32 + a33 * a33 + a34 * a43 + rho5 * rho6 - a33 * (rho5 + rho6)))
	       + (a22 - rho2) * (a13 * a21 * a32 + a43 * (a24 * a32 + a35 * a54 + a34 * (a33 + a44 - rho5 - rho6)) + a23 * a32 * (a22 + a33 - rho5 - rho6) + a12 * a21 * (a11 + a22 + a33 - rho4 - rho5 - rho6) + (a33 - rho4) * (a23 * a32 + a33 * a33 + a34 * a43 + rho5 * rho6 - a33 * (rho5 + rho6)) + (a22 - rho3) * (a12 * a21 + a23 * a32 + a33 * a33 + a34 * a43 + (a22 - rho4) * (a22 + a33 - rho5 - rho6) + rho5 * rho6 - a33 * (rho5 + rho6)))
	       + (a11 - rho1)
	       * (a13 * a21 * a32 + a43 * (a24 * a32 + a35 * a54 + a34 * (a33 + a44 - rho5 - rho6)) + a23 * a32 * (a22 + a33 - rho5 - rho6) + a12 * a21 * (a11 + a22 + a33 - rho4 - rho5 - rho6) + (a33 - rho4) * (a23 * a32 + a33 * a33 + a34 * a43 + rho5 * rho6 - a33 * (rho5 + rho6)) + (a22 - rho3) * (a12 * a21 + a23 * a32 + a33 * a33 + a34 * a43 + (a22 - rho4) * (a22 + a33 - rho5 - rho6) + rho5 * rho6 - a33 * (rho5 + rho6)) + (a11 - rho2) * (a12 * a21 + a23 * a32 + a33 * a33 + a34 * a43 + (a22 - rho4) * (a22 + a33 - rho5 - rho6) + (a11 - rho3) * (a11 + a22 + a33 - rho4 - rho5 - rho6) + rho5 * rho6 - a33 * (rho5 + rho6)))
	      )
	 )[1],
	 (
	    a21
	    * a32
	    * a43
	    * (
	       a13 * a21 * a32
	       + a24 * a32 * a43
	       + a54 * (a35 * a43 + a46 * a65 + a45 * (a44 + a55 - rho5 - rho6))
	       + a34 * a43 * (a33 + a44 - rho5 - rho6)
	       + a23 * a32 * (a22 + a33 + a44 - rho4 - rho5 - rho6)
	       + a12 * a21 * (a11 + a22 + a33 + a44 - rho3 - rho4 - rho5 - rho6)
	       + (a44 - rho4) * (a34 * a43 + a44 * a44 + a45 * a54 + rho5 * rho6 - a44 * (rho5 + rho6))
	       + (a33 - rho3) * (a23 * a32 + a34 * a43 + a44 * a44 + a45 * a54 + (a33 - rho4) * (a33 + a44 - rho5 - rho6) + rho5 * rho6 - a44 * (rho5 + rho6))
	       + (a22 - rho2) * (a12 * a21 + a23 * a32 + a34 * a43 + a44 * a44 + a45 * a54 + (a33 - rho4) * (a33 + a44 - rho5 - rho6) + (a22 - rho3) * (a22 + a33 + a44 - rho4 - rho5 - rho6) + rho5 * rho6 - a44 * (rho5 + rho6))
	       + (a11 - rho1)
	       * (a12 * a21 + a23 * a32 + a34 * a43 + a44 * a44 + a45 * a54 + (a33 - rho4) * (a33 + a44 - rho5 - rho6) + (a22 - rho3) * (a22 + a33 + a44 - rho4 - rho5 - rho6) + (a11 - rho2) * (a11 + a22 + a33 + a44 - rho3 - rho4 - rho5 - rho6) + rho5 * rho6 - a44 * (rho5 + rho6))
	      )
	 )[1],
	 (
	    a21
	    * a32
	    * a43
	    * a54
	    * (
	       a12 * a21
	       + a23 * a32
	       + a34 * a43
	       + a45 * a54
	       + a55 * a55
	       + a56 * a65
	       + (a44 - rho4) * (a44 + a55 - rho5 - rho6)
	       + (a33 - rho3) * (a33 + a44 + a55 - rho4 - rho5 - rho6)
	       + (a22 - rho2) * (a22 + a33 + a44 + a55 - rho3 - rho4 - rho5 - rho6)
	       + (a11 - rho1) * (a11 + a22 + a33 + a44 + a55 - rho2 - rho3 - rho4 - rho5 - rho6)
	       + rho5 * rho6
	       - a55 * (rho5 + rho6)
	      )
	 )[1],
	 (
	    a21
	    * a32
	    * a43
	    * a54
	    * a65
	    * (a11 + a22 + a33 + a44 + a55 + a66 - rho1 - rho2 - rho3 - rho4 - rho5 - rho6)
	 )[1],
	 (a21 * a32 * a43 * a54 * a65 * a76)[1],
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
	 print(
	    "Francis Six failed to converge in "
	    .. tostring(maxIters)
	    .. " iterations! Breaking on "
	    .. tostring(min)
	    .. "."
	 )
	 return Tools.list.join(
	    matrix:submatrix(1, minimizer, 1, minimizer):francisSix(tol),
	    matrix:submatrix(minimizer + 1, n, minimizer + 1, n):francisSix(tol)
	 )
      end
   end
   return matrix
end

function Matrix:eigenvalues()
   local matrix = self:hessenbergForm()
   local eigList = matrix:francisSix()
   Complex:sort(eigList, "d")
   return Tools.list.sublist(eigList, math.min(matrix.length, matrix.width))
end

Matrices.Matrix = Matrix

local ComplexMatrix = {
   length = 0,
   width = 0,
}

function ComplexMatrix:new(list, realQ)
   if type(list[1]) == "table" then
      if realQ == true then
	 for i, v in ipairs(list) do
	    for j = 1, #v do
	       list[i][j] = Complex:new(list[i][j])
	    end
	 end
	 setmetatable(list, self)
	 self.__index = self
	 list.length = #list
	 list.width = #list[1]
	 return list
      else
	 setmetatable(list, self)
	 self.__index = self
	 list.length = #list
	 list.width = #list[1]
	 return list
      end
   else
      local matrix = {}
      for key, value in ipairs(list) do
	 matrix[key] = { value }
      end
      setmetatable(matrix, self)
      self.__index = self
      matrix.length = #list
      matrix.width = 1
      return matrix
   end
end

--[[
   +--------------------------------------------------+
   |                 Matrix Utilities                 |
   +--------------------------------------------------+
   |This section contains many useful functions       |
   |including those to copy and manipulate matrices.  |
   |Functions of the form "to_____" or "set_____" will|
   |change the underlying matrix, while others will   |
   |return a shallow copy.                            |
   +--------------------------------------------------+
]]

function ComplexMatrix:copy()
   local data = {}
   for i = 1, self.length do
      data[i] = {}
      for j = 1, self.width do
	 data[i][j] = self[i][j]
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:submatrix(rowStart, rowEnd, columnStart, columnEnd)
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
   return ComplexMatrix:new(data)
end

function ComplexMatrix:setSubmatrix(rowStart, rowEnd, columnStart, columnEnd, matrix, realQ)
   local row, column = 1, 1
   for i = rowStart, rowEnd, 1 do
      column = 1
      for j = columnStart, columnEnd, 1 do
	 if realQ == true then
	    self[i][j] = Complex:new(matrix[row][column])
	 else
	    self[i][j] = matrix[row][column]
	 end
	 column = column + 1
      end
      row = row + 1
   end
   return self
end

function ComplexMatrix:toScaled(lambda)
   for i = 1, self.length do
      for j = 1, self.width do
	 self[i][j] = self[i][j]:scale(lambda)
      end
   end
   return self
end

function ComplexMatrix:scaled(lambda)
   local copy = self:copy()
   for i = 1, copy.length do
      for j = 1, copy.width do
	 copy[i][j] = copy[i][j]:scale(lambda)
      end
   end
   return copy
end

function ComplexMatrix:padTo(length, width)
   local data = {}
   for i = 1, length, 1 do
      data[i] = {}
      for j = 1, width, 1 do
	 data[i][j] = self[i][j] or Complex:new(0)
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:strassenSubdivide()
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
	 data1[i][j], data2[i][j], data3[i][j], data4[i][j] =
	    self[i][j] or Complex:new(0),
	    self[i][jPlus] or Complex:new(0),
	    self[iPlus][j] or Complex:new(0),
	    self[iPlus][jPlus] or Complex:new(0)
      end
   end
   return {
      ComplexMatrix:new(data1),
      ComplexMatrix:new(data2),
      ComplexMatrix:new(data3),
      ComplexMatrix:new(data4),
   }
end

function ComplexMatrix:column(i)
   if i > self.width then
      error("Matrix doesn't have " .. tostring(i) .. " columns.")
   end
   local column = {}
   for j = 1, self.width do
      column[j] = self[i][j]
   end
   return column
end

function _padStringToLength(string, length)
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

function ComplexMatrix:pretty(n, m)
   local length = self.length
   local width = self.width

   n = n or 20
   m = m or 20

   local result = ""

   if length == 1 then
      result = "("
      for i = 1, width - 1 do
	 result = result .. _padStringToLength(string.sub(tostring(self[1][i]), 1, n), 20) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self[1][width]), 20), 1, n) .. ")"
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
	 result = result .. _padStringToLength(string.sub(tostring(self[i][j]), 1, n), 20) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self[i][width]), 1, n), 20)
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
   +--------------------------------------------------+
   |                Matrix Metamethods                |
   +--------------------------------------------------+
   |This section contains all of the metamethods for  |
   |matrices. Addition and subtraction are relatively |
   |standard, but multiplication is an implementation |
   |of Strassen's method. The size of matrix at which |
   |Strassen multiplication will be used is set in    |
   |Matrices.constant.STRASSENLIMIT.                  |
   +--------------------------------------------------+
]]

function ComplexMatrix.__add(left, right)
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

   return ComplexMatrix:new(data)
end

function ComplexMatrix.__sub(left, right)
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

   return ComplexMatrix:new(data)
end

function ComplexMatrix.__mul(left, right)
   if left.width ~= right.length then
      error("Attempting to multiply matrices of incompatible dimension!", -1)
   end
   local data
   if math.max(left.length, left.width, right.length, right.width) >= Matrices.constants.COMPLEXSTRASSENLIMIT then
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
      local strassenResult = { M1 + M4 - M5 + M7, M3 + M5, M2 + M4, M1 - M2 + M3 + M6 }
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
      return ComplexMatrix:new(data):submatrix(1, left.length, 1, right.width)
   else
      data = {}
      for i = 1, left.length do
	 local row = {}
	 for j = 1, right.width do
	    local val = Complex:new(0)
	    for k = 1, left.width do
	       val = val + left[i][k] * right[k][j]
	    end
	    row[j] = val
	 end
	 data[i] = row
      end
      return ComplexMatrix:new(data)
   end
end

function ComplexMatrix.__tostring(matrix)
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

function ComplexMatrix.__eq(left, right)
   if left.length ~= right.length or left.width ~= right.width then
      return false
   else
      for i = 1, left.length do
	 for j = 1, left.width do
	    if left[i][j] ~= right[i][j] then
	       return false
	    end
	 end
      end
   end
   return true
end

--[[
   +--------------------------------------------------+
   |                 Common Matrices                  |
   +--------------------------------------------------+
   |This section contains methods for constructing    |
   |many common matrices such as the identity and zero|
   |matrices.                                         |
   +--------------------------------------------------+
]]

function ComplexMatrix:zero(n, m)
   m = m or n
   local data = {}
   for i = 1, n do
      data[i] = {}
      for j = 1, m do
	 data[i][j] = Complex:new(0)
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:random(n, m, a, b)
   m = m or n
   a, b = a or 0, b or 1
   local data = {}
   for i = 1, n do
      data[i] = {}
      for j = 1, m do
	 data[i][j] = Complex:new((b - a) * math.random() + a, (b - a) * math.random() + a)
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:identity(n)
   local data = {}
   for i = 1, n do
      data[i] = {}
      for j = 1, n do
	 if i == j then
	    data[i][j] = Complex:new(1)
	 else
	    data[i][j] = Complex:new(0)
	 end
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:permutation(permutation)
   local matrix = ComplexMatrix:identity(#permutation)
   local data = {}
   for key, value in ipairs(permutation) do
      data[key] = matrix[value]
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:permuted(permutation)
   local data = {}
   for key, value in ipairs(permutation) do
      data[key] = self[value]
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:toPermuted(permutation)
   local matrix = self:copy()
   for key, value in ipairs(permutation) do
      self[key] = matrix[value]
   end
   return self
end

--[[
   +--------------------------------------------------+
   |                   Matrix Maps                    |
   +--------------------------------------------------+
   |This section contains the common maps that send   |
   |matrices to matrices. This includes functions like|
   |the transpose and inverse.                        |
   +--------------------------------------------------+
]]

function ComplexMatrix:transpose()
   local data = {}
   for j = 1, self.width do
      data[j] = {}
      for i = 1, self.length do
	 data[j][i] = self[i][j]
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:conjugateTranspose()
   local data = {}
   for j = 1, self.width do
      data[j] = {}
      for i = 1, self.length do
	 data[j][i] = self[i][j]:conjugate()
      end
   end
   return ComplexMatrix:new(data)
end

function ComplexMatrix:toTranspose()
   local data = {}
   for j = 1, self.width do
      data[j] = {}
      for i = 1, self.length do
	 data[j][i] = self[i][j]
      end
   end
   self = data
   return ComplexMatrix:new(data)
end

function ComplexMatrix:toConjugateTranspose()
   local data = {}
   for j = 1, self.width do
      data[j] = {}
      for i = 1, self.length do
	 data[j][i] = self[i][j]:conjugate()
      end
   end
   self = data
   return ComplexMatrix:new(data)
end

function ComplexMatrix:inverse()
   local length, width = self.length, self.width
   local matrix = self:copy()

   if length ~= width then
      error("Cannot compute inverse of rectangular matrix.", -1)
   end

   local result = matrix:identity(length)

   for i = 1, width - 1 do
      local maxRow = i
      local max = matrix[i][i] or Complex:new(0)

      for j = i, length do
	 local maxCandidate = matrix[j][i]
	 maxCandidate = maxCandidate:norm()
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
	 local valOverMax = val:scale(1 / max)
	 matrix[j][i] = Complex:new(0)
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

function ComplexMatrix:toInverse()
   local length, width = self.length, self.width
   local matrix = self:copy()

   if length ~= width then
      error("Cannot compute inverse of rectangular matrix.", -1)
   end

   local result = matrix:identity(length)

   for i = 1, width - 1 do
      local maxRow = i
      local max = matrix[i][i] or Complex:new(0)

      for j = i, length do
	 local maxCandidate = matrix[j][i]
	 maxCandidate = maxCandidate:norm()
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
	 local valOverMax = val:scale(1 / max)
	 matrix[j][i] = Complex:new(0)
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
   +--------------------------------------------------+
   |                  Linear Systems                  |
   +--------------------------------------------------+
   |This section contains methods pertaining to       |
   |solving systems of linear equations. This includes|
   |linear solve and the LU factorization.            |
   +--------------------------------------------------+
]]

function ComplexMatrix:solve(vector)
   local matrix = self:copy()
   local numberOfRows = #matrix
   local numberOfColumns = #matrix[1]

   if numberOfRows ~= numberOfColumns then
      error("Cannot solve rectangular system with this function.")
   end

   local columnVector = ComplexMatrix:new(vector, true)

   for i = 1, numberOfColumns - 1 do
      local maxRow = i
      local max = matrix[i][i]

      for j = i, numberOfRows do
	 local maxCandidate = matrix[j][i]
	 maxCandidate = maxCandidate:norm()
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
	 local valOverMax = val:scale(1 / max)
	 local columnVal1, columnVal2 = columnVector[j][1], columnVector[i][1]
	 columnVector[j][1] = columnVal1 - valOverMax * columnVal2
	 matrix[j][i] = Complex:new(0)
	 for k = i + 1, numberOfColumns do
	    matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
	 end
      end
   end

   local result = {}

   for i = numberOfRows, 1, -1 do
      local temp = Complex:new(0)
      for j = i + 1, numberOfColumns, 1 do
	 temp = temp + matrix[i][j] * columnVector[j][1]
      end
      columnVector[i][1] = columnVector[i][1]
      if matrix[i][i]:norm() == 0 then
	 error("Matrix system is not solvable")
      end
      columnVector[i][1] = (columnVector[i][1] - temp) / matrix[i][i]
   end

   for i = 1, numberOfRows do
      result[i] = columnVector[i][1]
   end

   return result
end

function ComplexMatrix:LUDecomposition()
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

   local l = ComplexMatrix:identity(numberOfRows)

   for i = 1, numberOfColumns - 1 do
      local maxRow = i
      local max = matrix[i][i]

      for j = i, numberOfRows do
	 local maxCandidate = matrix[j][i]
	 maxCandidate = maxCandidate:norm()
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
	 local valOverMax = val:scale(1 / max)
	 l[j][i] = valOverMax
	 matrix[j][i] = Complex:new(0)
	 for k = i + 1, numberOfColumns do
	    matrix[j][k] = matrix[j][k] - matrix[i][k] * valOverMax
	 end
      end
   end

   return { l, matrix, permutation }
end

--[[
   +--------------------------------------------------+
   |                   Scalar Maps                    |
   +--------------------------------------------------+
   |This section contains many common scalar maps     |
   |including the trace and determinant.              |
   +--------------------------------------------------+
]]

function ComplexMatrix:determinant()
   if self.length == 2 and self.width == 2 then
      return self[1][1] * self[2][2] - self[1][2] * self[2][1]
   end
   local LUDecomposition
   local test = pcall(function()
	 LUDecomposition = self:LUDecomposition()
   end)
   if test == false then
      return Complex:new(0)
   end
   local determinant = 1
   for i = 1, self.length do
      determinant = determinant * LUDecomposition[2][i][i]
   end
   return determinant * math.pow(-1, Tools.combinatorics.inversionNumber(LUDecomposition[3]))
end

function ComplexMatrix:trace()
   local n = math.min(self.length, self.width)
   local sum = Complex:new(0)
   for i = 1, n do
      sum = sum + self[i][i]
   end
   return sum
end

function ComplexMatrix:dot(matrix)
   return (self:conjugateTranspose() * matrix):trace()
end

--[[
   +--------------------------------------------------+
   |             Eigenvalue Computations              |
   +--------------------------------------------------+
   |This section contains the math needed to compute  |
   |the eigenvalues of a matrix. The primary tool for |
   |this is the Francis algorithm. This has been      |
   |implemented for the Rayleigh shifts of degree 1,  |
   |2, 3, and 6. Much of the algorithms in this       |
   |section come from the book Fundamentals of Matrix |
   |Computation by Watkins.                           |
   +--------------------------------------------------+
]]

function ComplexMatrix:hessenbergForm()
   local a
   if self.length ~= self.width then
      a = self:padTo(math.max(self.length, self.width))
   else
      a = self:copy()
   end
   local n = a.length
   local b = ComplexMatrix:zero(n, 1)
   for k = 1, n - 2 do
      local gamma = 0
      local sum = 0
      for i = k + 1, n do
	 sum = sum + a[i][k]:norm() ^ 2
      end
      local tau = math.sqrt(sum)
      a[k + 1][k] = a[k + 1][k] - Complex:new(tau)
      local sum = 0
      for i = k + 1, n do
	 sum = sum + a[i][k]:norm() ^ 2
      end
      gamma = 2 / sum

      print(tau)

      local temp = a:submatrix(k + 1, n, k + 1, n)
      local temp2 = a:submatrix(k + 1, n, k, k)
      a:setSubmatrix(k + 1, n, k + 1, n, temp + temp2 * (temp2:conjugateTranspose() * temp):toScaled(-gamma))

      temp = a:submatrix(1, n, k + 1, n)
      temp2 = a:submatrix(k + 1, n, k, k)
      b:setSubmatrix(1, n, 1, 1, (temp * temp2):toScaled(-gamma))
      a:setSubmatrix(1, n, k + 1, n, temp + b:submatrix(1, n, 1, 1) * temp2:conjugateTranspose())
      a:setSubmatrix(k + 1, n, k, k, temp2:toScaled(0))

      a[k + 1][k] = Complex:new(tau)
   end
   return a
end

function ComplexMatrix:toHessenbergForm()
   local a
   if self.length ~= self.width then
      a = self:padTo(math.max(self.length, self.width))
      self = a
   else
      a = self
   end
   local n = a.length
   local b = ComplexMatrix:zero(n, 1)
   for k = 1, n - 2 do
      local maxList = {}
      for i = k + 1, n do
	 maxList[#maxList + 1] = a[i][k]:norm()
      end
      local beta = math.max(table.unpack(maxList))
      local gamma = 0
      if beta ~= 0 then
	 local sum = Complex:new(0)
	 for i = k + 1, n do
	    a[i][k] = a[i][k]:scale(1 / beta)
	    sum = sum + a[i][k]:norm() ^ 2
	 end
	 local tau = math.sqrt(sum)
	 local eta = a[k + 1][k]:norm() + tau
	 a[k + 1][k] = Complex:new(1)
	 for i = k + 2, n do
	    a[i][k] = a[i][k]:scale(1 / eta)
	 end
	 gamma = eta / tau
	 tau = tau * beta

	 local temp = a:submatrix(k + 1, n, k + 1, n)
	 local temp2 = a:submatrix(k + 1, n, k, k)
	 b:setSubmatrix(k + 1, n, 1, 1, (temp2:conjugateTranspose() * temp):toScaled(-gamma):conjugateTranspose())
	 a:setSubmatrix(k + 1, n, k + 1, n, temp + temp2 * b:submatrix(k + 1, n, 1, 1):conjugateTranspose())

	 temp = a:submatrix(1, n, k + 1, n)
	 temp2 = a:submatrix(k + 1, n, k, k)
	 b:setSubmatrix(1, n, 1, 1, (temp * temp2):toScaled(-gamma))
	 a:setSubmatrix(1, n, k + 1, n, temp + b:submatrix(1, n, 1, 1) * temp2:conjugateTranspose())
	 a:setSubmatrix(k + 1, n, k, k, temp2:toScaled(0))

	 a[k + 1][k] = Complex:new(-tau)
      end
   end
   return a
end

Matrices.ComplexMatrix = ComplexMatrix

--[[
   +--------------------------------------------------+
   |                 Sparse Matrices                  |
   +--------------------------------------------------+
   |Sparse matrices are the go-to option if your      |
   |problem involves a matrix with a lot of zero      |
   |entries. Sparse matrices are significantly faster |
   |if the number of non-zero entries is low, but are,|
   |in general, slower to use than regular matrices.  |
   |In terms of memory, sparse matrices are always    |
   |cheaper to store than their traditional           |
   |counterparts, so it may be worth the hit in CPU   |
   |performance.                                      |
   +--------------------------------------------------+
   |                 !!!WARNING!!!                    |
   +--------------------------------------------------+
   |While Matrix can be treated roughly like a table  |
   |of tables in terms of data access, SparseMatrix   |
   |uses a different data structure to store the      |
   |values. As such, users should use SparseMatrix:set|
   |and SparseMatrix:get to get and set values.       |
   +--------------------------------------------------+
]]

local SparseMatrix = {
   length = 0,
   rowConstant = 0,
   width = 0,
   data = nil,
}

function SparseMatrix:set(i: number, j: number, val: number): Object
   if val == 0 then
      self.data[self.width * (i - 1) + j] = nil
   else
      self.data[self.width * (i - 1) + j] = val
   end
   return self
end

function SparseMatrix:tableSet(table: Tensor): Object
   for i, _ in pairs(table) do
      for j, w in pairs(table[i]) do
	 if w == 0 then
	    self.data[self.width * (i - 1) + j] = nil
	 else
	    self.data[self.width * (i - 1) + j] = w
	 end
      end
   end
   return self
end

function SparseMatrix:rowSet(i: number, table: Vector): Object
   local rowConstant = self.width * (i - 1)
   for j, w in pairs(table) do
      if w == 0 then
	 self.data[rowConstant + j] = nil
      else
	 self.data[rowConstant + j] = w
      end
   end
end

function SparseMatrix:columnSet(i: number, table: Vector): Object
   for j, w in pairs(table) do
      if w == 0 then
	 self.data[self.width * (j - 1) + i] = nil
      else
	 self.data[self.width * (j - 1) + i] = w
      end
   end
   return self
end

function SparseMatrix:get(i: number, j: number): number
   return self.data[self.width * (i - 1) + j] or 0
end

function SparseMatrix:new(list: Tensor, length: number, width: number): Object
   local sparseMatrix = {}
   
   if type(list[1]) == "table" then
      setmetatable(sparseMatrix, self)
      self.__index = self
      sparseMatrix.length = #list
      sparseMatrix.width = width or #list[1]
      sparseMatrix.data = {}
      sparseMatrix:tableSet(list)
      return sparseMatrix
   elseif length == nil and width == nil then
      setmetatable(sparseMatrix, self)
      self.__index = self
      sparseMatrix.length = #list
      sparseMatrix.width = 1
      sparseMatrix.data = {}
      sparseMatrix:columnSet(list)
      return sparseMatrix
   elseif width == nil then
      setmetatable(sparseMatrix, self)
      self.__index = self
      sparseMatrix.length = length
      sparseMatrix.width = length
      sparseMatrix.data = {}
      sparseMatrix:tableSet(list)
      return sparseMatrix
   else
      setmetatable(sparseMatrix, self)
      self.__index = self
      sparseMatrix.length = length
      sparseMatrix.width = width
      sparseMatrix.data = {}
      for i, v in pairs(list) do
	 sparseMatrix.data[i] = v
      end
      return sparseMatrix
   end
end

function SparseMatrix:clean(): Object
   for i, v in pairs(self.data) do
      if v == 0 then
	 self.data[i] = nil
      end
   end
   return self
end

function SparseMatrix:numericalClean(tol: number): Object
   for i, v in pairs(self.data) do
      if math.abs(v) < tol then
	 self.data[i] = nil
      end
   end
   return self
end

function SparseMatrix:copy(): Object
   local data = {}
   for i, v in pairs(self.data) do
      data[i] = v
   end
   return SparseMatrix:new(data, self.length, self.width)
end

function SparseMatrix:addBand(value: number, position: number)
   local cp = self:copy()
   if position == 0 then
      for i = 1, math.min(self.length, self.width) do
	 cp.data[self.width * (i - 1) + i] = value
      end
   elseif position > 0 then
      for i = 1, math.min(self.length, self.width - position) do
	 cp.data[self.width * (i - 1) + i + position] = value
      end
   else
      for i = 1, math.min(self.length + position, self.width) do
	 cp.data[self.width * (i - 1 - position) + i] = value
      end
   end
   return cp
end

function SparseMatrix:toAddedBand(value: number, position: number)
   if position == 0 then
      for i = 1, math.min(self.length, self.width) do
	 local index = self.width * (i - 1) + i
	 if self.data[index] then
	    self.data[index] += value
	 else
	    self.data[index] = value
	 end
      end
   elseif position > 0 then
      for i = 1, math.min(self.length, self.width - position) do
	 local index = self.width * (i - 1) + i + position
	 if self.data[index] then
	    self.data[index] += value
	 else
	    self.data[index] = value
	 end
      end
   else
      for i = 1, math.min(self.length + position, self.width) do
	 local index = self.width * (i - 1 - position) + i
	 if self.data[index] then
	    self.data[index] += value
	 else
	    self.data[index] = value
	 end
      end
   end
   return self
end

function SparseMatrix:sparsity(): number
   local denominator = self.length * self.width
   local numerator = 0
   for _, _ in pairs(self.data) do
      numerator += 1
   end
   return numerator / denominator
end

function SparseMatrix:transpose()
   local data = {}
   for i = 1, self.length do
      for j = 1, self.width do
	 data[(j - 1) * self.length + i] = self.data[(i - 1) * self.width + j]
      end
   end
   return SparseMatrix:new(data, self.width, self.length)
end

function SparseMatrix:getRow(n: number): Vector
   local data = {}
   for i = 1, self.width do
      data[i] = self:get(n, i)
   end
   return data
end

function SparseMatrix:getColumn(n: number): Vector
   local data = {}
   for i = 1, self.length do
      data[i] = self:get(i, n)
   end
   return data
end

function SparseMatrix:getDiagonal()
   local data = {}
   for i = 1, math.min(self.length, self.width) do
      data[i] = self:get(i, i)
   end
   return data
end

function SparseMatrix:getSubdiagonal()
   local data = {}
   for i = 1, math.min(self.length - 1, self.width - 1) do
      data[i] = self:get(i + 1, i)
   end
   return data
end

function SparseMatrix:toDense()
   local data  = {}
   for i = 1, self.length do
      data[i] = {}
      for j = 1, self.width do
	 data[i][j] = self:get(i, j)
      end
   end
   
   return Matrix:new(data)
end
      
function SparseMatrix.__tostring(matrix: Object): string
   local result = "{"

   for i = 1, matrix.length - 1 do
      result = result .. "{"
      for j = 1, matrix.width - 1 do
	 result = result .. tostring(matrix:get(i, j)) .. ", "
      end
      result = result .. tostring(matrix:get(i, matrix.width)) .. "}, "
   end

   result = result .. "{"
   for j = 1, matrix.width - 1 do
      result = result .. tostring(matrix:get(matrix.length, j)) .. ", "
   end
   result = result .. tostring(matrix:get(matrix.length, matrix.width)) .. "}"

   result = result .. "}"

   return result
end

function SparseMatrix:pretty(n: number, m: number): string
   local length = self.length
   local width = self.width

   n = n or 20
   m = m or 20

   local result = ""

   if length == 1 then
      result = "("
      for i = 1, width - 1 do
	 result = result .. _padStringToLength(string.sub(tostring(self.get(1, i)), 1, n), 20) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self.get(1, width)), 1, n), 20) .. ")"
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
	 result = result .. _padStringToLength(string.sub(tostring(self.get(i, j)), 1, n), 20) .. " "
      end
      result = result .. _padStringToLength(string.sub(tostring(self.get(i, width)), 1, n), 20)
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

function SparseMatrix:padTo(length: number, width: number): Object
   local copy = SparseMatrix:copy()
   copy.rowConstant = length - 1
   copy.length = length
   copy.width = width
   return copy
end

function SparseMatrix:strassenSubdivide(): Array<Object>
   local size = math.max(self.length, self.width)
   if size % 2 == 1 then
      size = size + 1
   end
   size = size / 2
   local data1 = SparseMatrix:new({}, size, size)
   local data2 = SparseMatrix:new({}, size, size)
   local data3 = SparseMatrix:new({}, size, size)
   local data4 = SparseMatrix:new({}, size, size)
   for i = 1, size, 1 do
      local iPlus = i + size
      for j = 1, size, 1 do
	 local jPlus = j + size
	 data1:set(i, j, self:get(i, j))
	 data2:set(i, j, self:get(i, jPlus))
	 data3:set(i, j, self:get(iPlus, j))
	 data4:set(i, j, self:get(iPlus, jPlus))
      end
   end
   return {
      data1,
      data2,
      data3,
      data4,
   }
end

function SparseMatrix.__add(left: Object, right: Object): Object
   if left.length ~= right.length or left.width ~= right.width then
      error("Attempting to add matrices of incompatible dimension!", -1)
   end

   local data = SparseMatrix:new({}, left.length, left.width)
   for i, v in pairs(left.data) do
      data.data[i] = v
   end
   for i, v in pairs(right.data) do
      local val = (data.data[i] or 0) + v
      if val == 0 then
	 data.data[i] = nil
      else
	 data.data[i] = val
      end
   end

   return data
end

function SparseMatrix.__sub(left: Object, right: Object): Object
   if left.length ~= right.length or left.width ~= right.width then
      error("Attempting to add matrices of incompatible dimension!", -1)
   end

   local data = SparseMatrix:new({}, left.length, left.width)
   for i, v in pairs(left.data) do
      data.data[i] = v
   end
   for i, v in pairs(right.data) do
      local val = (data.data[i] or 0) - v
      if val == 0 then
	 data.data[i] = nil
      else
	 data.data[i] = val
      end
   end

   return data
end

function SparseMatrix.__mul(left: Object, right: Object): Object
   if left.width ~= right.length then
      error("Attempting to multiply matrices of incompatible dimension!", -1)
   end
   local data = {}
   local leftData = left.data
   local rightData = right.data
   local leftLength = left.length
   local leftWidth = left.width
   local rightLength = right.length
   local rightWidth = right.width
   for k, v in pairs(leftData) do
      local c = k % leftWidth
      if c == 0 then
	 c = leftWidth
      end
      local r = math.floor((k - c) / leftWidth)
      for kk, vv in pairs(rightData) do
	 local cc = kk % rightWidth
	 if cc == 0 then
	    cc = rightWidth
	 end
	 local rr = math.floor((kk - cc) / rightWidth)
	 if c == rr + 1 then
	    local rcc = cc + r * rightWidth
	    local temp = data[rcc]
	    local vvv = v * vv
	    if temp ~= nil and vvv ~= 0 then
	       data[rcc] = temp + vvv
	    elseif vvv ~= 0 then
	       data[rcc] = vvv
	    end
	 end
      end
   end
   return SparseMatrix:new(data, leftLength, rightWidth)
end

function SparseMatrix.__eq(left: Object, right: Object): boolean
   if left.length ~= right.length or left.width ~= right.width then
      return false
   else
      local test = left - right
      for _, v in pairs(test.data) do
	 return false
      end
   end
   return true
end

function SparseMatrix:zero(n: number, m: number): Object
   return SparseMatrix:new({}, n, m)
end

function SparseMatrix:identity(n: number): Object
   local matrix = SparseMatrix:new({}, n)
   for i = 1, n do
      matrix:set(i, i, 1)
   end
   return matrix
end

function SparseMatrix:scaledIdentity(n, mu)
   local matrix = SparseMatrix:new({}, n)
   for i = 1, n do
      matrix:set(i, i, mu)
   end
   return matrix
end

function SparseMatrix:conjugate(matrix)
   if getmetatable(matrix) == getmetatable(self) then
      return matrix:transpose() * self * matrix
   else
      local newmatrix = SparseMatrix:new(matrix)
      return newmatrix:transpose() * self * newmatrix
   end
end

function SparseMatrix:scale(mu)
   for i, v in self.data do
      self[i] = mu * v
   end
end

function SparseMatrix:scaled(mu)
   local matrix = self:copy()
   return matrix:scale(mu)
end

function SparseMatrix:rayleighQuotient(vector)
   local vector2 = self:apply(vector)
   local numerator = Vectors.dot(vector, vector2)
   local denominator = Vectors.dot(vector, vector)
   return numerator / denominator
end

function SparseMatrix:frobenius()
   local sum = 0
   for k, v in self.data do
      sum += v ^ 2
   end
   return math.sqrt(sum)
end

function SparseMatrix:submatrix(rowStart: number, rowEnd: number, columnStart: number, columnEnd: number): Object
   local data = {}
   for k, v in pairs(self.data) do
      local c = k % self.width
      if c == 0 then
	 c = self.width
      end
      local r = (k - c) / self.width + 1
      if rowStart <= r and r <= rowEnd and columnStart <= c and c <= columnEnd then
	 local rc = columnEnd - columnStart + 1
	 local nr = r - rowStart
	 local nc = c - columnStart + 1
	 rc = rc * nr + nc
	 data[rc] = v
      end
   end
   return SparseMatrix:new(data, rowEnd - rowStart + 1, columnEnd - columnStart + 1)
end

function SparseMatrix:setSubmatrix(rowStart: number, columnStart: number, matrix: Object): Object
   for k, v in pairs(matrix.data) do
      local c = k % matrix.width
      if c == 0 then
	 c = matrix.width
      end
      local r = (k - c) / matrix.width + 1
      if rowStart <= r and columnStart <= c then
	 local nr = r + rowStart - 2
	 local nc = c + columnStart - 1
	 local rc = self.width * nr + nc
	 self.data[rc] = v
      end
   end
   return self
end

function SparseMatrix:apply(vector: Vector): Vector
   local data = {}
   for i = 1, self.length do
      data[i] = 0
   end
   for k, v in pairs(self.data) do
      local c = k % self.width
      if c == 0 then
	 c = self.width
      end
      local r = (k - c) / self.width + 1
      if vector[c] then
	 data[r] = data[r] + v * vector[c]
      end
   end

   return data
end

function SparseMatrix:cascadeApply(matrixList: Array<Object>, vector: Vector): Vector
   local data = vector

   for i = #matrixList, 1, -1 do
      data = matrixList[i]:apply(data)
   end

   return data
end

--[[
   We are trying to find solutions to the problem Ax = lx, for l a scalar. The power method
   allows us to find the largest l such that there is a solution, and it also finds an 
   eigenvector.
]]
function SparseMatrix:powerMethod(vector: Vector, tolerance: number): (number, Vector)
   -- Choosing and initial vector for the power method. We want this function to be deterministic, but also don't want degenerate behavior.
   tolerance = tolerance or 10^(-13)
   if not vector then
      vector = {}
      local iter = 1
      for k, v in pairs(self.data) do
	 if iter > self.width then
	    break
	 else
	    vector[iter] = v
	    iter += 1
	 end
      end
      if iter == 1 then
	 for i = 1, self.width do
	    vector[i] = 0
	 end
	 return 0, vector
      else
	 for i = iter, self.width do
	    vector[i] = 0
	 end
      end
   end
   -- Repeatedly apply the matrix to the vector normalizing at each step until the change is minimal
   local newVector = Vectors.scale(1 / Vectors.norm(vector), self:apply(vector))
   local displacement = Vectors.sub(vector, newVector)
   while Vectors.norm(displacement) > tolerance do
      vector = newVector
      newVector = Vectors.scale(1 / Vectors.norm(vector), self:apply(vector))
      displacement = Vectors.sub(vector, newVector)
   end
   -- Right now, we have a vector that solves Ax = lx, but we don't know what l is!
   vector = newVector
   newVector = self:apply(vector)
   local sum = 0
   for i = 1, #vector do
      if(vector[i] ~= 0) then
	 sum += newVector[i] / vector[i]
      end
   end

   return (sum / #vector), vector
end

function SparseMatrix:inverseIteration(vector, mu, tolerance)
   tolerance = tolerance or 10 ^ -13
   mu = mu or 0
   if not vector then
      vector = {}
      local iter = 1
      for k, v in pairs(self.data) do
	 if iter > self.width then
	    break
	 else
	    vector[iter] = v
	    iter += 1
	 end
      end
      if iter == 1 then
	 for i = 1, self.width do
	    vector[i] = 0
	 end
	 return 0, vector
      else
	 for i = iter, self.width do
	    vector[i] = 0
	 end
      end
   end
   local matrix = self - SparseMatrix:scaledIdentity(self.length, mu)
   -- Repeatedly apply the matrix to the vector normalizing at each step until the change is minimal
   local newVector = Vectors.scale(1 / Vectors.norm(vector), matrix:solve(vector,
                                                                          tolerance))
   local displacement = Vectors.sub(vector, newVector)
   while Vectors.norm(displacement) > tolerance do
      vector = newVector
      newVector = Vectors.scale(1 / Vectors.norm(vector), matrix:solve(vector))
      displacement = Vectors.sub(vector, newVector)
   end
   -- Right now, we have a vector that solves Ax = lx, but we don't know what l is!
   vector = newVector
   newVector = self:apply(vector)
   local sum = 0
   for i = 1, #vector do
      if(vector[i] ~= 0) then
	 sum += newVector[i] / vector[i]
      end
   end

   return (sum / #vector), vector
end

function SparseMatrix:rayleighIteration(vector, mu, tolerance)
   tolerance = tolerance or 10 ^ -13
   -- Make a vector if one was not provided
   if not vector then
      vector = {}
      local iter = 1
      for k, v in pairs(self.data) do
	 if iter > self.width then
	    break
	 else
	    vector[iter] = v
	    iter += 1
	 end
      end
      if iter == 1 then
	 for i = 1, self.width do
	    vector[i] = 0
	 end
	 return 0, vector
      else
	 for i = iter, self.width do
	    vector[i] = 0
	 end
      end
   end
   mu = mu or self:rayleighQuotient(vector)
   local matrix = self - SparseMatrix:scaledIdentity(self.length, mu)
   -- Repeatedly apply the matrix to the vector normalizing at each step until the change is minimal
   local newVector = matrix:solve(vector)
   newVector = Vectors.scale(1 / Vectors.norm(newVector), newVector)
   local displacement = Vectors.sub(vector, newVector)
   while Vectors.norm(displacement) > tolerance do
      vector = newVector
      mu = self:rayleighQuotient(vector)
      matrix = self - SparseMatrix:scaledIdentity(self.length, mu)
      if Vectors.norm(matrix:apply(vector)) < tolerance then
	 break
      end
      newVector = matrix:solve(vector, tolerance)
      newVector = Vectors.scale(1 / Vectors.norm(newVector), newVector)
      displacement = Vectors.sub(vector, newVector)
   end
   -- Right now, we have a vector that solves Ax = lx, but we don't know what l is!
   vector = newVector
   newVector = self:apply(vector)
   local sum = 0
   for i = 1, #vector do
      if(vector[i] ~= 0) then
	 sum += newVector[i] / vector[i]
      end
   end

   return (sum / #vector), vector
end

function SparseMatrix:arnoldiProcess(n: number, x: Vector, tolerance: number): Array<Vector>
   tolerance = tolerance or 10 ^ -13
   n = n or math.min(self.length, self.width) + 1
   local t = 0
   local q = {}
   local h = {}
   q[1] = Vectors.scale(1 / Vectors.norm(x), x) or Vectors.randomVector(self.width, 2)
   for i = 2, n do
      t += 1
      q[i] = self:apply(q[i - 1])
      for j = 1, i - 1 do
	 if not h[j] then
	    h[j] = {}
	 end
	 local dot = Vectors.dot(q[j], q[i])
	 if math.abs(dot) > tolerance then
	    h[j][i - 1] = dot
	 end
	 q[i] = Vectors.sub(q[i], Vectors.scale(dot, q[j]))
      end
      if not h[i] then
	 h[i] = {}
      end
      local norm = Vectors.norm(q[i])
      if math.abs(norm) < tolerance then
	 q[i] = nil
	 break
      end
      h[i][i - 1] = norm
      q[i] = Vectors.scale(1 / norm, q[i])
   end
   return q, SparseMatrix:new(h, t + 1, t)
end

-- Our goal is to solve the equation Ax = lx for A our matrix, x some unknown
--   vector and l some real/complex number. This is called the eigenvalue problem.
--   With this function we answer the question "what are the values of l for which
--   an answer exists?"
function SparseMatrix:eigenvalues(n: number, x: Vector, tolerance: number)
   tolerance = tolerance or 10 ^ -13
   n = n or 1
   if n == 1 then
      local temp, tempVec = self:powerMethod(x, tolerance)
      return {temp}
   end
   local m = math.min(self.length, 2 * n)
   x = x or Vectors.randomVector(self.width, 2)
   local t
   local q = {}
   local oldEigenvalues = {tolerance + 1}
   local eigenvalues = {0}
   local Q, H
   q[1] = Vectors.scale(1 / Vectors.norm(x), x)
   while Vectors.norm(Vectors.sub(oldEigenvalues, eigenvalues)) > tolerance do
      t=0
      for i = 2, m do
	 t += 1
	 q[i] = self:apply(q[i - 1])
	 for j = 1, i - 1 do
	    local dot = Vectors.dot(q[j], q[i])
	    q[i] = Vectors.sub(q[i], Vectors.scale(dot, q[j]))
	 end
	 local norm = Vectors.norm(q[i])
	 if math.abs(norm) < tolerance then
	    q[i] = nil
	    break
	 end
	 q[i] = Vectors.scale(1 / norm, q[i])
      end
      Q = SparseMatrix:new(q)
      -- print(Q)
      -- print(self)
      H = Q * self * Q:transpose()
      oldEigenvalues = Tools.list.copy(eigenvalues)
      eigenvalues = Tools.list.sublist(H:numericalClean(tolerance):toDense():eigenvalues(tolerance), n)
      if #oldEigenvalues < #eigenvalues then
	 oldEigenvalues = Tools.list.copy(eigenvalues)
	 oldEigenvalues[1] *= 2
	 oldEigenvalues[1] += tolerance
      end
      local temp = q[#q]
      q = {}
      q[1] = temp
   end
   return eigenvalues
end

function SparseMatrix:symmetricCopy()
   local n = math.min(self.length, self.width)
   local data = {}
   for i = 1, n do
      for j = 1, i do
	 data[n * (i - 1) + j] = self:get(i, j)
	 data[n * (j - 1) + i] = self:get(i, j)
      end
   end
   return SparseMatrix:new(data, n, n)
end

function SparseMatrix:tridiagonalize()
   -- We will assume that the matrix is symmetric and take the lower
   -- triangle as authoritative
   local matrix = self:symmetricCopy()
   for i = 1, self.length - 2 do
      -- i is the current column
      for j = i + 2, self.length do
	 local a = matrix:get(i + 1, i)
	 local b = matrix:get(j, i)
	 if b == 0 then
	    continue
	 end
	 local r = Numerics.hypot(a, b)
	 local c = a / r
	 local s = -b / r
	 matrix:set(i + 1, i, r)
	 matrix:set(i, i + 1, r)
	 matrix:set(j, i, 0)
	 matrix:set(i, j, 0)
	 for k = i + 1, self.length do
	    a = matrix:get(i + 1, k)
	    b = matrix:get(j, k)
	    matrix:set(i + 1, k, c * a - s * b)
	    matrix:set(j, k, s * a + c * b)
	    a = matrix:get(k, i + 1)
	    b = matrix:get(k, j)
	    matrix:set(k, i + 1, c * a - s * b)
	    matrix:set(k, j, s * a + c * b)
	 end
      end
   end
   return matrix
end

function SparseMatrix:QRDecomposition(tol)
   tol = tol or 10 ^ -13
   local QT = SparseMatrix:identity(self.length)
   local R = self:copy()
   for i = 1, self.width do
      for j = self.length, i + 1, -1 do
	 local a = R:get(i, i)
	 local b = R:get(j, i)
	 if math.abs(b) < tol then
	    continue
	 end
	 local r = Numerics.hypot(a, b)
	 R:set(i, i, r)
	 R:set(j, i, 0)
	 local c = a / r
	 local s = -b / r
	 for k = i + 1, self.width do
	    local x = R:get(i, k)
	    local y = R:get(j, k)
	    R:set(i, k, c * x - s * y)
	    R:set(j, k, s * x + c * y)
	 end
	 for k = 1, self.width do
	    local x = QT:get(i, k)
	    local y = QT:get(j, k)
	    QT:set(i, k, c * x - s * y)
	    QT:set(j, k, s * x + c * y)
	 end
      end
   end
   return QT:transpose(), R
end	 

function SparseMatrix:QREigenvalues(tol)
   tol = tol or 10 ^ -13
   if self.length ~= self.width then
      error("Sparse system is not square")
   end
   local matrix = self:tridiagonalize()
   local sub = matrix:getSubdiagonal()
   local err = Vectors.norm(sub) / #sub
   local q = SparseMatrix:identity(self.length)
   while err >= tol do
      local Q, R = matrix:QRDecomposition(tol)
      matrix = R * Q
      q = Q * q
      sub = matrix:getSubdiagonal()
      err = Vectors.norm(sub) / #sub
   end
   return matrix:getDiagonal(), q
end	 

function SparseMatrix:QREigenvaluesAnalysis(tol)
   tol = tol or 10 ^ -13
   if self.length ~= self.width then
      error("Sparse system is not square")
   end
   local matrix = self:tridiagonalize()
   local sub = matrix:getSubdiagonal()
   local err = Vectors.norm(sub) / #sub
   local q = SparseMatrix:identity(self.length)
   local qVec, dVec = {}, {}
   while err >= tol do
      local Q, R = matrix:QRDecomposition(tol)
      matrix = R * Q
      q = q * Q
      qVec[#qVec + 1] = q
      dVec[#dVec + 1] = matrix:getDiagonal()
      sub = matrix:getSubdiagonal()
      err = Vectors.norm(sub) / #sub
   end
   return matrix:getDiagonal(), q, qVec, dVec
end

function SparseMatrix:upperTriangularSolve(b: Vector)
   -- This overwrites the b vector
   if self.length ~= self.width then
      error("Sparse system is not square")
   end
   if self.length ~= #b then
      error("Sparse system has incompatibly sized b vector")
   end
   for i = #b, 1, -1 do
      temp = self:get(i, i)
      if temp == 0 then
	 error("Sparse system is not solvable")
      end
      for j = #b, i + 1, -1 do
	 temp2 = self:get(i, j)
	 if temp2 ~= 0 then
	    b[i] -= temp2 * b[j]
	 end
      end
      b[i] /= temp
   end
   return b
end

function SparseMatrix:hessenbergLeastSquare(b0: Vector)
   -- This overwrites the b vector and self
   -- We will compute this using Givens rotations!
   local length = self.length
   local width = self.width
   for i = 1, length - 1 do
      local index = (i - 1) * width + i
      local indexPlus = i * width + i
      local a, b = self.data[index], self.data[indexPlus]
      if b == nil then
	 continue
      end
      a = a or 0
      local r = Numerics.hypot(a, b)
      self.data[index], self.data[indexPlus] = r, 0
      local c, s = a / r, -b / r
      local iPlus = i + 1
      b0[i], b0[iPlus] = c * b0[i] - s * b0[iPlus], s * b0[i] + c * b0[iPlus]
      for j = i + 1, width do
	 index = (i - 1) * width + j
	 indexPlus = i * width + j
	 local temp1 = self.data[index]
	 local temp2 = self.data[indexPlus]
	 local sum1 = 0
	 local sum2 = 0
	 if temp1 ~= nil then
	    if temp2 ~= nil then
	       sum1 = c * temp1 - s * temp2
	       sum2 = s * temp1 + c * temp2
	    else
	       sum1 = c * temp1
	       sum2 = s * temp1
	    end
	 else
	    if temp2 ~= nill then
	       sum1 = -s * temp2
	       sum2 = c * temp2
	    end
	 end
	 self.data[index], self.data[indexPlus] = sum1, sum2
      end
   end
   local newLength = math.min(length, width)
   return self:submatrix(1, newLength, 1, newLength):upperTriangularSolve(Tools.list.sublist(b0, newLength))
end

function SparseMatrix:GMRES(b, tolerance)
   tolerance = tolerance or 10 ^ -13
   local x0 = Vectors.zeros(self.width) --Initial Guess
   local r0 = Vectors.sub(b, self:apply(x0))
   local beta = Vectors.norm(r0)
   local q, y, e
   while beta >= tolerance do
      q, h = self:arnoldiProcess(nil, r0)
      e = Vectors.zeros(#q + 1)
      e[1] = beta
      local Q = SparseMatrix:new(q, #q, #q[1]):transpose()
      y = h:hessenbergLeastSquare(e)
      x0 = Vectors.add(x0, Q:apply(y))
      r0 = Vectors.sub(b, self:apply(x0))
      beta = Vectors.norm(r0)
   end
   return x0
end

function SparseMatrix:solve(b, tol)
   return self:GMRES(b, tol)
end

function SparseMatrix:inverse(tol)
   tol = tol or 10 ^ -13
   if self.length ~= self.width then
      error("Cannot invert a non-square matrix")
   end
   local iden = Matrix:identity(self.length)
   local data = {}
   for i = 1, self.length do
      local e = iden:row(i)
      local r = self:solve(e)
      for j = 1, #r do
	 if math.abs(r[j]) >= tol then
	    data[self.width * (j - 1) + i] = r[j]
	 end
      end
   end
   return SparseMatrix:new(data, self.length, self.length)
end
	 
Matrices.SparseMatrix = SparseMatrix

--[[
   +--------------------------------------------------+
   |                 Matrix Constants                 |
   +--------------------------------------------------+
   |This section contains some tools for computing    |
   |constants related to matrix computation.          |
   +--------------------------------------------------+
]]

function Matrices.constants.ComputeStrassenLimit()
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

function Matrices.constants.ComputeSparseStrassenLimit()
   local strassenLimit = Matrices.constants.SPARSESTRASSENLIMIT
   for i = 1, 100, 1 do
      local matrix = SparseMatrix:identity(math.pow(2, i))
      Matrices.constants.SPARSESTRASSENLIMIT = math.floor(math.pow(2, i))
      local tic = os.clock()
      for j = 1, 1000, 1 do
	 local store = matrix * matrix
      end
      local time = (os.clock() - tic)
      Matrices.constants.SPARSESTRASSENLIMIT = math.floor(math.pow(2, i + 1))
      tic = os.clock()
      for j = 1, 1000, 1 do
	 local store = matrix * matrix
      end
      local time2 = (os.clock() - tic)
      if time < time2 then
	 print(time, time2)
	 Matrices.constants.SPARSESTRASSENLIMIT = math.floor(math.pow(2, i))
	 return nil
      end
   end
end

function Matrices.constants.ComputeComplexStrassenLimit()
   local strassenLimit = Matrices.constants.COMPLEXSTRASSENLIMIT
   for i = 1, 100, 1 do
      local matrix = ComplexMatrix:identity(math.pow(2, i))
      Matrices.constants.COMPLEXSTRASSENLIMIT = math.floor(math.pow(2, i))
      local tic = os.clock()
      for j = 1, 100, 1 do
	 local store = matrix * matrix
      end
      local time = (os.clock() - tic)
      Matrices.constants.COMPLEXSTRASSENLIMIT = math.floor(math.pow(2, i + 1))
      tic = os.clock()
      for j = 1, 100, 1 do
	 local store = matrix * matrix
      end
      local time2 = (os.clock() - tic)
      if time < time2 then
	 print(time, time2)
	 Matrices.constants.COMPLEXSTRASSENLIMIT = math.floor(math.pow(2, i))
	 return nil
      end
   end
end

return Matrices
