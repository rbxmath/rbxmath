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
            (self[i] or {})[j] or 0, (self[i] or {})[jPlus] or 0, (self[iPlus] or {})[j] or 0, (self[iPlus] or {})[jPlus] or 0
      end
   end
   return {
      Matrix:new(data1),
      Matrix:new(data2),
      Matrix:new(data3),
      Matrix:new(data4),
   }
end

function Matrix:column(i: number, start: number): Vector
   if i > self.width then
      error("Matrix doesn't have " .. tostring(i) .. " columns.")
   end
   start = start or 1
   local column = {}
   local ii = 1
   for j = start, self.length do
      column[ii] = self[j][i]
      ii += 1
   end
   return column
end

function Matrix:getColumn(i: number, start: number): Vector
   return self:column(i, start)
end

function Matrix:getColumnVector(i: number, start: number): Vector
   return self:column(i, start)
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
   if math.min(left.length, left.width, right.length, right.width) >= Matrices.constants.STRASSENLIMIT then
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
   |             Linear Algebra Subroutines           |
   +--------------------------------------------------+
   |Some common matrix operations can be dramatically |
   |improved if they are implemented all at once.     |
   |These are implemented here for use by "power      |
   |users" and for internal efficiency.               |
   +--------------------------------------------------+
   |As a guide to the reader the following symbols are|
   |used in function names:                           |
   |                                                  |
   | s   --- Scalar value                             |
   | m   --- Matrix                                   |
   | p   --- Plus                                     |
   | v   --- Vector                                   |
   | r1u --- Rank One Update                          |
   | b   --- Block                                    |
   +--------------------------------------------------+
--]]

--[[
   +--------------------------------------------------+
   |              Matrix Vector Operations            |
   +--------------------------------------------------+
--]]

--[=[
   Computes the matrix vector product Ax.

   @param matrix --- The matrix A in the matrix vector product Ax
   @param vector --- The vector x in the matrix vector product Ax
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The matrix vector product Ax.
--]=]
function Matrix.mv(matrix, vector, unsafe)
   if not unsafe and matrix.width ~= #vector then
      error("Matrix and vector are of incompatible size!")
   end
   local output = {}
   for i = 1, matrix.length do
      local row = matrix[i]
      local sum = 0
      for k, v in row do
	 sum += v * vector[k]
      end
      output[i] = sum
   end
   return output
end

--[=[
   Computes the vector matrix product y^T A.

   @param vector --- The vector y in the vector matrix product y^T A
   @param matrix --- The matrix A in the vector matrix product y^T A
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The vector matrix product y^T A.
--]=]
function Matrix.vm(vector, matrix, unsafe)
   if not unsafe and matrix.length ~= #vector then
      error("Matrix and vector are of incompatible size!")
   end
   local output = {}
   for k, v in ipairs(vector) do
      local sum = 0
      for i = 1, matrix.length do
	 sum += v * matrix[k][i]
      end
      output[k] = sum
   end
   return output
end

--[=[
   Computes the iterative matrix product ABx.

   @param matrix1 --- The matrix A in the matrix product ABx
   @param matrix2 --- The matrix B in the matrix product ABx
   @param vector  --- The vector x in the matrix product ABx
   @param unsafe  --- If true, no check will be done on the dimensions 

   @return --- The iterative matrix product ABx
--]=]
function Matrix.mmv(matrix1, matrix2, vector, unsafe)
   return Matrix.mv(matrix1, Matrix.mv(matrix2, vector, unsafe), unsafe)
end

--[=[
   Computes the iterative matrix product A^n x.

   @param matrix --- The matrix A in the matrix product A^n x
   @param power  --- The integer n in the matrix product A^n x
   @param vector --- The vector x in the matrix product A^n x
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The iterative matrix product A^n x
--]=]
function Matrix.imv(matrix, power, vector, unsafe)
   local output = Matrix.mv(matrix, vector, unsafe)
   for i = 2, power do
      output = Matrix.mv(matrix, output, unsafe)
   end
   return output
end

--[=[
   Computes the scalar-matrix-vector product cAx.

   @param scalar --- The scalar c in the scalar-matrix-vector product cAx
   @param matrix --- The matrix A in the scalar-matrix-vector product cAx
   @param vector --- The vector x in the scalar-matrix-vectro product cAx
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The scalar-matrix-vector product cAx
--]=]
function Matrix.smv(scalar, matrix, vector, unsafe)
   if not unsafe and matrix.width ~= #vector then
      error("Matrix and vector are of incompatible size!")
   end
   local output = {}
   for i = 1, matrix.length do
      local row = matrix[i]
      local sum = 0
      for k, v in row do
	 sum += v * vector[k]
      end
      output[i] = scalar * sum
   end
   return output
end

--[=[
   Computes the vector sum Ax + y.

   @param matrix  --- The matrix A in the sum Ax + y
   @param vector1 --- The vector x in the sum Ax + y
   @param vector2 --- The vector y in the sum Ax + y
   @param unsage  --- If true, no checks will be done on the dimensions

   @return --- The sum Ax + y
--]=]
function Matrix.mvpv(matrix, vector1, vector2, unsafe)
   if not unsafe and matrix.width ~= #vector1 then
      error("Matrix and vector are of incompatible size!")
   elseif not unsafe and matrix.length ~= #vector2 then
      error("Attempting to add vectors of incompatible size!")
   end
   local output = {}
   for i = 1, matrix.length do
      local row = matrix[i]
      local sum = vector2[i]
      for k, v in row do
	 sum += v * vector1[k]
      end
      output[i] = sum
   end
   return output
end

--[=[
   Computes the vector sum cAx + y.

   @param scalar  --- The scalar c in the sum cAx + y
   @param matrix  --- The matrix A in the sum cAx + y
   @param vector1 --- The vector x in the sum cAx + y
   @param vector2 --- The vector y in the sum cAx + y
   @param unsage  --- If true, no checks will be done on the dimensions

   @return --- The sum cAx + y
--]=]
function Matrix.smvpv(scalar, matrix, vector1, vector2, unsafe)
   if not unsafe and matrix.width ~= #vector1 then
      error("Matrix and vector are of incompatible size!")
   elseif not unasfe and matrix.length ~= #vector2 then
      error("Attempting to add vectors of incompatible size!")
   end
   local output = {}
   for i = 1, matrix.length do
      local row = matrix[i]
      local sum = 0
      for k, v in row do
	 sum += v * vector1[k]
      end
      output[i] = vector2[i] + scalar * sum
   end
   return output
end

--[=[
   Computes the vector sum cAx + dy.

   @param scalar1 --- The scalar c in the sum cAx + dy
   @param matrix  --- The matrix A in the sum cAx + dy
   @param vector1 --- The vector x in the sum cAx + dy
   @param scalar2 --- The scalar d in the sum cAx + dy
   @param vector2 --- The vector y in the sum cAx + dy
   @param unsage  --- If true, no checks will be done on the dimensions

   @return --- The sum cAx + dy
--]=]
function Matrix.smvpsv(scalar1, matrix, vector1, scalar2, vector2, unsafe)
   if not unsafe and matrix.width ~= #vector1 then
      error("Matrix and vector are of incompatible size!")
   elseif not unsafe and matrix.length ~= #vector2 then
      error("Attempting to add vectors of incompatible size!")
   end
   local output = {}
   for i = 1, matrix.length do
      local row = matrix[i]
      local sum = 0
      for k, v in row do
	 sum += v * vector1[k]
      end
      output[i] = scalar2 * vector2[i] + scalar1 * sum
   end
   return output
end

--[=[
   Computes the vector sum cAx + dBy.

   @param scalar1 --- The scalar c in the sum cAx + dBy
   @param matrix1 --- The matrix A in the sum cAx + dBy
   @param vector1 --- The vector x in the sum cAx + dBy
   @param scalar2 --- The scalar d in the sum cAx + dBy
   @param matrix2 --- The matrix B in the sum cAx + dBy
   @param vector2 --- The vector y in the sum cAx + dBy
   @param unsage  --- If true, no checks will be done on the dimensions

   @return --- The sum cAx + dBy
--]=]
function Matrix:smvpsmv(scalar1, matrix1, vector1, scalar2, matrix2, vector2, unsafe)
   if not unsafe and matrix1.width ~= #vector1 then
      error("Matrix and vector are of incompatible size!")
   elseif not unsafe and matrix2.width ~= #vector2 then
      error("Matrix and vector are of incompatible size!")
   elseif not unsafe and matrix1.length ~= matrix2.length then
      error("Attempting to add vectors of incompatible size!")
   end
   local output = {}
   for i = 1, #vector1 do
      local sum1 = 0
      local sum2 = 0
      local row1 = matrix1[i]
      local row2 = martix2[i]
      for k, v in row1 do
	 sum1 += v * vector1[k]
      end
      for k, v in row2 do
	 sum2 += v * vector2[k]
      end
      output[i] = scalar1 * sum1 + scalar2 * sum2
   end
   return output
end

--[=[
   Performs the rank one update A + cxy^T. Overwrites the matrix A.

   @param self   --- The matrix A in the rank one update A + cxy^T
   @param scalar --- The scalar c in the rank one update A + cxy^T
   @param column --- The vector x in the rank one update A + cxy^T
   @param row    --- The vector y in the rank one update A + cxy^T
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T
--]=]
function Matrix.r1u(self, scalar, column, row, unsafe)
   if not unsafe and (matrix.length ~= #column or matrix.width ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = 1, self.length do
      local row2 = self[i]
      for j = 1, self.width do
	 row2[j] += col[i] * row[j]
      end
   end
   return self
end

--[=[
   Performs the rank one update dA + cxy^T. Overwrites the matrix A.

   @param self    --- The matrix A in the rank one update dA + cxy^T
   @param scalar1 --- The scalar d in the rank one update dA + cxy^T
   @param scalar2 --- The scalar c in the rank one update dA + cxy^T
   @param column  --- The vector x in the rank one update dA + cxy^T
   @param row     --- The vector y in the rank one update dA + cxy^T
   @param unsafe  --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.sr1u(self, scalar1, scalar2, column, row, unsafe)
   if not unsafe and (matrix.length ~= #column or matrix.width ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = 1, self.length do
      local row2 = self[i]
      for j = 1, self.width do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return self
end

--[=[
   Performs the rank one update A + cxy^T.

   @param matrix --- The matrix A in the rank one update A + cxy^T
   @param scalar --- The scalar c in the rank one update A + cxy^T
   @param column --- The vector x in the rank one update A + cxy^T
   @param row    --- The vector y in the rank one update A + cxy^T
   @param unsafe --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T
--]=]
function Matrix.mpr1u(matrix, scalar, column, row, unsafe)
   if not unsafe and (matrix.length ~= #column or matrix.width ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix:copy()
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = 1, output.length do
      local row2 = output[i]
      for j = 1, output.width do
	 row2[j] += col[i] * row[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update dA + cxy^T.

   @param scalar1 --- The scalar d in the rank one update dA + cxy^T
   @param matrix  --- The matrix A in the rank one update dA + cxy^T
   @param scalar2 --- The scalar c in the rank one update dA + cxy^T
   @param column  --- The vector x in the rank one update dA + cxy^T
   @param row     --- The vector y in the rank one update dA + cxy^T
   @param unsafe  --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.smpr1u(scalar1, matrix, scalar2, column, row, unsafe)
   if not unsafe and (matrix.length ~= #column or matrix.width ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix:copy()
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = 1, output.length do
      local row2 = output[i]
      for j = 1, output.width do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update B + cxy^T on a submatrix of the matrix A.

   @param self     --- The matrix A above
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar   --- The scalar c in the rank one update B + cxy^T
   @param column   --- The vector x in the rank one update B + cxy^T
   @param row      --- The vector y in the rank one update B + cxy^T
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.br1u(self, rowStart, rowEnd, colStart, colEnd, scalar, column, row, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or colEnd - colStart + 1 ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   local rowIndex = 1
   local colIndex = 1
   for i = rowStart, rowEnd do
      local matRow = self[i]
      local rowConstant = row[rowIndex]
      for j = colStart, colEnd do
	 matRow[j] += col[colIndex] * rowConstant
	 colIndex += 1
      end
      colIndex = 1
      rowConstant += 1
   end
   return self
end

--[=[
   Performs the rank one update B + cxy^T on a submatrix of the matrix A.

   @param matrix   --- The matrix A above
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar   --- The scalar c in the rank one update B + cxy^T
   @param column   --- The vector x in the rank one update B + cxy^T
   @param row      --- The vector y in the rank one update B + cxy^T
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.bpr1u(matrix, rowStart, rowEnd, colStart, colEnd, scalar, column, row, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or colEnd - colStart + 1 ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = self:copy()
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   local rowIndex = 1
   local colIndex = 1
   for i = rowStart, rowEnd do
      local matRow = output[i]
      local rowConstant = row[rowIndex]
      for j = colStart, colEnd do
	 matRow[j] += col[colIndex] * rowConstant
	 colIndex += 1
      end
      colIndex = 1
      rowConstant += 1
   end
   return output
end

--[=[
   Performs the rank one update dB + cxy^T on a submatrix of the matrix A.

   @param self     --- The matrix A above
   @param scalar1  --- The scalar d in the rank one update dB + cxy^T
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar2  --- The scalar c in the rank one update dB + cxy^T
   @param column   --- The vector x in the rank one update dB + cxy^T
   @param row      --- The vector y in the rank one update dB + cxy^T
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.sbr1u(self, scalar1, rowStart, rowEnd, colStart, colEnd, scalar2, column, row, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or colEnd - colStart + 1 ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   local rowIndex = 1
   local colIndex = 1
   for i = rowStart, rowEnd do
      local matRow = self[i]
      local rowConstant = row[rowIndex]
      for j = colStart, colEnd do
	 matRow[j] = scalar1 * matRow[j] + col[colIndex] * rowConstant
	 colIndex += 1
      end
      colIndex = 1
      rowConstant += 1
   end
   return self
end

--[=[
   Performs the rank one update dB + cxy^T on a submatrix of the matrix A.

   @param scalar1  --- The scalar d in the rank one update dB + cxy^T
   @param matrix   --- The matrix A above
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar2  --- The scalar c in the rank one update dB + cxy^T
   @param column   --- The vector x in the rank one update dB + cxy^T
   @param row      --- The vector y in the rank one update dB + cxy^T
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T
--]=]
function Matrix.sbpr1u(scalar1, matrix, rowStart, rowEnd, colStart, colEnd, scalar2, column, row, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or colEnd - colStart + 1 ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = self:copy()
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   local rowIndex = 1
   local colIndex = 1
   for i = rowStart, rowEnd do
      local matRow = output[i]
      local rowConstant = row[rowIndex]
      for j = colStart, colEnd do
	 matRow[j] = scalar1 * matRow[j] + col[colIndex] * rowConstant
	 colIndex += 1
      end
      colIndex = 1
      rowConstant += 1
   end
   return output
end

--[=[
   Performs the rank one update A + cxy^T B.

   @param matrix1  --- The matrix A in the rank one update A + cxy^T B
   @param scalar   --- The scalar c in the rank one update A + cxy^T B
   @param column   --- The vector x in the rank one update A + cxy^T B
   @param row      --- The vector y in the rank one update A + cxy^T B
   @param matrix2  --- The matrix B in the rank one update A + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T B.
--]=]
function Matrix.mpr1um(matrix1, scalar, column, row, matrix2, unsafe)
   if not unsafe and (matrix1.length ~= #column or
		      matrix1.width ~= matrix2.width or
		      matrix2.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix1:copy()
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix2.length do
	 sum += v * matrix2[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = 1, output.length do
      local row2 = output[i]
      for j = 1, output.width do
	 row2[j] += col[i] * rol[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update dA + cxy^T B.

   @param scalar1  --- The scalar d in the rank one update dA + cxy^T B
   @param matrix1  --- The matrix A in the rank one update dA + cxy^T B
   @param scalar2  --- The scalar c in the rank one update dA + cxy^T B
   @param column   --- The vector x in the rank one update dA + cxy^T B
   @param row      --- The vector y in the rank one update dA + cxy^T B
   @param matrix2  --- The matrix B in the rank one update dA + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T B.
--]=]
function Matrix.smpr1um(scalar1, matrix1, scalar2, column, row, matrix2, unsafe)
   if not unsafe and (matrix1.length ~= #column or
		      matrix1.width ~= matrix2.width or
		      matrix2.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix1:copy()
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix2.length do
	 sum += v * matrix2[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = 1, output.length do
      local row2 = output[i]
      for j = 1, output.width do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update A + cxy^T B.

   @param self     --- The matrix A in the rank one update A + cxy^T B
   @param scalar   --- The scalar c in the rank one update A + cxy^T B
   @param column   --- The vector x in the rank one update A + cxy^T B
   @param row      --- The vector y in the rank one update A + cxy^T B
   @param matrix   --- The matrix B in the rank one update A + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T B.
--]=]
function Matrix.r1um(self, scalar, column, row, matrix, unsafe)
   if not unsafe and (self.length ~= #column or
		      self.width ~= matrix.width or
		      matrix.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix.length do
	 sum += v * matrix[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = 1, self.length do
      local row2 = self[i]
      for j = 1, self.width do
	 row2[j] += col[i] * rol[j]
      end
   end
   return self
end

--[=[
   Performs the rank one update dA + cxy^T B.

   @param self     --- The matrix A in the rank one update dA + cxy^T B
   @param scalar1  --- The scalar d in the rank one update dA + cxy^T B
   @param scalar2  --- The scalar c in the rank one update dA + cxy^T B
   @param column   --- The vector x in the rank one update dA + cxy^T B
   @param row      --- The vector y in the rank one update dA + cxy^T B
   @param matrix   --- The matrix B in the rank one update dA + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T B.
--]=]
function Matrix.sr1um(self, scalar1, scalar2, column, row, matrix, unsafe)
   if not unsafe and (self.length ~= #column or
		      self.width ~= matrix.width or
		      matrix.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix.length do
	 sum += v * matrix[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = 1, self.length do
      local row2 = self[i]
      for j = 1, self.width do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return self
end

--[=[
   Performs the block rank one update A + cxy^T B.

   @param matrix1  --- The matrix A in the rank one update A + cxy^T B
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar   --- The scalar c in the rank one update A + cxy^T B
   @param column   --- The vector x in the rank one update A + cxy^T B
   @param row      --- The vector y in the rank one update A + cxy^T B
   @param matrix2  --- The matrix B in the rank one update A + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T B.
--]=]
function Matrix.bpr1um(matrix1, rowStart, rowEnd, colStart, colEnd, scalar, column, row, matrix2, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or
		      colEnd - colStart + 1 ~= matrix2.width or
		      matrix2.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix1:copy()
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix2.length do
	 sum += v * matrix2[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = rowStart, rowEnd do
      local row2 = output[i]
      for j = colStart, colEnd do
	 row2[j] += col[i] * rol[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update dA + cxy^T B.

   @param scalar1  --- The scalar d in the rank one update dA + cxy^T B
   @param matrix1  --- The matrix A in the rank one update dA + cxy^T B
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar2  --- The scalar c in the rank one update dA + cxy^T B
   @param column   --- The vector x in the rank one update dA + cxy^T B
   @param row      --- The vector y in the rank one update dA + cxy^T B
   @param matrix2  --- The matrix B in the rank one update dA + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T B.
--]=]
function Matrix.sbpr1um(scalar1, matrix1, rowStart, rowEnd, colStart, colEnd, scalar2, column, row, matrix2, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or
		      colEnd - colStart + 1 ~= matrix2.width or
		      matrix2.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local output = matrix1:copy()
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix2.length do
	 sum += v * matrix2[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = rowStart, rowEnd do
      local row2 = output[i]
      for j = colStart, colEnd do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return output
end

--[=[
   Performs the rank one update A + cxy^T B.

   @param self     --- The matrix A in the rank one update A + cxy^T B
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar   --- The scalar c in the rank one update A + cxy^T B
   @param column   --- The vector x in the rank one update A + cxy^T B
   @param row      --- The vector y in the rank one update A + cxy^T B
   @param matrix   --- The matrix B in the rank one update A + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update A + cxy^T B.
--]=]
function Matrix.br1um(self, rowStart, rowEnd, colStart, colEnd, scalar, column, row, matrix, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or
		      colEnd - colStart + 1 ~= matrix.width or
		      matrix.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix.length do
	 sum += v * matrix[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar * column[i]
   end
   for i = rowStart, rowEnd do
      local row2 = self[i]
      for j = colStart, colEnd do
	 row2[j] += col[i] * rol[j]
      end
   end
   return self
end

--[=[
   Performs the rank one update dA + cxy^T B.

   @param self     --- The matrix A in the rank one update dA + cxy^T B
   @param rowStart --- The start row for the submatrix
   @param rowEnd   --- The end row for the submatrix
   @param colStart --- The start column for the submatrix
   @param colEnd   --- The end column for the submatrix
   @param scalar1  --- The scalar d in the rank one update dA + cxy^T B
   @param scalar2  --- The scalar c in the rank one update dA + cxy^T B
   @param column   --- The vector x in the rank one update dA + cxy^T B
   @param row      --- The vector y in the rank one update dA + cxy^T B
   @param matrix   --- The matrix B in the rank one update dA + cxy^T B
   @param unsafe   --- If true, no checks will be done on the dimensions

   @return --- The rank one update dA + cxy^T B.
--]=]
function Matrix.sbr1um(self, rowStart, rowEnd, colStart, colEnd, scalar1, scalar2, column, row, matrix, unsafe)
   if not unsafe and (rowEnd - rowStart + 1 ~= #column or
		      colEnd - colStart + 1 ~= matrix.width or
		      matrix.length ~= #row) then
      error("Incompatible dimensions for rank one update!")
   end
   local rol = {}
   for k, v in ipairs(row) do
      local sum = 0
      for i = 1, matrix.length do
	 sum += v * matrix[k][i]
      end
      rol[k] = sum
   end
   local col = {}
   -- Cheaper to do O(n) multiplies than n^2 multiplies
   for i = 1, #column do
      col[i] = scalar2 * column[i]
   end
   for i = rowStart, rowEnd do
      local row2 = self[i]
      for j = colStart, colEnd do
	 row2[j] = scalar1 * row2[j] + col[i] * row[j]
      end
   end
   return self
end

function Matrix.lr(self, r1, r2, c, s)
   for i = 1, self.width do
      local a, b = self[r1][i], self[r2][i]
      if a == 0 and b == 0 then
	 continue
      end
      self[r1][i], self[r2][i] = c * a - s * b, s * a + c * b
   end
   return self
end

function Matrix.tsol(self, vec, tol, unsafe)
   if not unsafe and self.length ~= #vec then
      error("Incompatible dimensions in triangular system vec!")
   end
   tol = tol or 10^-13
   for i = #self, 1, -1 do
      for j = i + 1, #self do
	 vec[i] -= vec[j] * self[i][j]
      end
      local diag = self[i][i]
      if math.abs(diag) <= tol then
	 if math.abs(vec[i]) <= tol then
	    continue
	 else
	    error("Triangular system is not solvable!")
	 end
      end
      vec[i] /= diag
   end
   return vec
end

function Matrix.tinv(self, tol, unsafe)
   if not unsafe and self.length ~= self.width then
      error("Cannot invert nonsquare triangular matrix!")
   end
   tol = tol or 10^-13
   for i = 1, self.length do
      if math.abs(self[i][i]) <= tol then
	 error("Upper triangular matrix is not invertible!")
      end
      self[i][i] = 1 / self[i][i]
   end
   for i = 1, self.length - 1 do
      for j = i + 1, self.width do
	 self[i][j] *= -self[i][i] * self[j][j]
      end
   end
   return self
end

function Matrix.stpstTinv(self, scalar1, scalar2, tol, unsafe)
   if not unsafe and self.length ~= self.width then
      error("Cannot invert nonsquare triangular matrix!")
   end
   local n = self.length
   tol = tol or 10^-13
   local diag = table.create(n, 0)
   for i = 1, n do
      diag[i] = self[i][i]
   end
   for i = 1, self.length do
      if math.abs(self[i][i]) <= tol then
	 error("Upper triangular matrix is not invertible!")
      end
      self[i][i] = 1 / self[i][i]
   end
   for i = 1, self.length - 1 do
      for j = i + 1, self.width do
	 self[j][i] = -scalar2 * self[i][i] * self[j][j] * self[i][j]
      end
   end
   for i = 1, self.length - 1 do
      for j = i + 1, self.width do
	 self[i][j] *= scalar1
      end
   end
   for i = 1, n do
      self[i][i] = scalar1 * diag[i] + scalar2 * self[i][i]
   end
   return self
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
   |                  Boolean Functions               |
   +--------------------------------------------------+
   |Functions that return boolean values about        |
   |matrices
   +--------------------------------------------------+
]]

function Matrix:qrIsOrthogonal(tol)
   for i = 1, self.length do
      for j = i, self.width do
	 if i == j then
	    if math.abs(math.abs(self[i][j]) - 1) >= tol then
	       return false
	    end
	 else
	    if math.abs(self[i][j]) >= tol then
	       return false
	    end
	 end
      end
   end
   return true
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

--[=[
   Computes an LU factorization of the matrix. Partial pivoting will be used to guarantee the factorization works.
   The matrix will be returned in the condensed LU form, namely the entire LU will be stored in a single matrix.
   This will overwrite the initial matrix.
--]=]
function Matrix:LU(tol)
   if self.length ~= self.width then
      error("Cannot compute LU factorization of rectangular matrix!")
   end

   tol = tol or 10^-13

   local n = self.length

   local permutation = {}
   for i = 1, n do
      permutation[i] = i
   end

   for i = 1, n - 1 do
      -- We now must find the largest element in the ith column
      local max = 0
      local ind = 0
      for j = i, n do
	 local val = math.abs(self[j][i])
	 if val > max then
	    max = val
	    ind = j
	 end
      end

      if math.abs(max) <= tol then
	 for j = i, n do
	    self[j][i] = 0
	 end
	 continue
      end

      self[i], self[ind] = self[ind], self[i]
      permutation[i], permutation[ind] = permutation[ind], permutation[i]

      local pivot = self[i][i]
      local pivml = 1 / pivot

      -- We will now eliminate everything below the pivot!
      for j = i + 1, n do
	 local target = self[j][i]
	 if target ~= 0 then
	    local l = target * pivml
	    self[j][i] = l
	    for k = i + 1, n do
	       self[j][k] -= l * self[i][k]
	    end
	 end
      end
   end
   return self, permutation
end

--[=[
   Solves the linear system Ax = b assuming that the matrix is already in the condensed LU form of A.

   @param b    --- The vector b in the linear system Ax = b
   @param perm --- The permutation vector returned in the computation of the LU
   @param tol  --- The tolerance value for determining if something is zero
--]=]
function Matrix:LUSolve(b, perm, tol)
   tol = tol or 10^-13
   local bPerm = {}
   for k, v in ipairs(perm) do
      bPerm[v] = b[k]
   end
   local bInvs = {bPerm[1]}
   for i = 2, #bPerm do
      local sum = bPerm[i]
      for j = 1, i - 1 do
	 sum -= self[i][j] * bPerm[j]
      end
      bInvs[i] = sum
   end
   for i = #self, 1, -1 do
      bPerm[i] = bInvs[i]
      for j = i + 1, #self do
	 bPerm[i] -= bPerm[j] * self[i][j]
      end
      local diag = self[i][i]
      if math.abs(diag) <= tol then
	 if math.abs(bPerm[i]) <= tol then
	    continue
	 else
	    error("System is not solvable!")
	 end
      end
      bPerm[i] /= diag
   end
   return bPerm
end

--[=[
   Computes an LDU factorization of the matrix. Partial pivoting will be used to guarantee the factorization works.
   The matrix will be returned in the condensed LDU form, namely the entire LDU will be stored in a single matrix.
   This will overwrite the initial matrix.
--]=]
function Matrix:LDU()
   if self.length ~= self.width then
      error("Cannot compute LU factorization of rectangular matrix!")
   end

   local n = self.length

   local permutation = {}
   for i = 1, n do
      permutation[i] = i
   end

   for i = 1, n - 1 do
      -- We now must find the largest element in the ith column
      local max = 0
      local ind = 0
      for j = i, n do
	 local val = math.abs(self[j][i])
	 if val > max then
	    max = val
	    ind = j
	 end
      end

      if max == 0 then
	 error("Matrix is not invertible!")
      end

      self[i], self[ind] = self[ind], self[i]
      permutation[i], permutation[ind] = permutation[ind], permutation[i]

      local pivot = self[i][i]
      local pivml = 1 / pivot

      -- We will now eliminate everything below the pivot!
      for j = i + 1, n do
	 local target = self[j][i]
	 if target ~= 0 then
	    local l = target * pivml
	    self[j][i] = l
	    for k = i + 1, n do
	       self[j][k] -= l * self[i][k]
	    end
	 end
      end
   end

   for i = 1, n - 1 do
      for j = i + 1, n do
	 self[i][j] /= self[i][i]
      end
   end
   
   return self, permutation
end

function Matrix:CholeskyDecomposition(tol)
   tol = tol or 10^-13
   local rtt = math.sqrt(tol)
   local matrix = self:copy()
   local numberOfRows = #matrix
   local numberOfColumns = #matrix[1]

   if numberOfRows ~= numberOfColumns then
      error("Cannot compute Cholesky of rectangular matrix.")
   end

   for i = 1, numberOfColumns do
      for j = 1, i - 1 do
	 local entry = matrix[j][i]
	 if math.abs(entry) >= rtt then
	    matrix[i][i] -= entry ^ 2
	 end
      end
      if matrix[i][i] < -tol then
         error("Matrix is not positive definite.")
      end
      matrix[i][i] = math.sqrt(math.max(matrix[i][i], tol))
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
      local ithCol = a:column(i, i)

      a:br1um(i, a.length, i + 1, a.width, -gamma, ithCol, ithCol, everythingButIthCol)

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
      qSub = Q:submatrix(i, #R, i, #R)
      Q:br1um(i, #R, i, #R, -gammaVector[i], u, u, qSub)
   end
   return Q, R
end

function Matrix:FullHouseholderQR(tol)
   local a, g = self:HouseholderQR(tol)
   return Matrix.ExpandHouseholderQR(a, g)
end

function Matrix:getGivens()
   -- Output will be of the form {length, width, theta1, theta2, ...}
   -- column indexed and then row
   local givensList = {self.length}
   for i = 2, self.length do
      for j = 1, i - 1 do
	 givensList[#givensList + 1] = self[i][j]
      end
   end
   return givensList
end

function Matrix.applyGivensTransposeToVector(givens, vector)
   local length = givens[1]
   local iter = 2
   -- i is the current column
   for i = 1, length - 1 do
      -- j is the current row
      for j = i + 1, length do
	 local theta = givens[iter]
	 iter += 1
	 local c = math.cos(theta)
	 local s = math.sin(theta)
	 vector[i], vector[j] = c * vector[i] - s * vector[j],
	    s * vector[i] + c * vector[j]
      end
   end
   return vector
end

function Matrix.applyGivensTransposeToMatrix(givens, matrix)
   local length = givens[1]
   local iter = 2
   -- i is the current column
   for i = 1, length - 1 do
      -- j is the current row
      for j = i + 1, length do
	 local theta = givens[iter]
	 iter += 1
	 local c = math.cos(theta)
	 local s = math.sin(theta)
	 matrix:lr(i, j, c, s)
      end
   end
   return matrix
end

function Matrix.applyGivensToMatrix(givens, matrix)
   local length = givens[1]
   local iter = #givens
   -- i is the current column
   for i = length - 1, 1 do
      -- j is the current row
      for j = length, i + 1 do
	 local theta = givens[iter]
	 iter -= 1
	 local c = math.cos(theta)
	 local s = -math.sin(theta)
	 matrix:lr(i, j, c, s)
      end
   end
   return matrix
end

function Matrix:GivensQR(tol)
   -- Overwrites self
   tol = tol or 10 ^ -13
   local R = self
   for i = 1, self.width do
      for j = self.length, i + 1, -1 do
         local a = R[i][i]
         local b = R[j][i]
         if math.abs(b) < tol then
            continue
         end
         local r = Numerics.hypot(a, b)
         local c = a / r
         local s = -b / r
         R[i][i], R[j][i] = r, math.atan2(s, c)
         for k = i + 1, self.width do
            local x = R[i][k]
            local y = R[j][k]
            R[i][k] = c * x - s * y
            R[j][k] = s * x + c * y
         end
      end
   end
   return R
end

function Matrix:ExpandGivensQR()
   local n = math.min(self.length, self.width)
   local Q = Matrix.identity(n)
   for col = n - 1, 1, -1 do
      for row = n, col - 1, -1 do
	 local c = math.cos(self[row][col])
	 self[row][col] = 0
	 local s = -math.sin(self[row][col])
	 Q:lr(col, row, c, s)
      end
   end
   return Q, R
end

function Matrix:FullGivensQR(tol)
   -- Overwrites self
   local R = self:GivensQR(tol)
   return R:ExpandGivensQR()
end

function Matrix:GivensQRSolve(vec, tol)
   -- Doesn't overwrite self
   local n = math.min(self.length, self.width)
   for col = n - 1, 1, -1 do
      for row = n, col + 1, -1 do
	 local c = math.cos(self[row][col])
	 local s = math.sin(self[row][col])
	 vec[col], vec[row] = c * vec[col] - s * vec[row],
	    s * vec[col] + c * vec[row]
      end
   end
   return self:tsol(vec)
end

-- Not working yet, need to to right multiplication by Q
function Matrix:GivensQRInvert(tol)
   local n = self.length
   self:tinv(tol)
   for col = n - 1, 1, -1 do
      for row = n, col + 1, -1 do
	 local c = math.cos(self[row][col])
	 local s = math.sin(self[row][col])
	 self[row][col] = 0
	 for i = col, n do
	    self[col][i], self[row][i] = c * self[col][i] - s * self[row][i], s * self[col][i] + c * self[row][i]
	 end
      end
   end
   return self
end

function Matrix.GivensQ(givens)
   local length = givens[1]
   local matrix = Matrix:identity(length)
   return Matrix.applyGivensToMatrix(givens, matrix)
end

function Matrix:polarTerm(tol)
   -- Doesn't overwrite self
   if self.length ~= self.width then
      error("Nonsquare matrices currently unsupported for polar decomposition!")
   end
   tol = tol or 10^-13
   local matrix = self:copy()
   matrix:GivensQR(tol)
   local givens = Matrix.getGivens(matrix)
   while not matrix:qrIsOrthogonal(tol) do
      local max, min = math.abs(matrix[1][1]), math.abs(matrix[1][1])
      for i = 2, self.length do
	 local val = math.abs(matrix[i][i])
	 if val > max then
	    max = val
	 elseif val < min then
	    min = val
	 end
      end
      local gamma = math.sqrt(1 / (max * min))
      matrix:stpstTinv(0.5 * gamma, 0.5 / gamma, tol, true)
      Matrix.applyGivensToMatrix(givens, matrix)
      matrix:GivensQR(tol)
      givens = Matrix.getGivens(matrix)
   end
   matrix = Matrix:identity(self.length)
   return Matrix.applyGivensToMatrix(givens, matrix)
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

function Matrix:entrywiseMax()
   local max = 0
   for i = 1, self.length do
      for j = 1, self.width do
	 max = math.max(max, math.abs(self[i][j]))
      end
   end
   return max
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

function Matrix:diagonal(list)
   local matrix = Matrix:identity(#list)
   for i = 1, #list do
      matrix[i][i] = list[i]
   end
   return matrix
end

function Matrix:getDiagonal()
   local data = {}
   for i = 1, math.min(self.length, self.width) do
      data[i] = self[i][i]
   end
   return data
end

function Matrix:bidiagonalize(u, v, tol)
   tol = tol or 10 ^ -13
   local matrix = self:copy()
   u = u or Matrix:identity(self.length)
   v = v or Matrix:identity(self.width)
   -- We will proceed in two reductions
   for i = 1, math.max(self.length - 1, self.width - 1) do
      -- We will roll up the ith column, if there is one!
      if i < self.length then
	 for j = i + 1, self.length do
	    local a = matrix[i][i]
	    local b = matrix[j][i]
	    if math.abs(b) <= tol then
	       continue
	    end
	    local r = Numerics.hypot(a, b)
	    if r <= tol then
	       matrix[i][i], matrix[j][i] = 0, 0
	    end
	    -- I want my diagonal to be positive because
	    -- I want positive singular values
	    if a < 0 then
	       r = -r
	    end
	    local c = a / r
	    local s = -b / r
	    matrix[i][i], matrix[j][i] = r, 0
	    for k = 1, self.length do
	       u[i][k], u[j][k] = c * u[i][k] - s * u[j][k], s * u[i][k] + c * u[j][k]
	    end
	    for k = i + 1, self.width do
	       matrix[i][k], matrix[j][k] = c * matrix[i][k] - s * matrix[j][k],
		  s * matrix[i][k] + c * matrix[j][k]
	    end
	 end
      end
      -- We will roll up the ith row, if there is one!
      if i < self.width then
	 local col = i + 1
	 for j = col + 1, self.width do
	    local a = matrix[i][col]
	    local b = matrix[i][j]
	    if math.abs(b) <= tol then
	       continue
	    end
	    local r = Numerics.hypot(a, b)
	    if r <= tol then
	       matrix[i][col], matrix[i][j] = 0, 0
	    end
	    local c = a / r
	    local s = -b / r
	    matrix[i][col], matrix[i][j] = r, 0
	    for k = 1, self.width do
	       v[col][k], v[j][k]
		  = c * v[col][k] - s * v[j][k], s * v[col][k] + c * v[j][k]
	    end
	    for k = i + 1, self.length do
	       matrix[k][col], matrix[k][j] = c * matrix[k][col] - s * matrix[k][j],
		  s * matrix[k][col] + c * matrix[k][j]
	    end
	 end
      end
   end
   return matrix, u, v
end

function Matrix:gkSVDOne(u, v, tol, maxIters)
   u = u or Matrix:identity(self.length)
   v = v or Matrix:identity(self.width)
   tol = tol or 10 ^ -13
   maxIters = maxIters or 10000
   local n = math.min(self.length, self.width)
   local matrix = self:submatrix(1, n, 1, n)
   for i = 1, maxIters do
      for i = 1, n - 1 do
	 if math.abs(matrix[i][i + 1]) < tol then
	    local sig1, u1, v1
	       = matrix:submatrix(1, i, 1, i):gkSVDOne(nil,
						       nil,
						       tol,
						       maxIters - i)
	    local sig2, u2, v2
	       = matrix:submatrix(i + 1, n, i + 1, n):gkSVDOne(nil,
							       nil,
							       tol,
							       maxIters - i)
	    u:setSubmatrix(1, i, 1, n, u1 * u:submatrix(1, i, 1, n))
	    u:setSubmatrix(i + 1, n, 1, n, u2 * u:submatrix(i + 1, n, 1, n))
	    v:setSubmatrix(1, i, 1, n, v1 * v:submatrix(1, i, 1, n))
	    v:setSubmatrix(i + 1, n, 1, n, v2 * v:submatrix(i + 1, n, 1, n))
	    return Tools.list.join(sig1, sig2), u, v
	 end
      end
      local lambda1, lambda2, lambda
      if n == 1 then
	 return {matrix[1][1]}, Matrix:new({{1}}), Matrix:new({{1}})
      elseif n == 2 then
	 local a = matrix[1][1]^2
	 local b = matrix[1][2]^2
	 local c = matrix[2][2]^2
	 lambda1 = 0.5 * (a + b + c - math.sqrt(-4 * a * c + (a + b + c)^2))
	 lambda2 = 0.5 * (a + b + c + math.sqrt(-4 * a * c + (a + b + c)^2))
	 if math.abs(lambda1 - b - c) < math.abs(lambda2 - b - c) then
	    lambda = lambda1
	 else
	    lambda = lambda2
	 end
      else
	 local a = matrix[n - 2][n - 1]^2
	 local b = matrix[n - 1][n - 1]^2
	 local c = matrix[n - 1][n]^2
	 local d = matrix[n][n]^2
	 lambda1 = 0.5 * (a + b + c + d +
			  math.sqrt((a + b + c + d)^2 - 4 * (a * c + a * d + b * d)))
	 lambda2 = 0.5 * (a + b + c + d -
			  math.sqrt((a + b + c + d)^2 - 4 * (a * c + a * d + b * d)))
	 if math.abs(lambda1 - c - d) < math.abs(lambda2 - c - d) then
	    lambda = lambda1
	 else
	    lambda = lambda2
	 end
      end
      local y, z = matrix[1][1]^2 - lambda, matrix[1][1] * matrix[1][2]
      local r = Numerics.hypot(y, z)
      
      local c, s = y / r, -z / r
      for k = 1, n do
	 matrix[k][1], matrix[k][2] = c * matrix[k][1] - s * matrix[k][2],
	    s * matrix[k][1] + c * matrix[k][2]
      end
      for k = 1, v.width do
	 v[1][k], v[2][k] = c * v[1][k] - s * v[2][k], s * v[1][k] + c * v[2][k]
      end
      matrix, u, v = matrix:bidiagonalize(u, v, tol)
   end
   return matrix:getDiagonal(), u, v
end

function Matrix:tridiagonalize(que, tol)
   tol = tol or 10 ^ -13
   local matrix = self:copy()
   que = que or Matrix:identity(self.length)
   -- This is the top row
   for i = 2, self.length - 1 do
      local column = i - 1
      -- This is the bottom row
      for j = i + 1, self.length do
         local a = matrix[i][column]
         local b = matrix[j][column]
	 if(math.abs(b) <= tol) then
	    continue
	 end
         local r = Numerics.hypot(a, b)
         -- Anything less than tolerance is assumed to be zero,
	 -- so just set things to zero!
         if r <= tol then
            matrix[i][column], matrix[j][column], matrix[column][i], matrix[column][j]
	       = 0, 0, 0, 0
            continue
         end
         local c = a / r
         local s = -b / r
         matrix[i][column], matrix[j][column], matrix[column][i], matrix[column][j]
	    = r, 0, r, 0
	 for k = 1, self.width do
	    que[i][k], que[j][k]
	       = c * que[i][k] - s * que[j][k], s * que[i][k] + c * que[j][k]
	 end
	 -- This is the active column
         for k = column + 1, self.width do
            if i == k then
               local a1 = matrix[i][i]
               local a3 = matrix[j][j]
               local a2 = matrix[i][j]
               matrix[i][i] = c * (a1 * c - a2 * s) - s * (a2 * c - a3 * s)
               matrix[j][j] = s * (a2 * c + a1 * s) + c * (a3 * c + a2 * s)
               local kappa = s * (a1 * c - a2 * s) + c * (a2 * c - a3 * s)
               matrix[i][j], matrix[j][i] = kappa, kappa
            elseif j == k then
               continue
            else
               matrix[i][k], matrix[j][k]
		  = c * matrix[i][k] - s * matrix[j][k],
		  s * matrix[i][k] + c * matrix[j][k]
               matrix[k][i], matrix[k][j] = matrix[i][k], matrix[j][k]
            end
         end
      end
   end
   return matrix, que
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

function Matrix:symmetricFrancisOne(que, tol, maxIters)
   que = que or Matrix:identity(self.length)
   tol = tol or 10 ^ -7
   maxIters = maxIters or 1000
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }, Matrix:new({{1}})
   elseif self.length ~= self.width then
      matrix = self:padTo(math.max(self.length, self.width))
      self = matrix
   else
      matrix = self
   end
   local n = matrix.length
   if n == 2 then
      local temp1 = 4 * matrix[1][2] * matrix[2][1] +
	 math.pow((matrix[1][1] - matrix[2][2]), 2)
      if matrix[1][2] < tol then
	 return {matrix[1][1], matrix[2][2]}, Matrix:identity(2)
      else
	 local lambda1 = (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2
	 local lambda2 = (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2
	 local vect = {matrix[2][2] - lambda1, matrix[1][2]}
	 local norm = Numerics.hypot(vect[1], vect[2])
	 local c, d = vect[1] / norm, vect[2] / norm
	 return {lambda1, lambda2}, Matrix:new({{c, -d}, {d, c}})
      end
   end
   local iterations = 0
   local temp = 0
   local tempMin = 0
   local p = 1
   while iterations < maxIters do
      local min = 10 ^ 12
      local minimizer = 0
      for i = 1, n - 1 do
         if math.abs(matrix[i + 1][i]) < min then
            min = math.abs(matrix[i + 1][i])
            minimizer = i
         end
         if math.abs(matrix[i + 1][i]) < tol
            or (temp == minimizer and math.abs(tempMin - min) < tol) then
            local eig1, que1 =
               matrix:submatrix(1, i, 1, i):symmetricFrancisOne(Matrix:identity(i),
                                                                tol,
                                                                maxIters)
            local eig2, que2 =
               matrix:submatrix(i + 1, n,
				i + 1, n):symmetricFrancisOne(Matrix:identity(n - i),
							      tol,
							      maxIters)
            que:setSubmatrix(1, i, 1, n, que1 * que:submatrix(1, i, 1, n))
            que:setSubmatrix(i + 1, n, 1, n, que2 * que:submatrix(i + 1, n, 1, n))
            return Tools.list.join(eig1, eig2), que
         end
      end
      local shift = matrix[n][n]
      local u = {
         matrix[1][1] - shift + math.sqrt(math.pow(matrix[1][1] - shift, 2)
					  + math.pow(matrix[2][1], 2)),
         matrix[2][1],
      }
      local gamma = 2 / (math.pow(u[1], 2) + math.pow(u[2], 2))
      u = Matrix:new(u)
      local Q = Matrix:identity(2) - (u * u:transpose()):toScaled(gamma)
      for i = 1, n do
         matrix:setSubmatrix(1, p + 1, i, i, Q * matrix:submatrix(1, p + 1, i, i))
      end
      for i = 1, n do
         matrix:setSubmatrix(i, i, 1, p + 1,
			     matrix:submatrix(i, i, 1, p + 1) * Q:transpose())
      end
      que:setSubmatrix(1, p + 1, 1, n, Q * que:submatrix(1, p + 1, 1, n))
      matrix, que = matrix:tridiagonalize(que, tol)
      iterations = iterations + 1
      temp = minimizer
      tempMin = min
      if iterations == maxIters then
         local eig1, que1 =
            matrix:submatrix(1, minimizer, 1, minimizer)
	    :symmetricFrancisOne(Matrix:identity(minimizer),
				 tol,
				 maxIters)
         local eig2, que2 =
            matrix:submatrix(minimizer + 1, n, minimizer + 1, n):symmetricFrancisOne(Matrix:identity(n - minimizer),
                                                                                     tol,
                                                                                     maxIters)
         que:setSubmatrix(1,
			  minimizer,
			  1,
			  n,
			  que1 * que:submatrix(1, minimizer, 1, n))
         que:setSubmatrix(minimizer + 1,
			  n,
			  1,
			  n,
			  que2 * que:submatrix(minimizer + 1, n, 1, n))
         return Tools.list.join(eig1, eig2), que
      end
   end
   print("You shouldn't have gotten here")
   return matrix
end

function Matrix:symmetricFrancisTwo(que, tol, maxIters)
   tol = tol or 10 ^ -7
   maxIters = maxIters or 10000
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }, Matrix:new({{1}})
   elseif self.length ~= self.width then
      matrix = self:padTo(math.max(self.length, self.width))
      self = matrix
   else
      matrix = self
   end
   local n = matrix.length
   if n == 2 then
      local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
      if matrix[1][2] < tol then
	 return {matrix[1][1], matrix[2][2]}, Matrix:identity(2)
      else
	 local lambda1 = (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2
	 local lambda2 = (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2
	 local vect = {matrix[2][2] - lambda1, matrix[1][2]}
	 local norm = Numerics.hypot(vect[1], vect[2])
	 local c, d = vect[1] / norm, vect[2] / norm
	 return {lambda1, lambda2}, Matrix:new({{c, -d}, {d, c}})
      end
   elseif n == 3 then
      return self:symmetricFrancisOne(que, tol, maxIters)
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
         if math.abs(matrix[i + 1][i]) < tol
            or (temp == minimizer and math.abs(tempMin - min) < tol) then
            local eig1, que1 =
               matrix:submatrix(1, i, 1, i):symmetricFrancisTwo(Matrix:identity(i),
                                                                tol,
                                                                maxIters)
            local eig2, que2 =
               matrix:submatrix(i + 1, n, i + 1, n):symmetricFrancisTwo(Matrix:identity(n - i),
                                                                        tol,
                                                                        maxIters)
            que:setSubmatrix(1, i, 1, n, que1 * que:submatrix(1, i, 1, n))
            que:setSubmatrix(i + 1, n, 1, n, que2 * que:submatrix(i + 1, n, 1, n))
            return Tools.list.join(eig1, eig2), que
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
      que:setSubmatrix(1, p + 1, 1, n, Q * que:submatrix(1, p + 1, 1, n))
      matrix, que = matrix:tridiagonalize(que, tol)
      iterations = iterations + 1
      temp = minimizer
      tempMin = min
      if iterations == maxIters then
         local eig1, que1 =
            matrix:submatrix(1, minimizer, 1, minimizer):symmetricFrancisTwo(Matrix:identity(minimizer),
                                                                             tol,
                                                                             maxIters)
         local eig2, que2 =
            matrix:submatrix(minimizer + 1, n, minimizer + 1, n):symmetricFrancisTwo(Matrix:identity(n - minimizer),
                                                                                     tol,
                                                                                     maxIters)
         que:setSubmatrix(1,
			  minimizer,
			  1,
			  n,
			  que1 * que:submatrix(1, minimizer, 1, n))
         que:setSubmatrix(minimizer + 1,
			  n,
			  1,
			  n,
			  que2 * que:submatrix(minimizer + 1, n, 1, n))
         return Tools.list.join(eig1, eig2), que
      end
   end
   print("You shouldn't have gotten here")
   return matrix
end

function Matrix:symmetricFrancisSix(que, tol, maxIters)
   -- Set Defaults
   tol = tol or 10 ^ -13
   maxIters = maxIters or 10000
   -- Create Matrix
   local matrix
   if self.length == 1 and self.width == 1 then
      return { self[1][1] }, Matrix:new({{1}})
   elseif self.length ~= self.width then
      -- We will assume that the matrix was meant to be square
      matrix = self:padTo(math.max(self.length, self.width))
      self = matrix
   else
      matrix = self
   end
   -- Handle the case where n is too small
   local n = matrix.length
   if n == 2 then
      local temp1 = 4 * matrix[1][2] * matrix[2][1] + math.pow((matrix[1][1] - matrix[2][2]), 2)
      if matrix[1][2] < tol then
	 return {matrix[1][1], matrix[2][2]}, Matrix:identity(2)
      else
	 local lambda1 = (matrix[1][1] + matrix[2][2] + math.sqrt(temp1)) / 2
	 local lambda2 = (matrix[1][1] + matrix[2][2] - math.sqrt(temp1)) / 2
	 local vect = {matrix[2][2] - lambda1, matrix[1][2]}
	 local norm = Numerics.hypot(vect[1], vect[2])
	 local c, d = vect[1] / norm, vect[2] / norm
	 return {lambda1, lambda2}, Matrix:new({{c, -d}, {d, c}})
      end
   elseif n < 7 then
      return self:symmetricFrancisTwo(que, tol, maxIters)
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
         if math.abs(matrix[i + 1][i]) < tol or
            (temp == minimizer and math.abs(tempMin - min) < tol) then
            local eig1, que1 =
               matrix:submatrix(1, i, 1, i):symmetricFrancisSix(Matrix:identity(i),
                                                                tol,
                                                                maxIters)
            local eig2, que2 =
               matrix:submatrix(i + 1, n, i + 1, n):symmetricFrancisSix(Matrix:identity(n - i),
                                                                        tol,
                                                                        maxIters)
            que:setSubmatrix(1, i, 1, n, que1 * que:submatrix(1, i, 1, n))
            que:setSubmatrix(i + 1, n, 1, n, que2 * que:submatrix(i + 1, n, 1, n))
            return Tools.list.join(eig1, eig2), que
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
      que:setSubmatrix(1, p + 1, 1, n, Q * que:submatrix(1, p + 1, 1, n))
      matrix, que = matrix:tridiagonalize(que, tol)
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
         local eig1, que1 = matrix:submatrix(
            1,
            minimizer,
            1,
            minimizer):symmetricFrancisSix(Matrix:identity(minimizer), tol, maxIters)
         local eig2, que2 = matrix:submatrix(minimizer + 1,
                                             n,
                                             minimizer + 1,
                                             n):symmetricFrancisSix(Matrix:identity(minimizer), tol, maxIters)
         que:setSubmatrix(1,
                          minimizer,
                          1,
                          n,
                          que1 * que:submatrix(1, minimizer, 1, n))
         que:setSubmatrix(minimizer + 1,
                          n,
                          1,
                          n,
                          que2 * que:submatrix(minimizer + 1, n, 1, n))
         return Tools.list.join(eig1, eig2), que
      end
   end
   print("You shouldn't have gotten here")
   return matrix
end

function Matrix:eigenvalues()
   local matrix = self:hessenbergForm()
   local eigList = matrix:francisSix()
   Complex:sort(eigList, "d")
   return Tools.list.sublist(eigList, math.min(matrix.length, matrix.width))
end

function Matrix:symmetricEigendecomposition(tol)
   local matrix, que = self:tridiagonalize(tol)
   return matrix:symmetricFrancisSix(que, tol)
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
      local rowconstant = self.width * (i - 1)
      for j, w in pairs(table[i]) do
         if w == 0 then
            self.data[rowconstant + j] = nil
         else
            self.data[rowconstant + j] = w
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
   local rowconstant = i - self.width
   for j, w in pairs(table) do
      if w == 0 then
         self.data[self.width * j + rowconstant] = nil
      else
         self.data[self.width * j + rowconstant] = w
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
      local columnconstant = i - self.length
      local rowconstant = (i - 1) * self.width
      for j = 1, self.width do
         data[j * self.length + columnconstant] = self.data[rowconstant + j]
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
      self.data[i] = mu * v
   end
   return self
end

function SparseMatrix:scaled(mu)
   local matrix = self:copy()
   return matrix:scale(mu)
end

function SparseMatrix:kroneckerProduct(matrix)
   local result = SparseMatrix:new({}, self.length * matrix.length, self.width * matrix.width)
   for i = 1, self.length do
      for j = 1, self.width do
	 local scalar = self:get(i, j)
	 if scalar ~= 0 then
	    result:setSubmatrix((i - 1) * matrix.width + 1, (j - 1) * matrix.length + 1, matrix:scaled(scalar))
	 end
      end
   end
   return result
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
      local nr = r + rowStart - 2
      local nc = c + columnStart - 1
      local rc = self.width * nr + nc
      self.data[rc] = v
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

function SparseMatrix:conjugateGradient(b, x, tol)
   x = x or Vectors.zeros(A.width)
   local r = Vectors.sub(b, self:apply(x))
   local rr = r
   local p = r
   local res = Vectors.norm(r)
   while res > tol do
      local q = self:apply(p)
      local alpha = Vectors.dot(r, r) / Vectors.dot(p, q)
      x = Vectors.add(x, Vectors.scale(alpha, p))
      r = Vectors.sub(r, Vectors.scale(alpha, q))
      res = Vectors.norm(r)
      local beta = Vectors.dot(r, r) / Vectors.dot(rr, rr)
      rr = r
      p = Vectors.add(r, Vectors.scale(beta, p))
   end
   return x
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

--- WARNING EXPERIMENTAL!

local SparseSymmetricMatrix = {
   length = 0,
   data = nil
}

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
