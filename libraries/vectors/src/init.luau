--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require(script.Parent.Tools)
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarFunction = Tools.ScalarFunction
type Tensor = Tools.Tensor
type Object = Tools.Object
local Scalars = require(script.Parent.Scalars)
local Complex = Scalars.Complex

local Vectors = {}

function Vectors.copy(vector: Vector): Vector
   local data = {}
   for k, v in ipairs(vector) do
      data[k] = v
   end
   return data
end

function Vectors.zeros(n)
   local data = {}
   for i = 1, n do
      data[i] = 0
   end
   return data
end

function Vectors.dot(left: Vector, right: Vector): number
   if #left ~= #right then
      error("Attempting to compute the dot product of incompatible vectors!", 1)
   end
   local sum = 0
   for i = 1, #left do
      sum += left[i] * right[i]
   end
   return sum
end

function Vectors.scale(scalar: number, vector: Vector): Vector
   local data = {}
   for i = 1, #vector do
      data[i] = scalar * vector[i]
   end
   return data
end

function Vectors.norm(vector: Vector, norm: number | string)
   if not norm or norm == 2 then
      local sumLarge = 0
      local sumSmall = 0
      for i = 1, #vector do
	 local temp = vector[i] * vector[i]
	 if temp < sumSmall or (sumSmall == 0 and temp < 1) then
	    sumSmall += temp
	 else
	    sumLarge += temp
	 end
      end
      return math.sqrt(sumLarge + sumSmall)
   elseif norm == "infinity" then
      local max = 0
      for i = 1, #vector do
	 local maxCandidate = math.abs(vector[i])
	 if maxCandidate > max then
	    max = maxCandidate
	 end
      end
      return max
   elseif type(norm) == "number" then
      local sumLarge = 0
      local sumSmall = 0
      for i = 1, #vector do
	 temp = math.abs(vertor[i]) ^ norm
	 if temp < sumSmall or (sumSmall == 0 and temp < 1) then
	    sumSmall += temp
	 else
	    sumLarge += temp
	 end
      end
      return math.pow(sum, 1 / norm)
   else
      error("Improper norm supplied!")
   end
end

function Vectors.project(v, u, tol)
   tol = tol or 10 ^ -13
   local denominator = Vectors.dot(u, u)
   local numerator = Vectors.dot(u, v)
   if math.abs(denominator) < tol then
      return Vectors.scale(0, u)
   else
      return Vectors.scale(numerator / denominator, u)
   end
end

function Vectors.gramSchmidt(listOfVectors, tol)
   tol = tol or 10 ^ -13
   local listOfU = {}
   for i = 1, #listOfVectors do
      a = listOfVectors[i]
      u = Vectors.copy(a)
      for j = 1, i - 1 do
	 u = Vectors.sub(u, Vectors.project(a, listOfU[j]))
      end
      listOfU[i] = u
   end
   local basis = {}
   for i = 1, #listOfVectors do
      local scalar = Vectors.norm(listOfU[i])
      if scalar >= tol then
	 basis[#basis + 1] = Vectors.scale(1 / scalar, listOfU[i])
      end
   end
   return basis
end

function Vectors.randomVector(n: number, normalize: boolean | number | string): Vector
   local vector = {}
   for i = 1, n do
      local rand = math.random()
      vector[i] = rand
   end
   if normalize then
      if normalize == "infinity" then
	 local max = math.max(table.unpack(vector))
	 for i = 1, n do
	    vector[i] /= max
	 end
      elseif type(normalize) == "number" then
	 local sum = 0
	 for i = 1, n do
	    sum += math.abs(vector[i]) ^ normalize
	 end
	 sum = math.pow(sum, 1 / normalize)
	 for i = 1, n do
	    vector[i] /= sum
	 end
      end
   end

   return vector
end

function Vectors.add(left: Vector, right: Vector)
   if #left ~= #right then
      error("Attempting to subtract incompatible vectors!", 1)
   end

   local data = {}
   for i = 1, #left do
      data[i] = left[i] + right[i]
   end
   return data
end

function Vectors.sub(left: Vector, right: Vector)
   if #left ~= #right then
      error("Attempting to subtract incompatible vectors!", 1)
   end

   local data = {}
   for i = 1, #left do
      data[i] = left[i] - right[i]
   end
   return data
end

function Vectors.standardBasis(n: number, i: number)
   local data = {}
   for j = 1, n do
      if i == j then
	 data[j] = 1
      else
	 data[j] = 0
      end
   end
   return data
end

return Vectors
