--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Scalars = require(script.Parent.Scalars)
local Complex = Scalars.Complex
local Tools = require(script.Parent.Tools)
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object

local Numerics = {}

function Numerics.hypot(a: number, b: number)
   a, b = math.abs(a), math.abs(b)
   if a < b then
      a, b = b, a
   end
   return a * math.sqrt(1 + (b / a) ^ 2)
end

function Numerics.kahanHypot(...)
   local list = {...}
   local max = 0
   local ind = 0
   for i = 1, #list do
      list[i] = math.abs(list[i])
      if list[i] > max then
	 max = list[i]
	 ind = i
      end
   end
   local oneOverMax = 1 / max
   local list2 = {}
   for i = 1, #list do
      if i ~= ind then
	 list2[#list2 + 1] = (list[i] * oneOverMax)^2
      end
   end
   return math.sqrt(Numerics.kahan(table.unpack(list2)))
end

function Numerics.quadratic(a, b, c, disc)
   -- The user may have already computed the discriminant
   -- in order to figure out if the roots are real. The whole
   -- point of this library is to have fast, efficient, and stable
   -- methods for solving problems, so we will allow the user to pass
   -- the discriminant in as well.
   disc = disc or math.sqrt(b * b - 4 * a * c)
   local twoc = 2 * c
   if a == 0 then
      if b == 0 then
	 error("Quadratic isn't solvable")
      else
	 return -c / b
      end
   end
   if disk >= 0 then
      if b > 0 then
	 return -twoc / (b + disc), twoc / (-b + disc)
      else
	 return twoc / (-b - disc), -twoc / (b + disc)
      end
   else
      local overTwoA = 1 / (2 * a)
      local bOverTwoA = b * overTwoA
      local discOverTwoA = disc * overTwoA
      return Complex:new(bOverTwoA, discOverTwoA),
	 Complex:new(bOverTwoA, -discOverTwoA)
   end
end

-- Kahan implementation taken from Wikipedia
-- https://en.wikipedia.org/wiki/Kahan_summation_algorithm
function Numerics.kahan(...)
   local sum = 0
   local c = 0
   local list = {...}
   for k, v in ipairs(list) do
      local y = v - c
      local t = sum + y
      c = (t - sum) - y
      sum = t
   end
   return sum
end

return Numerics
