--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

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

return Numerics
