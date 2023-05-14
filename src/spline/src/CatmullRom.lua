local PositionSpline = require(script.Parent.PositionSpline)
local SplineUtils = require(script.Parent.SplineUtils)
local Types = require(script.Parent.Types)

local CatmullRom = setmetatable({}, PositionSpline)
CatmullRom.__index = CatmullRom

function CatmullRom.new(points: { Point }) end

return CatmullRom
