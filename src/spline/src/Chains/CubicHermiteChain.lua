local CubicHermite = require(script.Parent.Parent.Splines.CubicHermite)
local PositionSplineChain = require(script.Parent.PositionSplineChain)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicHermiteChain = setmetatable({}, PositionSplineChain)
CubicHermiteChain.__index = CubicHermiteChain

function CubicHermiteChain.new(points: { Point }, vectors: { Vector })
	assert(#points % 2 == 0 and #points >= 2)
	assert(#points == #vectors)

	local numSplines = #points - 1
	local splines = table.create(numSplines)
	local totalLength = 0

	for i = 1, numSplines do
		local startIndex = (i - 1) * 2
		local spline = CubicHermite.new(
			points[startIndex + 1],
			points[startIndex + 2],
			vectors[startIndex + 1],
			vectors[startIndex + 2]
		)

		splines[i] = spline
		totalLength += spline.Length
	end

	local self = setmetatable(PositionSplineChain.new(), CubicHermiteChain)

	self.Length = totalLength
	self.Splines = splines
	self.SplineDomains = SplineUtils.GetSplineDomains(splines)

	return self
end

return CubicHermiteChain
