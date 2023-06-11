local CubicBezier = require(script.Parent.Parent.Splines.CubicBezier)
local PositionSplineChain = require(script.Parent.PositionSplineChain)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point

local CubicBezierChain = setmetatable({}, PositionSplineChain)
CubicBezierChain.__index = CubicBezierChain

function CubicBezierChain.new(points: { Point })
	assert(#points % 3 == 1 and #points >= 4)

	local numSplines = (#points - 1) / 3
	local splines = table.create(numSplines)
	local totalLength = 0

	for i = 1, numSplines do
		local startIndex = (i - 1) * 3
		local spline = CubicBezier.new(
			points[startIndex + 1],
			points[startIndex + 2],
			points[startIndex + 3],
			points[startIndex + 4]
		)

		splines[i] = spline
		totalLength += spline.Length
	end

	local self = setmetatable(PositionSplineChain.new(), CubicBezierChain)

	self.Length = totalLength
	self.Splines = splines
	self.SplineDomains = SplineUtils.GetSplineDomains(splines)

	return self
end

return CubicBezierChain
