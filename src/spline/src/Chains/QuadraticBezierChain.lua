local QuadraticBezier = require(script.Parent.Parent.Splines.QuadraticBezier)
local PositionSplineChain = require(script.Parent.PositionSplineChain)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point

local QuadraticBezierChain = setmetatable({}, PositionSplineChain)
QuadraticBezierChain.__index = QuadraticBezierChain

function QuadraticBezierChain.new(points: { Point })
	assert(#points % 2 == 1 and #points >= 3)

	local numSplines = (#points - 1) / 2
	local splines = table.create(numSplines)
	local totalLength = 0

	for i = 1, numSplines do
		local startIndex = (i - 1) * 2
		local spline = QuadraticBezier.new(points[startIndex + 1], points[startIndex + 2], points[startIndex + 3])

		splines[i] = spline
		totalLength += spline.Length
	end

	local self = setmetatable(PositionSplineChain.new(), QuadraticBezierChain)

	self.Length = totalLength
	self.Splines = splines
	self.SplineDomains = SplineUtils.GetSplineDomains(splines)

	return self
end

return QuadraticBezierChain
