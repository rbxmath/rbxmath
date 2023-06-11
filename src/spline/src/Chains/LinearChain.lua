local Linear = require(script.Parent.Parent.Splines.Linear)
local PositionSplineChain = require(script.Parent.PositionSplineChain)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point

local LinearChain = setmetatable({}, PositionSplineChain)
LinearChain.__index = LinearChain

function LinearChain.new(points: { Point })
	assert(#points > 0)

	if #points == 1 then
		local point = points[1]
		local self = setmetatable(PositionSplineChain.new(), LinearChain)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(point)
		self.Length = 0
		self.Splines = { Linear.new(point, point) }
		self.SplineDomains = { 0 }

		return self
	end

	local prevPoint = points[1]
	local splines = table.create(#points - 1)
	local totalLength = 0

	for i = 2, #points do
		local point = points[i]
		local spline = Linear.new(point, prevPoint)

		splines[i - 1] = spline
		totalLength += spline.Length
		prevPoint = point
	end

	local self = setmetatable(PositionSplineChain.new(), LinearChain)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(prevPoint)
	self.Length = totalLength
	self.Splines = splines
	self.SplineDomains = SplineUtils.GetSplineDomains(splines)

	return self
end
