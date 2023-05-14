local CubicPolynomial = require(script.Parent.CubicPolynomial)
local SplineUtils = require(script.Parent.SplineUtils)
local Types = require(script.Parent.Types)

local CubicBezier = setmetatable({}, CubicPolynomial)
CubicBezier.__index = CubicBezier

function CubicBezier.new(points: { Point }, closed: boolean?)
	local self = setmetatable(CubicPolynomial.new(), CubicBezier)
	self._closed = closed
	self._codimension = SplineUtils.GetCodimensionFromPoint(points[1])

	local numCurves = (#points - 1) / 3
	assert(numCurves % 1 == 0, "Needs 3n + 1 points, for some positive integer n")

	local curves = table.create(numCurves)
	local p0 = points[1]

	for i = 1, numCurves do
		local idx = 3 * (i - 1) + 1
		local p1 = points[i + 1]
		local p2 = points[i + 2]
		local p3 = points[i + 3]

		curves[i] = {
			a = p0,
			b = 3 * (p1 - p0),
			c = 3 * (p0 + p2) - 6 * p1,
			d = p3 + 3 * (p1 - p2) - p0,
		}

		p0 = p3
	end

	self.Curves = curves
	self._numCurves = numCurves
	self:_cacheLengths()

	return self
end

return CubicBezier
