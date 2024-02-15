local QuadraticPolynomial = require(script.Parent.QuadraticPolynomial)
local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticBezier = {}

QuadraticBezier.Interpolant = setmetatable({}, QuadraticPolynomial.Interpolant)
QuadraticBezier.Interpolant.__index = QuadraticBezier.Interpolant

QuadraticBezier.Spline = setmetatable({}, Position.Spline)
QuadraticBezier.Spline.__index = QuadraticBezier.Spline

function QuadraticBezier.Interpolant.new(p0: Point, p1: Point, p2: Point)
	local a = p0
	local b = 2 * (p1 - p0)
	local c = p0 - 2 * p1 + p2

	return setmetatable(QuadraticPolynomial.Interpolant.new(a, b, c), QuadraticBezier.Interpolant)
end

function QuadraticBezier.Spline.new(points: { Point })
	assert(#points % 2 == 1 and #points >= 3, "Needs 2n+1 points for n >= 1")

	local numInterpolants = (#points - 1) / 2
	local interpolants = table.create(numInterpolants)
	local totalLength = 0

	for i = 1, numInterpolants do
		local startIndex = (i - 1) * 2
		local interpolant =
			QuadraticBezier.Interpolant.new(points[startIndex + 1], points[startIndex + 2], points[startIndex + 3])

		interpolants[i] = interpolant
		totalLength += interpolant.Length
	end

	local self = setmetatable(Position.Spline.new(), QuadraticBezier.Spline)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(points[1])
	self.Length = totalLength
	self.Interpolants = interpolants
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return QuadraticBezier
