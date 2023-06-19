local CubicPolynomial = require(script.Parent.CubicPolynomial)
local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicBezier = {}

CubicBezier.Interpolant = setmetatable({}, CubicPolynomial.Interpolant)
CubicBezier.Interpolant.__index = CubicBezier.Interpolant

CubicBezier.Spline = setmetatable({}, Position.Spline)
CubicBezier.Spline.__index = CubicBezier.Spline

function CubicBezier.Interpolant.new(p0: Point, p1: Point, p2: Point, p3: Point)
	local a = p0
	local b = 3 * (p1 - p0)
	local c = 3 * (p0 + p2) - 6 * p1
	local d = p3 + 3 * (p1 - p2) - p0

	return setmetatable(CubicPolynomial.Interpolant.new(a, b, c, d), CubicBezier.Interpolant)
end

function CubicBezier.Spline.new(points: { Point })
	assert(#points % 3 == 1 and #points >= 4, "Needs 3n+1 points for n >= 1")

	local numInterpolants = (#points - 1) / 3
	local interpolants = table.create(numInterpolants)
	local totalLength = 0

	for i = 1, numInterpolants do
		local startIndex = (i - 1) * 3
		local interpolant = CubicBezier.Interpolant.new(
			points[startIndex + 1],
			points[startIndex + 2],
			points[startIndex + 3],
			points[startIndex + 4]
		)

		interpolants[i] = interpolant
		totalLength += interpolant.Length
	end

	local self = setmetatable(Position.Spline.new(), CubicBezier.Spline)

	self.Length = totalLength
	self.Interpolants = interpolants
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return CubicBezier
