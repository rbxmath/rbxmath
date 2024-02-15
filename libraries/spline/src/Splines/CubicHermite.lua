local CubicPolynomial = require(script.Parent.CubicPolynomial)
local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicHermite = {}

CubicHermite.Interpolant = setmetatable({}, CubicPolynomial.Interpolant)
CubicHermite.Interpolant.__index = CubicHermite.Interpolant

CubicHermite.Spline = setmetatable({}, Position.Spline)
CubicHermite.Spline.__index = CubicHermite.Spline

function CubicHermite.Interpolant.new(p0: Point, p1: Point, v0: Vector, v1: Vector)
	local a = p0
	local b = v0
	local c = 3 * (p1 - p0) - 2 * v0 - v1
	local d = v1 + v0 + 2 * (p0 - p1)

	return setmetatable(CubicPolynomial.Interpolant.new(a, b, c, d), CubicHermite.Interpolant)
end

function CubicHermite.Spline.new(points: { Point }, vectors: { Vector })
	assert(#points == #vectors, "Needs equal number of points and vectors")

	local numInterpolants = #points - 1
	local interpolants = table.create(numInterpolants)
	local totalLength = 0

	local prevPoint = points[1]
	local prevVector = vectors[1]

	for i = 1, numInterpolants do
		local point = points[i + 1]
		local vector = vectors[i + 1]
		local interpolant = CubicHermite.Interpolant.new(prevPoint, point, prevVector, vector)

		interpolants[i] = interpolant
		totalLength += interpolant.Length

		prevPoint = point
		prevVector = vector
	end

	local self = setmetatable(Position.Spline.new(), CubicHermite.Spline)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(prevPoint)
	self.Length = totalLength
	self.Interpolants = interpolants
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return CubicHermite
