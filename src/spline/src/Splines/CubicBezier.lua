local CubicPolynomial = require(script.Parent.CubicPolynomial)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicBezier = setmetatable({}, CubicPolynomial)
CubicBezier.__index = CubicBezier

function CubicBezier.new(p0: Point, p1: Point, p2: Point, p3: Point)
	local a = p0
	local b = 3 * (p1 - p0)
	local c = 3 * (p0 + p2) - 6 * p1
	local d = p3 + 3 * (p1 - p2) - p0

	local self = setmetatable(CubicPolynomial.new(a, b, c, d), CubicBezier)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)

	return self
end

return CubicBezier
