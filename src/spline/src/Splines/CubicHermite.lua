local CubicPolynomial = require(script.Parent.CubicPolynomial)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicHermite = setmetatable({}, CubicPolynomial)
CubicHermite.__index = CubicHermite

function CubicHermite.new(p0: Point, p1: Point, v0: Vector, v1: Vector)
	local a = p0
	local b = v0
	local c = 3 * (p1 - p0) - 2 * v0 - v1
	local d = v1 + v0 + 2 * (p0 - p1)

	return setmetatable(CubicPolynomial.new(a, b, c, d), CubicHermite)
end

return CubicHermite
