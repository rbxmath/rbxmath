local QuadraticPolynomial = require(script.Parent.QuadraticPolynomial)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticBezier = setmetatable({}, QuadraticPolynomial)
QuadraticBezier.__index = QuadraticBezier

function QuadraticBezier.new(p0: Point, p1: Point, p2: Point)
	local a = p0
	local b = 2 * (p1 - p0)
	local c = p0 - 2 * p1 + p2

	return setmetatable(QuadraticPolynomial.new(a, b, c), QuadraticPolynomial)
end

return QuadraticBezier
