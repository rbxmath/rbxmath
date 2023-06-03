local PositionSpline = require(script.Parent.PositionSpline)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticPolynomial = setmetatable({}, PositionSpline)
QuadraticPolynomial.__index = QuadraticPolynomial

function QuadraticPolynomial.new(a: Vector, b: Vector, c: Vector)
	local self = setmetatable(PositionSpline.new(), QuadraticPolynomial)

	self.a = a
	self.b = b
	self.c = c

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function QuadraticPolynomial:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * self.c)
end

function QuadraticPolynomial:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * 2 * self.c
end

function QuadraticPolynomial:SolveAcceleration(): Vector
	-- r''(t)
	return 2 * self.c
end

function QuadraticPolynomial:SolveJerk(): Vector
	-- r'''(t)
	return 0 * self.a
end

return QuadraticPolynomial
