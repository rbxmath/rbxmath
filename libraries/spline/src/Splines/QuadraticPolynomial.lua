local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticPolynomial = {}

QuadraticPolynomial.Interpolant = setmetatable({}, Position.Interpolant)
QuadraticPolynomial.Interpolant.__index = QuadraticPolynomial.Interpolant

function QuadraticPolynomial.Interpolant.new(a: Vector, b: Vector, c: Vector)
	local self = setmetatable(Position.Interpolant.new(), QuadraticPolynomial.Interpolant)

	self.a = a
	self.b = b
	self.c = c

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function QuadraticPolynomial.Interpolant:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * self.c)
end

function QuadraticPolynomial.Interpolant:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * 2 * self.c
end

function QuadraticPolynomial.Interpolant:SolveAcceleration(): Vector
	-- r''(t)
	return 2 * self.c
end

function QuadraticPolynomial.Interpolant:SolveJerk(): Vector
	-- r'''(t)
	return 0 * self.a
end

return QuadraticPolynomial
