local PositionSpline = require(script.Parent.PositionSpline)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicPolynomial = setmetatable({}, PositionSpline)
CubicPolynomial.__index = CubicPolynomial

function CubicPolynomial.new(a: Vector, b: Vector, c: Vector, d: Vector)
	local self = setmetatable(PositionSpline.new(), CubicPolynomial)

	self.a = a
	self.b = b
	self.c = c
	self.d = d

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function CubicPolynomial:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * (self.c + t * self.d))
end

function CubicPolynomial:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * (2 * self.c + t * 3 * self.d)
end

function CubicPolynomial:SolveAcceleration(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r''(t)
	return 2 * self.c + t * 6 * self.d
end

function CubicPolynomial:SolveJerk(): Vector
	-- r'''(t)
	return 6 * self.d
end

return CubicPolynomial
