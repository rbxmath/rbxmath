local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicPolynomial = {}

CubicPolynomial.Interpolant = setmetatable({}, Position.Interpolant)
CubicPolynomial.Interpolant.__index = CubicPolynomial.Interpolant

function CubicPolynomial.Interpolant.new(a: Vector, b: Vector, c: Vector, d: Vector)
	local self = setmetatable(Position.Spline.new(), CubicPolynomial.Interpolant)

	self.a = a
	self.b = b
	self.c = c
	self.d = d

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function CubicPolynomial.Interpolant:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * (self.c + t * self.d))
end

function CubicPolynomial.Interpolant:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * (2 * self.c + t * 3 * self.d)
end

function CubicPolynomial.Interpolant:SolveAcceleration(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r''(t)
	return 2 * self.c + t * 6 * self.d
end

function CubicPolynomial.Interpolant:SolveJerk(): Vector
	-- r'''(t)
	return 6 * self.d
end

return CubicPolynomial
