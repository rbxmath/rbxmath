local QuadraticPolynomial = require(script.Parent.Parent.SplineMethods.QuadraticPolynomial)
local PositionSpline = require(script.Parent.Parent.SplineMethods.PositionSpline)
local Types = require(script.Parent.Parent.Types)

local QuadraticBezier = {}

function QuadraticBezier.new(p0: Point, p1: Point, p2: Point)
	local self = {}

	self.a = p0
	self.b = 2 * (p1 - p0)
	self.c = p0 - 2 * p1 + p2

	-- Derivatives
	self.SolvePosition = QuadraticPolynomial.SolvePosition
	self.SolveVelocity = QuadraticPolynomial.SolveVelocity
	self.SolveAcceleration = QuadraticPolynomial.SolveAcceleration
	self.SolveJerk = QuadraticPolynomial.SolveJerk

	-- Frenet-Serret frame
	self.SolveTangent = PositionSpline.SolveTangent
	self.SolveNormal = PositionSpline.SolveNormal
	self.SolveBinormal = PositionSpline.SolveBinormal
	self.SolveCurvature = PositionSpline.SolveCurvature
	self.SolveTorsion = PositionSpline.SolveTorsion

	self.ToUnitSpeed = PositionSpline.ToUnitSpeed

	return self
end

return QuadraticBezier
