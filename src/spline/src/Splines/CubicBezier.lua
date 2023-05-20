local CubicPolynomial = require(script.Parent.Parent.SplineMethods.CubicPolynomial)
local PositionSpline = require(script.Parent.Parent.SplineMethods.PositionSpline)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicBezier = {}

function CubicBezier.new(p0: Point, p1: Point, p2: Point, p3: Point)
	local self = {}

	self.a = p0
	self.b = 3 * (p1 - p0)
	self.c = 3 * (p0 + p2) - 6 * p1
	self.d = p3 + 3 * (p1 - p2) - p0

	-- Derivatives
	self.SolvePosition = CubicPolynomial.SolvePosition
	self.SolveVelocity = CubicPolynomial.SolveVelocity
	self.SolveAcceleration = CubicPolynomial.SolveAcceleration
	self.SolveJerk = CubicPolynomial.SolveJerk

	-- Frenet-Serret frame
	self.SolveTangent = PositionSpline.SolveTangent
	self.SolveNormal = PositionSpline.SolveNormal
	self.SolveBinormal = PositionSpline.SolveBinormal
	self.SolveCurvature = PositionSpline.SolveCurvature
	self.SolveTorsion = PositionSpline.SolveTorsion

	self.ToUnitSpeed = PositionSpline.ToUnitSpeed

	return self
end

return CubicBezier
