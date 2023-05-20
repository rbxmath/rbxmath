local CubicPolynomial = require(script.Parent.Parent.SplineMethods.CubicPolynomial)
local PositionSpline = require(script.Parent.Parent.SplineMethods.PositionSpline)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicHermite = {}

function CubicHermite.new(p0: Point, p1: Point, v0: Vector, v1: Vector)
	local self = {}

	self.a = p0
	self.b = v0
	self.c = 3 * (p1 - p0) - 2 * v0 - v1
	self.d = v1 + v0 + 2 * (p0 - p1)

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

return CubicHermite
