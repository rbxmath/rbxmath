local CubicPolynomial = require(script.Parent.Parent.SplineMethods.CubicPolynomial)
local PositionSpline = require(script.Parent.Parent.SplineMethods.PositionSpline)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local DEFAULT_ALPHA = 0.5
local DEFAULT_TENSION = 0

local CatmullRom = {}

function CatmullRom.new(p0: Point, p1: Point, p2: Point, p3: Point, alpha: number?, tension: number?)
	alpha = alpha or DEFAULT_ALPHA
	tension = tension or DEFAULT_TENSION
	assert(alpha)
	assert(tension)

	-- https://qroph.github.io/2018/07/30/smooth-paths-using-catmull-rom-splines.html
	local p1_p0 = p1 - p0
	local p2_p1 = p2 - p1
	local p3_p2 = p3 - p2

	local t0 = 0
	local t1 = p1_p0.Magnitude ^ alpha + t0
	local t2 = p2_p1.Magnitude ^ alpha + t1
	local t3 = p3_p2.Magnitude ^ alpha + t2

	local scalar = (1 - tension) * (t2 - t1)
	local m1 = scalar * (p1_p0 / (t1 - t0) - (p2 - p0) / (t2 - t0) + p2_p1 / (t2 - t1))
	local m2 = scalar * (p2_p1 / (t2 - t1) - (p3 - p1) / (t3 - t1) + p3_p2 / (t3 - t2))

	local self = {}

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)

	self.a = p1
	self.b = m1
	self.c = 3 * p2_p1 - 2 * m1 - m2
	self.d = -2 * p2_p1 + m1 + m2

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

return CatmullRom
