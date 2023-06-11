--[[
	A generic spline class with spline-independent methods
	
	Notation
	r(t): Position
	T(t): Tangent vector
	N(t): Normal vector
	B(t): Binormal vector
	κ(t): Curvature
	τ(t): Torsion
]]

local Interpolation = require(script.Parent.Parent.Vendor.Interpolation)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local PositionSpline = {}
PositionSpline.__index = PositionSpline

function PositionSpline.new()
	local self = setmetatable({}, PositionSpline)

	self._chebyshevInterpolant = nil
	self.Codimension = nil
	self.IsUnitSpeed = false
	self.Length = nil

	return self
end

function PositionSpline:_accountForUnitSpeed(t: number)
	return self.IsUnitSpeed and self._chebyshevInterpolant:evaluate(t * self.Length) or t
end

function PositionSpline:ToUnitSpeed()
	local numGridPoints = 16
	local grid = Interpolation.Chebyshev.grid(numGridPoints - 1)
	local linearRescalingFunction = Interpolation.Chebyshev.linearRescalingFunction(0, 1)
	local shiftedGrid = table.create(numGridPoints)

	-- Rescale Chebyshev grid from [-1, 1] -> [0, 1]
	for i = 1, numGridPoints do
		shiftedGrid[i] = linearRescalingFunction(grid[i])
	end

	-- Compute arc length at Chebyshev grid points
	local gridValues = table.create(numGridPoints)
	gridValues[1] = 0

	for i = 1, numGridPoints - 1 do
		gridValues[i + 1] = gridValues[i]
			+ SplineUtils.GaussLegendre(function(t)
				return self:SolveVelocity(t).Magnitude
			end, shiftedGrid[i], shiftedGrid[i + 1])
	end

	local interpolant = Interpolation.ChebyshevInterpolant:new(gridValues, 0, 1, numGridPoints - 1)
	interpolant.solveMethod = "Monotone"
	interpolant = interpolant:inverse()

	self.IsUnitSpeed = true
	self.Length = gridValues[numGridPoints]
	self._chebyshevInterpolant = interpolant
end

function PositionSpline:SolveTangent(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- T(t) = r'(t) / |r'(t)|
	return self:SolveVelocity(t).Unit
end

function PositionSpline:SolveNormal(t: number): Vector
	if self.Codimension == 0 then
		error("SolveNormal is restricted from splines in 1 dimension")
	else
		if self.Codimension == 1 then
			local tangent = self:SolveTangent()
			if typeof(tangent) == "Vector2" then
				return Vector2.new(-tangent.Y, tangent.X)
			else
				return Vector.new({ -tangent[2], tangent[1] })
			end
		else
			t = self:_accountForUnitSpeed(t)

			-- N(t) = T'(t) / |T'(t)|
			-- The return is equivalent to N(t) when the derivatives are carried out.
			-- In particular, the vector being unitized is T'(t) * |r'(t)| ^ 3, but
			-- the |r'(t)| ^ 3 scaling doesn't affect the result because we unitize it
			-- anyway. This scaled version is faster to compute.
			local vel = self:SolveVelocity(t)
			local accel = self:SolveAcceleration(t)
			local speed = vel.Magnitude

			return (accel * speed ^ 2 - vel * accel:Dot(vel)).Unit
		end
	end
end

function PositionSpline:SolveBinormal(t: number): Vector
	if self.Codimension == 2 then
		t = self:_accountForUnitSpeed(t)

		-- T(t) x N(t)
		return self:SolveTangent(t):Cross(self:SolveNormal(t))
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function PositionSpline:SolveCurvature(t: number): number
	t = self:_accountForUnitSpeed(t)

	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local speed = vel.Magnitude
	local dTangent = accel / speed - vel * vel:Dot(accel) / speed ^ 3

	-- κ(t) = |T'(t)| / |r'(t)|
	return dTangent.Magnitude / speed
end

function PositionSpline:SolveTorsion(t: number): number
	t = self:_accountForUnitSpeed(t)

	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local jerk = self:SolveJerk(t)
	local cross = vel:Cross(accel)

	-- τ(t) = ((r'(t) x r''(t)) • r'''(t)) / |r'(t) x r''(t)|^2
	return cross:Dot(jerk) / cross.Magnitude ^ 2
end

function PositionSpline:SolveLength(from: number?, to: number?): number
	from = from or 0
	to = to or 1

	if from > to then
		from, to = to, from
	end

	return SplineUtils.GaussLegendre(function(t)
		return self:SolveVelocity(t).Magnitude
	end, from, to)
end

return PositionSpline
