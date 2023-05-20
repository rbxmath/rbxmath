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

local GaussLegendre = require(script.Parent.Parent.GaussLegendre)
local Types = require(script.Parent.Parent.Types)
local Vector = require(script.Parent.Parent.Vector)

type Point = Types.Point
type Vector = Types.Vector

local PositionSpline = {}

function PositionSpline.ToUnitSpeed(self) end

function PositionSpline.SolveTangent(self, t: number): Vector
	-- T(t) = r'(t) / |r'(t)|
	return self:SolveVelocity(t).Unit
end

function PositionSpline.SolveNormal(self, t: number): Vector
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

function PositionSpline.SolveBinormal(self, t: number): Vector
	if self.Codimension == 2 then
		-- T(t) x N(t)
		return self:SolveTangent(t):Cross(self:SolveNormal(t))
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function PositionSpline.SolveCurvature(self, t: number): number
	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local speed = vel.Magnitude
	local dTangent = accel / speed - vel * vel:Dot(accel) / speed ^ 3

	-- κ(t) = |T'(t)| / |r'(t)|
	return dTangent.Magnitude / speed
end

function PositionSpline.SolveTorsion(self, t: number): number
	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local jerk = self:SolveJerk(t)
	local cross = vel:Cross(accel)

	-- τ(t) = ((r'(t) x r''(t)) • r'''(t)) / |r'(t) x r''(t)|^2
	return cross:Dot(jerk) / cross.Magnitude ^ 2
end

function PositionSpline.SolveLength(self, a: number?, b: number?): number
	a = a or 0
	b = b or 1

	if a > b then
		a, b = b, a
	end

	return GaussLegendre.Ten(function(x)
		return self:SolveVelocity(x).Magnitude
	end, a, b)
end

return PositionSpline
