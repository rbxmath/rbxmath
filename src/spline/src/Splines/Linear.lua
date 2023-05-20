local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)
local Vector = require(script.Parent.Parent.Vector)

type Point = Types.Point
type Vector = Types.Vector

local Linear = {}

local function SolvePosition(self, t: number): Point
	return self.p0 + t * self.p1_p0
end

local function SolveVelocity(self): Vector
	return self.p1_p0
end

local function SolveAcceleration(self): Vector
	return self.p0 * 0
end

local function SolveJerk(self): Vector
	return self.p0 * 0
end

local function SolveTangent(self): Vector
	return self.p1_p0.Unit
end

local function SolveNormal(self): Vector
	if self.Codimension == 0 then
		error("SolveNormal is restricted from splines in 1 dimension")
	elseif self.Codimension == 1 then
		local tangent = self:SolveTangent()
		if typeof(tangent) == "Vector2" then
			return Vector2.new(-tangent.Y, tangent.X)
		else
			return Vector.new({ -tangent[2], tangent[1] })
		end
	else
		return self.p0 * 0
	end
end

local function SolveBinormal(self): Vector
	if self.Codimension == 2 then
		return self.p0 * 0
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

local function SolveCurvature(): number
	return 0
end

local function SolveTorsion(): number
	return 0
end

function Linear.new(p0: Point, p1: Point)
	local self = {}

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)

	self.p0 = p0
	self.p1 = p1
	self.p1_p0 = p1 - p0

	-- Derivatives
	self.SolvePosition = SolvePosition
	self.SolveVelocity = SolveVelocity
	self.SolveAcceleration = SolveAcceleration
	self.SolveJerk = SolveJerk

	-- Frenet-Serret frame
	self.SolveTangent = SolveTangent
	self.SolveNormal = SolveNormal
	self.SolveBinormal = SolveBinormal
	self.SolveCurvature = SolveCurvature
	self.SolveTorsion = SolveTorsion

	self.ToUnitSpeed = ToUnitSpeed

	return self
end

return Linear
