local PositionSpline = require(script.Parent.PositionSpline)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local Linear = setmetatable({}, PositionSpline)
Linear.__index = Linear

function Linear.new(p0: Point, p1: Point)
	local self = setmetatable(PositionSpline.new(), Linear)

	self.p0 = p0
	self.p1 = p1
	self.p1_p0 = p1 - p0

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)
	self.Length = self.p1_p0.Magnitude

	return self
end

function Linear:SolvePosition(t: number): Point
	return self.p0 + t * self.p1_p0
end

function Linear:SolveVelocity(): Vector
	return self.p1_p0
end

function Linear:SolveAcceleration(): Vector
	return self.p0 * 0
end

function Linear:SolveJerk(): Vector
	return self.p0 * 0
end

function Linear:SolveTangent(): Vector
	return self.p1_p0.Unit
end

function Linear:SolveNormal(): Vector
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

function Linear:SolveBinormal(): Vector
	if self.Codimension == 2 then
		return self.p0 * 0
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function Linear:SolveCurvature(): number
	return 0
end

function Linear:SolveTorsion(): number
	return 0
end

return Linear
