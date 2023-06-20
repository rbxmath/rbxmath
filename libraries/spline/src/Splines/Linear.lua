local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local Linear = {}

Linear.Interpolant = setmetatable({}, Position.Interpolant)
Linear.Interpolant.__index = Linear.Interpolant

Linear.Spline = setmetatable({}, Position.Spline)
Linear.Spline.__index = Linear.Spline

function Linear.Interpolant.new(p0: Point, p1: Point)
	local self = setmetatable(Position.Interpolant.new(), Linear.Interpolant)

	self.a = p0
	self.b = p1 - p0

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)
	self.Length = self.b.Magnitude

	return self
end

function Linear.Interpolant:SolvePosition(t: number): Point
	return self.a + self.b * t
end

function Linear.Interpolant:SolveVelocity(): Vector
	return self.b
end

function Linear.Interpolant:SolveAcceleration(): Vector
	return self.a * 0
end

function Linear.Interpolant:SolveJerk(): Vector
	return self.a * 0
end

function Linear.Interpolant:SolveTangent(): Vector
	return self.b.Unit
end

function Linear.Interpolant:SolveNormal(): Vector
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
		return self.a * 0
	end
end

function Linear.Interpolant:SolveBinormal(): Vector
	if self.Codimension == 2 then
		return self.a * 0
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function Linear.Interpolant:SolveCurvature(): number
	return 0
end

function Linear.Interpolant:SolveTorsion(): number
	return 0
end

function Linear.Spline.new(points: { Point })
	assert(#points > 0)

	if #points == 1 then
		local point = points[1]
		local self = setmetatable(Position.Spline.new(), Linear.Spline)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(point)
		self.Length = 0
		self.Interpolants = { Linear.Interpolant.new(point, point) }
		self.InterpolantDomains = { 0 }

		return self
	end

	local prevPoint = points[1]
	local interpolants = table.create(#points - 1)
	local totalLength = 0

	for i = 2, #points do
		local point = points[i]
		local interpolant = Linear.Interpolant.new(point, prevPoint)

		interpolants[i - 1] = interpolant
		totalLength += interpolant.Length
		prevPoint = point
	end

	local self = setmetatable(Position.Spline.new(), Linear.Spline)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(prevPoint)
	self.Length = totalLength
	self.Interpolants = interpolants
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return Linear
