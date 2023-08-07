local Position = require(script.Parent.Position)
local Linear = require(script.Parent.Linear)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local Catenary = {}

Catenary.Interpolant = setmetatable({}, Position.Interpolant)
Catenary.Interpolant.__index = Catenary.Interpolant

Catenary.Spline = setmetatable({}, Position.Spline)
Catenary.Spline.__index = Catenary.Spline

local function solveRootUpperBound(hDist: number, c: number): number
	-- Upper bound on root of 8th order Maclaurin polynomial of f(z)
	return 2 * math.pow(-362880 * (hDist - c) / hDist, 1 / 8)
end

local function solveRoot(hDist: number, vDist: number, length: number): number
	local c = math.sqrt(length ^ 2 - vDist ^ 2)
	local function f(z)
		return hDist * math.sinh(z) / z - c
	end

	-- Bisection method to solve root of f(x)
	local maxIters = 32
	local lower = 0
	local upper = solveRootUpperBound(hDist, c)

	for i = 1, maxIters do
		local mid = (lower + upper) / 2
		local f_mid = f(mid)

		if math.abs(f_mid) < 1e-4 then
			return mid
		end

		if f_mid < 0 then
			lower = mid
		else
			upper = mid
		end
	end

	-- Failed to converge; return approximate solution
	return (lower + upper) / 2
end

function Catenary.Interpolant.new(p0: Point, p1: Point, length: number, upVector: Vector?)
	if length < (p1 - p0).Magnitude then
		-- Length of rope is too short; degenerate to a line
		return Linear.Interpolant.new(p0, p1)
	end

	local yAxis = if upVector then upVector.Unit else Vector3.yAxis
	local vDist = (p1 - p0):Dot(yAxis)
	local ground = (p1 - p0) - yAxis * vDist
	local xAxis = ground.Unit
	local hDist = ground.Magnitude
	local root = solveRoot(hDist, vDist, length)
	local scale = hDist / (2 * root)

	local hShift = (hDist - scale * math.log((length + vDist) / (length - vDist))) / 2
	local vShift = (vDist - length / math.tanh(root)) / 2

	local self = setmetatable(Position.Interpolant.new(), Catenary.Interpolant)

	self.p0 = p0
	self.ground = ground
	self.xAxis = xAxis
	self.yAxis = yAxis
	self.scale = scale
	self.hShift = hShift
	self.vShift = vShift

	self.Codimension = SplineUtils.GetCodimensionFromPoint(p0)
	self.Length = length

	return self
end

local function asinh(x)
	return math.log(x + math.sqrt(1 + x ^ 2))
end

function Catenary.Interpolant:SolvePosition(t: number): Point
	local p0 = self.p0
	local ground = self.ground
	local scale = self.scale
	local hShift = self.hShift
	local vShift = self.vShift
	local arcLength = self.Length * t

	local x = if self.IsUnitSpeed
		then scale * asinh(arcLength / scale + math.sinh(-hShift / scale)) + hShift
		else ground * t
	local y = scale * math.cosh((x - hShift) / scale) + vShift

	return p0 + self.xAxis * x + self.yAxis * y
end

function Catenary.Spline.new(points: { Point }, lengths: { number })
	assert(#lengths > 0)
	assert(#points == #lengths + 1)

	-- TODO: Get unique points
	-- TODO: Throw if not enough points left

	local interpolants = table.create(#lengths)

	local prevPoint = points[1]
	local totalLength = 0

	for i = 2, #points do
		local point = points[i]
		local interpolant = Catenary.Interpolant.new(prevPoint, point, lengths[i - 1])
		interpolants[i - 1] = interpolant
		totalLength += interpolant.Length
		prevPoint = point
	end

	local self = setmetatable(Position.Spline.new(), Catenary.Spline)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(prevPoint)
	self.Interpolants = interpolants
	self.Length = totalLength
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return Catenary
