--[[
	A generic spline class with spline-independent methods
	
	Notation
	t: Time parameter
	s: Arc length parameter
	r: Position
	T: Tangent vector
	N: Normal vector
	B: Binormal vector
	κ: Curvature
	τ: Torsion
]]

local Interpolation = require(script.Parent.Parent.Vendor.Interpolation)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local Position = {}

Position.Interpolant = {}
Position.Interpolant.__index = Position.Interpolant

Position.Spline = {}
Position.Spline.__index = Position.Spline

local PositionInterpolant = {}
PositionInterpolant.__index = PositionInterpolant

function Position.Interpolant.new()
	local self = setmetatable({}, Position.Interpolant)

	self._chebyshevInterpolant = nil
	self.IsUnitSpeed = false

	-- Global properties
	self.Codimension = nil
	self.Length = nil

	return self
end

function Position.Interpolant:_accountForUnitSpeed(t: number)
	return self.IsUnitSpeed and self._chebyshevInterpolant:evaluate(t * self.Length) or t
end

function Position.Interpolant:ToUnitSpeed()
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

function Position.Interpolant:SolveLength(from: number?, to: number?): number
	from = from or 0
	to = to or 1

	if from > to then
		from, to = to, from
	end

	return SplineUtils.GaussLegendre(function(t)
		return self:SolveVelocity(t).Magnitude
	end, from, to)
end

function Position.Interpolant:SolveVelocity(t: number): Vector
	error("Must implement SolveVelocity()")
end

function Position.Interpolant:SolveAcceleration(t: number): Vector
	error("Must implement SolveAcceleration()")
end

function Position.Interpolant:SolveJerk(t: number): Vector
	error("Must implement SolveJerk()")
end

function Position.Interpolant:SolveTangent(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- T(t) = r'(t) / |r'(t)|
	return self:SolveVelocity(t).Unit
end

function Position.Interpolant:SolveNormal(t: number): Vector
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

function Position.Interpolant:SolveBinormal(t: number): Vector
	if self.Codimension == 2 then
		t = self:_accountForUnitSpeed(t)

		-- T(t) x N(t)
		return self:SolveTangent(t):Cross(self:SolveNormal(t))
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function Position.Interpolant:SolveCurvature(t: number): number
	t = self:_accountForUnitSpeed(t)

	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local speed = vel.Magnitude
	local dTangent = accel / speed - vel * vel:Dot(accel) / speed ^ 3

	-- κ(t) = |T'(t)| / |r'(t)|
	return dTangent.Magnitude / speed
end

function Position.Interpolant:SolveTorsion(t: number): number
	t = self:_accountForUnitSpeed(t)

	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local jerk = self:SolveJerk(t)
	local cross = vel:Cross(accel)

	-- τ(t) = ((r'(t) x r''(t)) • r'''(t)) / |r'(t) x r''(t)|^2
	return cross:Dot(jerk) / cross.Magnitude ^ 2
end

function Position.Spline.new()
	local self = setmetatable({}, Position.Spline)

	self.Interpolants = nil
	self.InterpolantDomains = nil
	self.IsUnitSpeed = false

	-- Global properties
	self.Closed = nil -- TODO
	self.Codimension = nil
	self.Convex = nil -- TODO
	self.Length = nil
	self.Regular = nil -- TODO
	self.RotationIndex = nil -- TODO
	self.WindingNumber = nil -- TODO

	return self
end

function Position.Spline:ToUnitSpeed()
	self.IsUnitSpeed = true

	for _, interpolant in self.Interpolants do
		interpolant:ToUnitSpeed()
	end
end

function Position.Spline:SolveLength(a: number?, b: number?): number
	a = a or 0
	b = b or 1

	if a > b then
		a, b = b, a
	end

	local interp0Index, interp0Time = self:_getInterpolantFromPercentArcLength(a)
	local interp1Index, interp1Time = self:_getInterpolantFromPercentArcLength(b)

	local length0 = self.Interpolants[interp0Index]:SolveLength(interp0Time, 1)
	local length1 = self.Interpolants[interp1Index]:SolveLength(0, interp1Time)

	local intermediateLengths = 0
	for i = interp0Index + 1, interp1Index - 1 do
		intermediateLengths += self.Interpolants[i].Length
	end

	return length0 + intermediateLengths + length1
end

--[=[
	Binary search for the interpolant in the spline containing the point s% of
	the arc length along the spline.

	@return number -- Index of interpolant
	@return number -- Percent arc length along interpolant
--]=]
function Position.Spline:_getInterpolantFromPercentArcLength(s: number): (number, number)
	local domains = self.InterpolantDomains
	local numInterpolants = #domains

	-- There is only one option if there is one interpolant
	if numInterpolants == 1 then
		return 1, s
	end

	-- Special cases for when s is on the boundary or outside of [0, 1]
	if s < 0 then
		return 1, s
	elseif s == 0 then
		return 1, 0
	elseif s == 1 then
		return numInterpolants, 1
	elseif s > 1 then
		return numInterpolants, (s - domains[numInterpolants]) / (1 - domains[numInterpolants])
	end

	-- Binary search for the interpolant containing s
	local left = 1
	local right = numInterpolants

	while left <= right do
		local mid = math.floor((left + right) / 2)
		local intervalStart = domains[mid]

		if s >= intervalStart then
			local intervalEnd = mid == numInterpolants and 1 or domains[mid + 1]

			if s <= intervalEnd then
				local interpolantTime = (s - intervalStart) / (intervalEnd - intervalStart)

				return mid, interpolantTime
			else
				left = mid + 1
			end
		else
			right = mid - 1
		end
	end

	-- This is theoretically impossible
	error("Failed to get interpolant from s")
end

function Position.Spline:SolvePosition(t: number): Point
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolvePosition(interpTime)
end

function Position.Spline:SolveVelocity(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveVelocity(interpTime)
end

function Position.Spline:SolveAcceleration(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveAcceleration(interpTime)
end

function Position.Spline:SolveJerk(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveJerk(interpTime)
end

function Position.Spline:SolveTangent(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveTangent(interpTime)
end

function Position.Spline:SolveNormal(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveNormal(interpTime)
end

function Position.Spline:SolveBinormal(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveBinormal(interpTime)
end

function Position.Spline:SolveCurvature(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveCurvature(interpTime)
end

function Position.Spline:SolveTorsion(t: number): Vector
	local interpIndex, interpTime = self:_getInterpolantFromPercentArcLength(t)
	return self.Interpolants[interpIndex]:SolveTorsion(interpTime)
end

return Position
