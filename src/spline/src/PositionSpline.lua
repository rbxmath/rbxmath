--[[
	A generic spline class with spline-independent methods
	
	Notation
	r(t): Position
	T(t): Tangent vector
	N(t): Normal vector
	B(t): Binormal vector
	K(t): Curvature
]]

local Types = require(script.Parent.Types)

local PositionSpline = {}
PositionSpline.__index = PositionSpline

function PositionSpline.new()
	return setmetatable({
		_unitSpeed = false,
		_curveDomains = {},
		_codimension = nil,
		_numCurves = nil,

		-- User-defined global properties
		Closed = nil,

		-- Computed global properties
		Regular = nil,
		Length = nil,
		Convex = nil,
		TotalCurvature = nil,

		-- Global properties for planar curves
		WindingNumber = nil,
		RotationIndex = nil,
	}, PositionSpline)
end

--[=[

	@return -- The same spline but with a unit-speed parametrization
--]=]
function PositionSpline:ToUnitSpeed() end

--[=[
	Binary search for the curve in the chain containing the point s% of the arc
	length along the chain.

	@return number -- Index of curve
	@return number -- Percent arc length along curve
--]=]
function PositionSpline:_getCurveFromPercentArcLength(s: number): (number, number)
	local domains = self._curveDomains
	local numCurves = #domains

	-- There is only one option if there is one curve
	if numCurves == 1 then
		return 1, s
	end

	-- Special cases for when s is on the boundary or outside of [0, 1]
	if s < 0 then
		return 1, s
	elseif s == 0 then
		return 1, 0
	elseif s == 1 then
		return numCurves, 1
	elseif s > 1 then
		return numCurves, (s - domains[numCurves]) / (1 - domains[numCurves])
	end

	-- Binary search for the spline containing s
	local left = 1
	local right = numCurves -- + 1

	while left <= right do
		local mid = math.floor((left + right) / 2)
		local intervalStart = domains[mid]

		if s >= intervalStart then
			local intervalEnd = mid == numCurves and 1 or domains[mid + 1]

			if s <= intervalEnd then
				local splineTime = (s - intervalStart) / (intervalEnd - intervalStart)
				return mid, splineTime
			else
				left = mid + 1
			end
		else
			right = mid - 1
		end
	end

	-- This is theoretically impossible
	error("Failed to get spline from s")
end

function PositionSpline:SolveTangent(t: number): Vector
	-- T(t) = r'(t) / |r'(t)|
	return self:SolveVelocity(t).Unit
end

function PositionSpline:SolveNormal(t: number): Vector
	if self._codimension == 0 then
		error("SolveNormal is restricted from splines in 1 dimension")
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

function PositionSpline:SolveBinormal(t: number): Vector
	if self._codimension == 2 then
		-- T(t) x N(t)
		return self:SolveTangent(t):Cross(self:SolveNormal(t))
	else
		error("SolveBinormal is restricted to splines in 3 dimensions")
	end
end

function PositionSpline:SolveCurvature(t: number): Vector
	local vel = self:SolveVelocity(t)
	local accel = self:SolveAcceleration(t)
	local speed = vel.Magnitude
	local dTangent = accel / speed - vel * vel:Dot(accel) / speed ^ 3

	-- K(t) = |T'(t)| / |r'(t)|
	-- N(t) is the direction of curvature
	local curvature = dTangent.Magnitude / speed
	local unitNormal = dTangent.Unit

	return curvature, unitNormal
end

return PositionSpline
