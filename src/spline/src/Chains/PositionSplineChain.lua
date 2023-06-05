local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)
local Vector = require(script.Parent.Parent.Vector)

type Point = Types.Point
type Vector = Types.Vector

local PositionSplineChain = {}

function PositionSplineChain.new(curves)
	local self = {}

	self._curveDomains = {}
end

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

function PositionSplineChain.SolvePosition(self, t: number): Vector end

return PositionSplineChain
