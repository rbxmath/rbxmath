local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local PositionSplineChain = {}

function PositionSplineChain.new()
	local self = setmetatable({}, PositionSplineChain)

	self.Codimension = nil
	self.IsUnitSpeed = false
	self.Length = nil
	self.Splines = nil
	self.SplineDomains = nil

	return self
end

function PositionSplineChain:ToUnitSpeed()
	self.IsUnitSpeed = true

	for _, spline in self.Splines do
		spline:ToUnitSpeed()
	end
end

function PositionSplineChain:SolveLength(a: number?, b: number?): number
	a = a or 0
	b = b or 1

	if a > b then
		a, b = b, a
	end

	local spline0Index, spline0Time = self:_getSplineFromPercentArcLength(a)
	local spline1Index, spline1Time = self:_getSplineFromPercentArcLength(b)

	local length0 = self.Splines[spline0Index]:SolveLength(spline0Time, 1)
	local length1 = self.Splines[spline1Index]:SolveLength(0, spline1Time)

	local intermediateLengths = 0
	for i = spline0Index + 1, spline1Index - 1 do
		intermediateLengths += self.Splines[i].Length
	end

	return length0 + intermediateLengths + length1
end

--[=[
	Binary search for the spline in the chain containing the point s% of the arc
	length along the chain.

	@return number -- Index of spline
	@return number -- Percent arc length along spline
--]=]
function PositionSplineChain:_getSplineFromPercentArcLength(s: number): (number, number)
	local domains = self.SplineDomains
	local numSplines = #domains

	-- There is only one option if there is one spline
	if numSplines == 1 then
		return 1, s
	end

	-- Special cases for when s is on the boundary or outside of [0, 1]
	if s < 0 then
		return 1, s
	elseif s == 0 then
		return 1, 0
	elseif s == 1 then
		return numSplines, 1
	elseif s > 1 then
		return numSplines, (s - domains[numSplines]) / (1 - domains[numSplines])
	end

	-- Binary search for the spline containing s
	local left = 1
	local right = numSplines -- + 1

	while left <= right do
		local mid = math.floor((left + right) / 2)
		local intervalStart = domains[mid]

		if s >= intervalStart then
			local intervalEnd = mid == numSplines and 1 or domains[mid + 1]

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

function PositionSplineChain:SolvePosition(self, t: number): Point
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolvePosition(splineTime)
end

function PositionSplineChain:SolveVelocity(self, t: number): Vector
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolveVelocity(splineTime)
end

function PositionSplineChain:SolveAcceleration(self, t: number): Vector
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolveAcceleration(splineTime)
end

function PositionSplineChain:SolveJerk(self, t: number): Vector
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolveJerk(splineTime)
end

function PositionSplineChain:SolveCurvature(self, t: number): Vector
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolveCurvature(splineTime)
end

function PositionSplineChain:SolveTorsion(self, t: number): Vector
	local splineIndex, splineTime = self:_getCurveFromPercentArcLength(t)
	return self.Splines[splineIndex]:SolveTorsion(splineTime)
end

return PositionSplineChain
