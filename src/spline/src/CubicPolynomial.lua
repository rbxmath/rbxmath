local PositionSpline = require(script.Parent.PositionSpline)
local SplineUtils = require(script.Parent.SplineUtils)
local Types = require(script.Parent.Types)

local CubicPolynomial = setmetatable({}, PositionSpline)
CubicPolynomial.__index = CubicPolynomial

function CubicPolynomial.new()
	return PositionSpline.new()
end

-- This isn't correct yet
function CubicPolynomial:_cacheLengths()
	local totalLength = 0
	local lengths = table.create(self._numCurves)

	for i = 1, self._numCurves do
		local length = 1
		lengths[i] = length
		totalLength += length
	end

	self._curveDomains = SplineUtils.GetCurveDomains(lengths)
	self.Length = totalLength
end

-- This isn't correct yet
function CubicPolynomial:SolveLength()
	return 1
end

function CubicPolynomial:SolvePosition(t: number, curveIndex: number, curveTime: number): Point
	if not curveIndex or not curveTime then
		curveIndex, curveTime = self:_getCurveFromPercentArcLength(t)
	end
	local curve = self.Curves[curveIndex]

	-- r(t) in Horner form
	return curve.a + curveTime * (curve.b + curveTime * (curve.c + curveTime * curve.d))
end

function CubicPolynomial:SolveVelocity(t: number, curveIndex: number, curveTime: number): Vector
	if not curveIndex or not curveTime then
		curveIndex, curveTime = self:_getCurveFromPercentArcLength(t)
	end
	local curve = self.Curves[curveIndex]

	-- r'(t) in Horner form
	return curve.b + curveTime * (2 * curve.c + curveTime * 3 * curve.d)
end

function CubicPolynomial:SolveAcceleration(t: number, curveIndex: number, curveTime: number): Vector
	if not curveIndex or not curveTime then
		curveIndex, curveTime = self:_getCurveFromPercentArcLength(t)
	end
	local curve = self.Curves[curveIndex]

	-- r''(t)
	return 2 * curve.c + curveTime * 6 * curve.d
end

function CubicPolynomial:SolveJerk(t: number, curveIndex: number, curveTime: number): Vector
	if not curveIndex or not curveTime then
		curveIndex, curveTime = self:_getCurveFromPercentArcLength(t)
	end
	local curve = self.Curves[curveIndex]

	-- r'''(t)
	return 6 * curve.d
end

function CubicPolynomial:SolveMany(method: (t: number) -> any, numSamples: number, start: number, stop: number)
	start = start or 0
	stop = stop or 1

	assert(numSamples > 0)
	assert(start >= 0 and start <= 1)
	assert(stop >= 0 and stop <= 1)
	assert(start <= stop)

	local domains = self._curveDomains
	local curveIndex, curveTime = self:_getCurveFromPercentArcLength(start)
	local values = table.create(numSamples)

	for i = 1, numSamples do
		local t = (i - 1) / numSamples
		local s = SplineUtils.Lerp(start, stop, t)
		local nextDomain = curveIndex == self._numCurves and math.huge or domains[curveIndex + 1]

		if s >= nextDomain then
			curveIndex += 1
		end

		if s < 0 then
			curveTime = s
		elseif s == 0 then
			curveTime = 0
		elseif s == 1 then
			curveTime = 1
		elseif s > 1 then
			curveTime = (s - domains[self._numCurves]) / (1 - domains[self._numCurves])
		elseif curveIndex == self._numCurves then
			curveTime = (s - domains[self._numCurves]) / (1 - domains[self._numCurves])
		else
			curveTime = (s - domains[curveIndex]) / (domains[curveIndex + 1] - domains[curveIndex])
		end

		values[i] = method(self, nil, curveIndex, curveTime)
	end

	return values
end

return CubicPolynomial
