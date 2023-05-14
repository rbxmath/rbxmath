local Types = require(script.Parent.Types)

local SplineUtils = {}

function SplineUtils.FuzzyEq(a: Point, b: Point): boolean
	local aType = typeof(a)

	if aType == "number" then
		return a == b or math.abs(a - b) <= (math.abs(a) + 1) * EPSILON
	elseif aType == "Vector2" then
		local aX, bX = a.X, b.X
		local aY, bY = a.Y, b.Y
		return aX == bX and aY == bY
			or math.abs(aX - bX) <= (math.abs(aX) + 1) * EPSILON
				and math.abs(aY - bY) <= (math.abs(aY) + 1) * EPSILON
	elseif aType == "Vector3" then
		return a:FuzzyEq(b, EPSILON)
	elseif aType == "table" then
		for i, coord in a do
			if math.abs(coord - b[i]) > (math.abs(coord) + 1) * EPSILON then
				return false
			end
		end
		return true
	else
		error("Bad a")
	end
end

function SplineUtils.Lerp(a: number, b: number, t): number
	return a + t * (b - a)
end

function SplineUtils.GetCodimensionFromPoint(point: Point): number
	local pointType = typeof(point)

	if pointType == "number" then
		return 0
	elseif pointType == "Vector2" then
		return 1
	elseif pointType == "Vector3" then
		return 2
	elseif pointType == "table" then
		return #point - 1
	else
		error("Bad point")
	end
end

--[=[
	Gets the percent arc length along the chain where each curve starts

	@param lengths -- Arc lengths of curves
--]=]
function SplineUtils.GetCurveDomains(lengths: { number }): { number }
	local totalLength = 0
	for _, length in lengths do
		totalLength += length
	end

	local curveDomains = table.create(#lengths)
	local runningLength = 0

	for i, length in lengths do
		curveDomains[i] = runningLength / totalLength
		runningLength += length
	end

	return curveDomains
end

function SplineUtils.GetUniquePoints(points: { Point }): { Point }
	local prevPoint = points[1]
	local uniquePoints = { prevPoint }
	local i = 2

	for j = 2, #points do
		local point = points[j]

		if not SplineUtils.FuzzyEq(point, prevPoint) then
			uniquePoints[i] = point
			prevPoint = point
			i += 1
		end
	end

	return uniquePoints
end

return SplineUtils
