local Types = require(script.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local EPSILON = 1e-13

local SplineUtils = {}

--[=[
	Performs a fuzzy equality check.

	@param a --- A point
	@param b --- A point
	@return
--]=]
function SplineUtils.FuzzyEq(a: Point, b: Point): boolean
	local aType = typeof(a)

	if aType == "Vector2" then
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
		error("Bad point type")
	end
end

--[=[
	Linearly interpolates between two points.

	@param from -- From
	@param to -- To
	@param t -- Time in [0, 1]
--]=]
function SplineUtils.Lerp(from, to, t: number): number
	return from + t * (to - from)
end

--[=[
	Infers the ambient space of the spline from the dimension of the given point
	and returns the codimension of the spline (i.e., the dimension of the
	ambient space minus the dimension of the spline).

	@param point -- Sample point on spline
--]=]
function SplineUtils.GetCodimensionFromPoint(point: Point): number
	local pointType = typeof(point)

	if pointType == "Vector2" then
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
	Gets the percent arc length along the spline where each interpolant starts.

	@param interpolants -- A list of interpolants
	@return
--]=]
function SplineUtils.GetInterpolantDomains(interpolants): { number }
	local totalLength = 0

	for _, interpolant in interpolants do
		totalLength += interpolant.Length
	end

	local interpolantDomains = table.create(#interpolants)
	local runningLength = 0

	for i, interpolant in interpolants do
		interpolantDomains[i] = runningLength / totalLength
		runningLength += interpolant.Length
	end

	return interpolantDomains
end

--[=[
	Removes non-unique adjacent points using fuzzy equality.

	@param points -- A list of points
	@return 
--]=]
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

--[=[
	Integrates a function using 5-point Gaussian quadrature. Legendre roots and
	quadrature weights sourced from https://pomax.github.io/bezierinfo/legendre-gauss.html
	with highest precision afforded by 64-bit floats.

	@param f --- An integrable function
	@param from --- The lower integral bound
	@param to --- The upper integral bound
--]=]
function SplineUtils.GaussLegendre(f: (number) -> number, from: number, to: number): number
	if from == to then
		return 0
	end

	if from > to then
		from, to = to, from
	end

	local halfInterval = (to - from) / 2
	local midpoint = (to + from) / 2

	local x2 = halfInterval * 0.5384693101056831
	local x3 = halfInterval * 0.906179845938664

	return halfInterval
		* (
			0.5688888888888889 * f(midpoint)
			+ 0.47862867049936647 * (f(midpoint - x2) + f(midpoint + x2))
			+ 0.23692688505618908 * (f(midpoint - x3) + f(midpoint + x3))
		)
end

return SplineUtils
