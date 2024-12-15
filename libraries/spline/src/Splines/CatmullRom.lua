local CubicPolynomial = require(script.Parent.CubicPolynomial)
local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local DEFAULT_ALPHA = 0.5
local DEFAULT_TENSION = 0

local CatmullRom = {}

CatmullRom.Interpolant = setmetatable({}, CubicPolynomial.Interpolant)
CatmullRom.Interpolant.__index = CatmullRom.Interpolant

CatmullRom.Spline = setmetatable({}, Position.Spline)
CatmullRom.Spline.__index = CatmullRom.Spline

function CatmullRom.Interpolant.new(p0: Point, p1: Point, p2: Point, p3: Point, alpha: number?, tension: number?)
	alpha = alpha or DEFAULT_ALPHA
	tension = tension or DEFAULT_TENSION
	assert(alpha)
	assert(tension)

	-- https://qroph.github.io/2018/07/30/smooth-paths-using-catmull-rom-splines.html
	local p1_p0 = p1 - p0
	local p2_p1 = p2 - p1
	local p3_p2 = p3 - p2

	local t0 = 0
	local t1 = vector.magnitude(p1_p0) ^ alpha + t0
	local t2 = vector.magnitude(p2_p1) ^ alpha + t1
	local t3 = vector.magnitude(p3_p2) ^ alpha + t2

	local scalar = (1 - tension) * (t2 - t1)
	local m1 = scalar * (p1_p0 / (t1 - t0) - (p2 - p0) / (t2 - t0) + p2_p1 / (t2 - t1))
	local m2 = scalar * (p2_p1 / (t2 - t1) - (p3 - p1) / (t3 - t1) + p3_p2 / (t3 - t2))

	local a = p1
	local b = m1
	local c = 3 * p2_p1 - 2 * m1 - m2
	local d = -2 * p2_p1 + m1 + m2

	return setmetatable(CubicPolynomial.Interpolant.new(a, b, c, d), CatmullRom.Interpolant)
end

function CatmullRom.Interpolant.fromPoint(p0: Point)
	local zero = p0 * 0

	return setmetatable(CubicPolynomial.Interpolant.new(zero, zero, zero, p0), CatmullRom.Interpolant)
end

function CatmullRom.Interpolant.fromLine(p0: Point, p1: Point)
	local zero = p0 * 0
	local p1_p0 = p1 - p0

	return setmetatable(CubicPolynomial.Interpolant.new(zero, zero, p1_p0, p0), CatmullRom.Interpolant)
end

function CatmullRom.Spline.new(points: { Point }, alpha: number?, tension: number?)
	assert(type(points) == "table")
	assert(#points > 0)

	points = SplineUtils.GetUniquePoints(points)
	local numPoints = #points

	-- Early exits
	if numPoints == 1 then
		local self = setmetatable(Position.Spline.new(), CatmullRom.Spline)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(points[1])
		self.Length = 0
		self.Interpolants = { CatmullRom.Interpolant.fromPoint(points[1]) }
		self.InterpolantDomains = { 0 }

		return self
	elseif numPoints == 2 then
		local self = setmetatable(Position.Spline.new(), CatmullRom.Spline)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(points[1])
		self.Length = vector.magnitude(points[2] - points[1])
		self.Interpolants = { CatmullRom.Interpolant.fromLine(points[1], points[2]) }
		self.InterpolantDomains = { 0 }

		return self
	end

	-- Extrapolate to get 0th and n+1th points
	local veryFirstPoint = points[2]:Lerp(points[1], 2)
	local veryLastPoint = points[numPoints - 1]:Lerp(points[numPoints], 2)

	-- Create interpolants
	local numInterpolants = numPoints - 1
	local interpolants = table.create(numInterpolants)
	local totalLength = 0

	local window1 = veryFirstPoint
	local window2 = points[1]
	local window3 = points[2]
	local window4 = points[3]

	-- Add the first interpolant, whose first point is not in the points table
	interpolants[1] = CatmullRom.Interpolant.new(window1, window2, window3, window4, alpha, tension)
	totalLength += interpolants[1].Length

	-- Add the middle interpolant
	for i = 1, numPoints - 3 do
		-- Shift the window
		window1 = window2
		window2 = window3
		window3 = window4
		window4 = points[i + 3]

		local interpolant = CatmullRom.Interpolant.new(window1, window2, window3, window4, alpha, tension)
		totalLength += interpolant.Length
		interpolants[i + 1] = interpolant
	end

	-- Add the last interpolant, whose fourth point is not in the points table
	interpolants[numInterpolants] = CatmullRom.Interpolant.new(window2, window3, window4, veryLastPoint, alpha, tension)
	totalLength += interpolants[numInterpolants].Length

	local self = setmetatable(Position.Spline.new(), CatmullRom.Spline)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(window1)
	self.Length = totalLength
	self.Interpolants = interpolants
	self.InterpolantDomains = SplineUtils.GetInterpolantDomains(interpolants)

	return self
end

return CatmullRom
