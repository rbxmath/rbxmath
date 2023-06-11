local CatmullRom = require(script.Parent.Parent.Splines.CatmullRom)
local PositionSplineChain = require(script.Parent.PositionSplineChain)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point

local CatmullRomChain = setmetatable({}, PositionSplineChain)
CatmullRomChain.__index = CatmullRomChain

function CatmullRomChain.new(points: { Point }, alpha: number?, tension: number?)
	assert(type(points) == "table")
	assert(#points > 0)

	points = SplineUtils.GetUniquePoints(points)
	local numPoints = #points

	-- Early exits
	if numPoints == 1 then
		local self = setmetatable(PositionSplineChain.new(), CatmullRomChain)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(points[1])
		self.Length = 0
		self.Splines = { CatmullRom.fromPoint(points[1]) }
		self.SplineDomains = { 0 }

		return self
	elseif numPoints == 2 then
		local self = setmetatable(PositionSplineChain.new(), CatmullRomChain)

		self.Codimension = SplineUtils.GetCodimensionFromPoint(points[1])
		self.Length = (points[2] - points[1]).Magnitude
		self.Splines = { CatmullRom.fromLine(points[1], points[2]) }
		self.SplineDomains = { 0 }

		return self
	end

	-- Extrapolate to get 0th and n+1th points
	local veryFirstPoint = points[2]:Lerp(points[1], 2)
	local veryLastPoint = points[numPoints - 1]:Lerp(points[numPoints], 2)

	-- Create splines
	local numSplines = numPoints - 1
	local splines = table.create(numSplines)
	local totalLength = 0

	local window1 = veryFirstPoint
	local window2 = points[1]
	local window3 = points[2]
	local window4 = points[3]

	-- Add the first spline, whose first point is not in the points table
	splines[1] = CatmullRom.new(window1, window2, window3, window4, alpha, tension)
	totalLength += splines[1].Length

	-- Add the middle splines
	for i = 1, numPoints - 3 do
		-- Shift the window
		window1 = window2
		window2 = window3
		window3 = window4
		window4 = points[i + 3]

		local spline = CatmullRom.new(window1, window2, window3, window4, alpha, tension)
		totalLength += spline.Length
		splines[i + 1] = spline
	end

	-- Add the last spline, whose fourth point is not in the points table
	splines[numSplines] = CatmullRom.new(window2, window3, window4, veryLastPoint, alpha, tension)
	totalLength += splines[numSplines].Length

	local self = setmetatable(PositionSplineChain.new(), CatmullRomChain)

	self.Codimension = SplineUtils.GetCodimensionFromPoint(window1)
	self.Length = totalLength
	self.Splines = splines
	self.SplineDomains = SplineUtils.GetSplineDomains(splines)

	return CatmullRomChain
end

return CatmullRomChain
