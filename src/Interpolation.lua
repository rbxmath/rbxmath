--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require("../src/Tools")
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object
local Matrices = require("../src/Matrices")
local Matrix = require("../src/Matrices").Matrix
local FFT = require("../src/FastFourierTransform")

local Interpolation = {}

local _chebyshevGrid = function(n: number): Vector
	local result = {}

	for i = 0, n, 1 do
		result[n + 1 - i] = math.cos(i * math.pi / n)
	end

	return result
end

local function _barycentricInterpolationInChebyshevPointsAtAPoint(
	fList: Vector,
	x: number,
	chebyshevGridPoints: Vector
): number
	local n = #fList
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)

	local xDiffList = {}
	local diff
	for i, v in ipairs(chebyshevGrid) do
		diff = x - v
		if diff == 0 then
			return fList[i]
		end
		xDiffList[i] = (-1) ^ (i - 1) / diff
	end

	local val = 0.5 * xDiffList[1]
	local denominator = val
	local numerator = fList[1] * val
	for i = 2, n - 1, 1 do
		val = xDiffList[i]
		numerator = numerator + fList[i] * val
		denominator = denominator + val
	end
	val = 0.5 * xDiffList[n]
	numerator = numerator + fList[n] * val
	denominator = denominator + val

	return numerator / denominator
end

local function _barycentricInterpolationInChebyshevPointsAtPointList(
	fList: Vector,
	xList: Vector,
	chebyshevGridPoints: Vector
): Vector
	local n = #fList
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
	local result = {}
	for j = 1, #xList, 1 do
		local x = xList[j]
		local xDiffList = {}
		local diff
		local index = -1
		for i, v in ipairs(chebyshevGrid) do
			diff = x - v
			if diff == 0 then
				index = i
			end
			xDiffList[i] = (-1) ^ (i - 1) / diff
		end

		if index ~= -1 then
			result[j] = fList[index]
		else
			local val = 0.5 * xDiffList[1]
			local denominator = val
			local numerator = fList[1] * val
			for i = 2, n - 1, 1 do
				val = xDiffList[i]
				numerator = numerator + fList[i] * val
				denominator = denominator + val
			end
			val = 0.5 * xDiffList[n]
			numerator = numerator + fList[n] * val
			denominator = denominator + val

			result[j] = numerator / denominator
		end
	end

	return result
end

local function _barycentricInterpolationSolve(
	fList: Vector,
	t: number,
	tol: number,
	gapTol: number,
	chebyshevGridPoints: Vector
): number
	local n = #fList
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n)
	local tolerance = tol or (10 ^ -13)
	local gapTolerance = gapTol or (10 ^ -13)

	local intervalList = {}
	local rootList = {}
	if math.abs(fList[1] - t) < tolerance then
		rootList[#rootList + 1] = chebyshevGrid[1]
	end
	for i = 1, n - 1, 1 do
		if
			t <= fList[i + 1] - tolerance and t >= fList[i] + tolerance
			or t <= fList[i] - tolerance and t >= fList[i + 1] + tolerance
		then
			intervalList[#intervalList + 1] = { chebyshevGrid[i], chebyshevGrid[i + 1] }
		end
		if math.abs(fList[i + 1] - t) < tolerance then
			rootList[#rootList + 1] = chebyshevGrid[i + 1]
		end
	end

	while #intervalList ~= 0 do
		local newIntervalList = {}
		local x0, x1, x2, f0, f1, f2, temp
		for i, v in ipairs(intervalList) do
			x0, x1 = v[1], v[2]
			temp = _barycentricInterpolationInChebyshevPointsAtPointList(fList, { x0, x1 }, chebyshevGrid)
			f0, f1 = temp[1] - t, temp[2] - t
			x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
			if x2 < x0 or x2 > x1 then
				x2 = (x0 + x1) / 2
			end
			f2 = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x2, chebyshevGrid) - t
			if math.abs(f2) < tolerance then
				rootList[#rootList + 1] = x2
			else
				local gap = gapTolerance + (x1 - x0) + 1
				while math.abs(f2) >= tolerance and gap >= gapTolerance and gap ~= x1 - x0 do
					gap = x1 - x0
					if
						((f0 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f0))
						and ((f1 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f1))
					then
						newIntervalList[#newIntervalList + 1] = { x2, x1 }
						x1 = x2
						f1 = f2
					elseif (f0 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f0) then
						x1 = x2
						f1 = f2
					elseif (f1 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f1) then
						x0 = x2
						f0 = f2
					else
						break
					end
					x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
					if x2 < x0 or x2 > x1 then
						x2 = (x0 + x1) / 2
					end
					f2 = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x2, chebyshevGrid) - t
				end
				if math.abs(f2) < tolerance or gap < gapTolerance then
					rootList[#rootList + 1] = x2
				end
			end
		end
		intervalList = newIntervalList
	end
	return rootList
end

local function _chebyshevSpectralDifferentionMatrix(n: number, chebyshevGridPoints: Vector): Tensor
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n)
	local d = { {}, {} }

	for i = 1, n + 1, 1 do
		d[1][i] = 1
		d[2][i] = 0
	end

	for i = 1, n + 1, 1 do
		for j = 1, n + 1, 1 do
			if i ~= j then
				d[2][i] = d[2][i] + d[1][i] / (chebyshevGrid[i] - chebyshevGrid[j])
			end
		end
	end

	local D = {}

	for i = 1, n + 1, 1 do
		D[i] = {}
	end

	local c = {}
	for i = 1, n + 1, 1 do
		c[i] = 1
		for j = 1, n + 1, 1 do
			if i ~= j then
				c[i] = c[i] * (chebyshevGrid[i] - chebyshevGrid[j])
			end
		end
	end

	for i = 1, n + 1, 1 do
		for j = 1, n + 1, 1 do
			if i == j then
				D[i][j] = d[2][i]
			else
				D[i][j] = c[i] / (c[j] * (chebyshevGrid[i] - chebyshevGrid[j]))
			end
		end
	end

	return Matrix:new(D)
end

local function _chebyshevHigherOrderSpectralDifferentionMatrix(
	n: number,
	p: number,
	chebyshevGridPoints: Vector
): Tensor
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n)
	local d = {}

	for i = 1, p + 1, 1 do
		d[i] = {}
	end

	for i = 1, n + 1, 1 do
		d[1][i] = 1
		for j = 2, p + 1, 1 do
			d[j][i] = 0
		end
	end

	for i = 1, n + 1, 1 do
		for j = 1, n + 1, 1 do
			if i ~= j then
				for k = p + 1, 2, -1 do
					d[k][i] = d[k][i] + (k - 1) * d[k - 1][i] / (chebyshevGrid[i] - chebyshevGrid[j])
				end
			end
		end
	end

	local D = {}

	for i = 1, p + 1, 1 do
		D[i] = {}
		for j = 1, n + 1, 1 do
			D[i][j] = {}
			for k = 1, n + 1, 1 do
				if j == k then
					D[i][j][k] = d[i][j]
				elseif i == 1 then
					D[i][j][k] = 0
				end
			end
		end
	end

	local c = {}
	for i = 1, n + 1, 1 do
		c[i] = 1
		for j = 1, n + 1, 1 do
			if i ~= j then
				c[i] = c[i] * (chebyshevGrid[i] - chebyshevGrid[j])
			end
		end
	end

	for i = 1, n + 1, 1 do
		for j = 1, n + 1, 1 do
			for k = 2, p + 1, 1 do
				if i ~= j then
					D[k][i][j] = (k - 1) * (d[k - 1][i] - D[k - 1][i][j]) / (chebyshevGrid[i] - chebyshevGrid[j])
				end
			end
		end
	end

	for i = 1, n + 1, 1 do
		for j = 1, n + 1, 1 do
			for k = 1, p + 1, 1 do
				if i ~= j then
					D[k][i][j] = D[k][i][j] * c[i] / c[j]
				end
			end
		end
	end

	local result = {}
	for i = 1, p + 1, 1 do
		result[i] = Matrix:new(D[i])
	end

	return result
end

local function _linearRescalingFunction(a: number, b: number): ScalarMap
	return function(x)
		return (b - a) / 2 * (x - 1) + b
	end
end

local function _inverseLinearRescalingFunction(a: number, b: number): ScalarMap
	return function(x)
		return 2 * (x - a) / (b - a) - 1
	end
end

local ChebyshevInterpolant = {
	grid = {},
	gridValues = {},
	degree = 0,
	linearRescalingFunction = nil,
	inverseLinearRescalingFunction = nil,
	leftBound = 0,
	rightBound = 0,
	solveMethod = "RegulaFalsi",
	evaluationFunction = nil,
}

-- Creates a Chebyshev interpolant, p(x), of a function/list of values f
function ChebyshevInterpolant:new(
	f: Vector | ScalarMap,
	a: number,
	b: number,
	n: number,
	chebyshevGrid: Vector,
	solveMethod: Vector
): Object
	solveMethod = solveMethod or "RegulaFalsi"
	local result = {}
	setmetatable(result, self)
	self.__index = self
	if type(f) == "table" then
		n = #f - 1
	end
	chebyshevGrid = chebyshevGrid or _chebyshevGrid(n)
	local fList = {}
	local rescalingFunction = _linearRescalingFunction(a, b)
	if type(f) == "function" then
		for i = 1, n + 1, 1 do
			fList[i] = f(rescalingFunction(chebyshevGrid[i]))
		end
	else
		fList = f
	end
	result.grid = chebyshevGrid
	result.gridValues = fList
	result.degree = n
	result.linearRescalingFunction = rescalingFunction
	result.inverseLinearRescalingFunction = _inverseLinearRescalingFunction(a, b)
	result.leftBound = a
	result.rightBound = b
	result.solveMethod = solveMethod
	result.evaluationFunction = function(x)
		return result:evaluate(x)
	end
	result.coefficientList = Tools.list.scale(FFT:FCT1(Tools.list.reverse(fList)), 2)

	return result
end

function ChebyshevInterpolant:adaptive(f: ScalarMap, a: number, b: number, tol: number): Object
	tol = tol or 10 ^ -10
	local n = 8
	local currentGrid = {}
	local currentGridValues = {}
	local norm = tol + 1
	local currentInterpolant = {}
	local rescalingFunction = _linearRescalingFunction(a, b)
	while norm >= tol do
		n *= 2
		currentGrid = _chebyshevGrid(n)
		for i = n + 1, 1, -1 do
			if i % 2 == 1 then
				currentGridValues[i] = currentGridValues[(i - 1) / 2 + 1] or f(rescalingFunction(currentGrid[i]))
			else
				currentGridValues[i] = f(rescalingFunction(currentGrid[i]))
			end
		end
		currentInterpolant = FFT:FCT1(Tools.list.reverse(currentGridValues))
		norm = 0
		for i = n / 2 + 1, n + 1, 1 do
			norm = math.max(norm, math.abs(currentInterpolant[i]))
		end
	end
	local index = n + 1
	while index >= 1 do
		if math.abs(currentInterpolant[index]) > tol then
			break
		else
			index -= 1
		end
	end
	return self:new(f, a, b, index - 1)
end

-- Computes p(x)
function ChebyshevInterpolant:evaluate(x: number): number
	local fList = self.gridValues
	local n = self.degree + 1
	local chebyshevGrid = self.grid
	x = self.inverseLinearRescalingFunction(x)

	local xDiffList = {}
	local diff
	for i, v in ipairs(chebyshevGrid) do
		diff = x - v
		if diff == 0 then
			return fList[i]
		end
		xDiffList[i] = (-1) ^ (i - 1) / diff
	end

	local val = 0.5 * xDiffList[1]
	local denominator = val
	local numerator = fList[1] * val
	for i = 2, n - 1, 1 do
		val = xDiffList[i]
		numerator = numerator + fList[i] * val
		denominator = denominator + val
	end
	val = 0.5 * xDiffList[n]
	numerator = numerator + fList[n] * val
	denominator = denominator + val

	return numerator / denominator
end

-- Solves the problem p(x) = t for x given t
function ChebyshevInterpolant:solve(t: number, tol: number, method: string, gapTol: number): number
	method = method or self.solveMethod

	if method == "Monotone" then
		return Tools.list.map(
			self.linearRescalingFunction,
			_barycentricInterpolationSolve(self.gridValues, t, tol, gapTol, self.grid)
		)
	elseif method == "RegulaFalsi" then
		tol = tol or 10 ^ -13
		local boundList = { self.leftBound, self.rightBound }
		local n = self.degree
		local derivativeList = self:derivativeList(n - 2)
		for i = n - 2, 1, -1 do
			local newBoundList = { self.leftBound }
			local derivative = derivativeList[i]
			if math.max(table.unpack(Tools.list.map(math.abs, derivative.gridValues))) >= tol then
				for j = 1, #boundList - 1, 1 do
					newBoundList[#newBoundList + 1] = Tools.solve.regulaFalsi(function(x)
						return derivative:evaluate(x)
					end, 0, boundList[j], boundList[j + 1], tol)
				end
				newBoundList[#newBoundList + 1] = self.rightBound
				boundList = newBoundList
			end
		end
		local pointList = {}
		for i = 1, #boundList - 1, 1 do
			pointList[#pointList + 1] = Tools.solve.regulaFalsi(function(x)
				return self:evaluate(x)
			end, t, boundList[i], boundList[i + 1], tol)
		end
		return pointList
	end
end

function ChebyshevInterpolant:derivative(n: number): Object
	n = n or 1
	local derivativeMatrix = _chebyshevHigherOrderSpectralDifferentionMatrix(self.degree, 1, self.grid)[2]
	local fVector = {}
	for key, value in ipairs(self.gridValues) do
		fVector[key] = { value }
	end
	local fColumn = Matrix:new(fVector):scaled(2 * (self.rightBound - self.leftBound))
	for i = 1, n, 1 do
		fColumn = derivativeMatrix * fColumn
	end
	local fList = fColumn:transpose()[1]
	return ChebyshevInterpolant:new(fList, self.leftBound, self.rightBound, self.degree, self.grid)
end

function ChebyshevInterpolant:derivativeList(n: number): { [number]: Object }
	n = n or 1
	local derivativeMatrix = _chebyshevHigherOrderSpectralDifferentionMatrix(self.degree, 1, self.grid)[2]
	local derivativeList = {}
	local fVector = {}
	for key, value in ipairs(self.gridValues) do
		fVector[key] = { value }
	end
	local fColumn = Matrix:new(fVector):scaled(2 * (self.rightBound - self.leftBound))
	for i = 1, n, 1 do
		fColumn = derivativeMatrix * fColumn
		derivativeList[i] =
			ChebyshevInterpolant:new(fColumn:transpose()[1], self.leftBound, self.rightBound, self.degree, self.grid)
	end

	return derivativeList
end

function ChebyshevInterpolant:max(method: string): number
	method = method or self.solveMethod
	local derivative = self:derivative()
	local extrema = derivative:solve(0, nil, method)
	extrema[#extrema + 1] = self.leftBound
	extrema[#extrema + 1] = self.rightBound
	local fList = {}
	for i = 1, #extrema, 1 do
		fList[i] = self:evaluate(extrema[i])
	end
	return math.max(table.unpack(fList))
end

function ChebyshevInterpolant:min(method: string): number
	method = method or self.solveMethod
	local derivative = self:derivative()
	local extrema = derivative:solve(0, nil, method)
	extrema[#extrema + 1] = self.leftBound
	extrema[#extrema + 1] = self.rightBound
	local fList = {}
	for i = 1, #extrema, 1 do
		fList[i] = self:evaluate(extrema[i])
	end
	return math.min(table.unpack(fList))
end

function ChebyshevInterpolant:extremes(method: string): Vector
	method = method or self.solveMethod
	local derivative = self:derivative()
	local extrema = derivative:solve(0, nil, method)
	extrema[#extrema + 1] = self.leftBound + 10 ^ -16
	extrema[#extrema + 1] = self.rightBound - 10 ^ -16
	local fList = {}
	for i = 1, #extrema, 1 do
		fList[i] = self:evaluate(extrema[i])
	end
	return math.max(table.unpack(fList)), math.min(table.unpack(fList))
end

function ChebyshevInterpolant:inverse(method: string): Object
	method = method or self.solveMethod
	local max, min = self:extremes(method)
	local fList = {}
	local linearRescalingFunction = _linearRescalingFunction(min, max)
	for key, value in ipairs(self.grid) do
		local save = self:solve(linearRescalingFunction(value), nil, method)
		if #save > 1 or #save == 0 then
			error("Function is not invertible!", -1)
		else
			fList[key] = save[1]
		end
	end
	return ChebyshevInterpolant:new(fList, min, max, self.degree, self.grid)
end

Interpolation.Chebyshev = {}

Interpolation.Chebyshev.grid = _chebyshevGrid

function Interpolation.Chebyshev.evaluate(f: ScalarMap, x: number, n: number, chebyshevGrid: Vector): number
	chebyshevGrid = chebyshevGrid or _chebyshevGrid(n)
	local fList = {}
	for i = 1, n + 1, 1 do
		fList[i] = f(chebyshevGrid[i])
	end
	return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x, chebyshevGrid)
end

function Interpolation.Chebyshev.rescaleAndEvaluate(
	f: ScalarMap,
	a: number,
	b: number,
	x: number,
	n: number,
	chebyshevGrid: Vector
): number
	chebyshevGrid = chebyshevGrid or _chebyshevGrid(n - 1)
	local rescalingFunction = _linearRescalingFunction(a, b)
	local fList = {}
	for i = 1, n, 1 do
		fList[i] = f(rescalingFunction(chebyshevGrid[i]))
	end
	return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x)
end

function Interpolation.Chebyshev.evaluateOnData(fList: Vector, x: number, chebyshevGrid: Vector): Vector
	return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x, chebyshevGrid)
end

function Interpolation.Chebyshev.evaluateAtPointList(
	f: ScalarMap,
	xList: Vector,
	n: number,
	chebyshevGridPoints: Vector
): Vector
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
	local fList = {}
	for i = 1, n, 1 do
		fList[i] = f(chebyshevGrid[i])
	end

	return _barycentricInterpolationInChebyshevPointsAtPointList(fList, xList, chebyshevGrid)
end

function Interpolation.Chebyshev.evaluateAtPointListOnData(fList: Vector, xList: Vector, chebyshevGrid: Vector): Vector
	chebyshevGrid = chebyshevGrid or _chebyshevGrid(#fList - 1)
	return _barycentricInterpolationInChebyshevPointsAtPointList(fList, xList, chebyshevGrid)
end

function Interpolation.Chebyshev.benchmark(
	f: ScalarMap,
	xList: Vector,
	deltaMax: number,
	maxIters: number
): (Vector, Vector)
	local delta = 0
	local iters = 0
	local mark, realAnswer, approximateAnswer
	local errorVector, timeVector = {}, {}

	while delta < deltaMax and iters < maxIters do
		iters = iters + 1
		realAnswer = Tools.list.map(f, xList)
		mark = os.clock()
		approximateAnswer = Interpolation.Chebyshev.evaluateAtPointList(f, xList, iters + 1)
		delta = os.clock() - mark
		errorVector[iters] = Tools.list.error(realAnswer, approximateAnswer)
		timeVector[iters] = delta
	end

	return errorVector, timeVector
end

function Interpolation.Chebyshev.solve(
	f: ScalarMap,
	t: number,
	n: number,
	method: string,
	tol: number,
	chebyshevGridPoints: Vector
): number
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
	local tolerance = tol or (10 ^ -13)
	local fList = {}
	for i = 1, n, 1 do
		fList[i] = f(chebyshevGrid[i])
	end

	local methodToUse = method or "RegulaFalsi"

	if methodToUse == "RegulaFalsi" then
		return _barycentricInterpolationSolve(fList, t, tolerance, nil, chebyshevGrid)
	end
end

function Interpolation.Chebyshev.solveFromData(
	fList: Vector,
	t: number,
	method: string,
	tol: number,
	chebyshevGridPoints: Vector
): number
	local n = #fList
	local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
	local tolerance = tol or (10 ^ -13)

	local methodToUse = method or "RegulaFalsi"

	if methodToUse == "RegulaFalsi" then
		return _barycentricInterpolationSolve(fList, t, tolerance, nil, chebyshevGrid)
	end
end

Interpolation.Chebyshev.derivativeMatrix = _chebyshevSpectralDifferentionMatrix

Interpolation.Chebyshev.linearRescalingFunction = _linearRescalingFunction

Interpolation.ChebyshevInterpolant = ChebyshevInterpolant

return Interpolation
