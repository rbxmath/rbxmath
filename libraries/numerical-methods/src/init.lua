--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require("../src/Tools")
local MA = require("../src/Matrices")

local NumericalMethods = {}

NumericalMethods.solvers = {}

function NumericalMethods.solvers.regulaFalsi(
	f: (number) -> number,
	t: number,
	a: number,
	b: number,
	tol: number
): number
	local tolerance = tol or 10 ^ -13
	local leftValue, rightValue, middleValue = f(a) - t, f(b) - t, 0
	local left, right, middle = a, b, (a + b) / 2
	if leftValue * rightValue > 0 then
		return nil
	elseif math.abs(leftValue) < tolerance then
		return left
	elseif math.abs(rightValue) < tolerance then
		return right
	end

	while math.abs(leftValue - rightValue) >= tolerance and math.abs(left - right) >= tolerance do
		middle = (right * leftValue - left * rightValue) / (leftValue - rightValue)
		middleValue = f(middle) - t
		if
			math.abs(middleValue) < tolerance
			or math.abs(left - middle) / math.abs(leftValue - rightValue) < tolerance
			or math.abs(right - middle) / math.abs(leftValue - rightValue) < tolerance
		then
			return middle
		elseif math.abs(leftValue) < tolerance then
			return left
		elseif math.abs(rightValue) < tolerance then
			return right
		elseif leftValue * middleValue > 0 then
			leftValue = middleValue
			left = middle
		else
			rightValue = middleValue
			right = middle
		end
	end

	return middle
end

function NumericalMethods.solvers.newtonsMethod(
	f: (number) -> number,
	fPrime: (number) -> number,
	t: number,
	x: number,
	tol: number
): number
	tol = tol or 10 ^ -13
	x = x or 0
	local value = f(x) - t
	local valuePrime = fPrime(x)
	while math.abs(value) >= tol do
		x = x - value / valuePrime
		value = f(x) - t
		valuePrime = fPrime(x)
	end
	return x
end

NumericalMethods.integration = {}

NumericalMethods.integration.fivePointGaussianGrid = {
	-1 * math.sqrt(5 + 2 * math.sqrt(10 / 7)) / 3,
	-1 * math.sqrt(5 - 2 * math.sqrt(10 / 7)) / 3,
	0,
	math.sqrt(5 - 2 * math.sqrt(10 / 7)) / 3,
	math.sqrt(5 + 2 * math.sqrt(10 / 7)) / 3,
}
NumericalMethods.integration.fivePointGaussianWeights = {
	(322 - 13 * math.sqrt(70)) / 900,
	(322 + 13 * math.sqrt(70)) / 900,
	128 / 225,
	(322 + 13 * math.sqrt(70)) / 900,
	(322 - 13 * math.sqrt(70)) / 900,
}

function NumericalMethods.integration.fivePointGaussianQuadrature(f: (number) -> number, a: number, b: number): number
	local sum = 0
	for i = 1, 5, 1 do
		local x = NumericalMethods.integration.fivePointGaussianGrid[i]
		local w = NumericalMethods.integration.fivePointGaussianWeights[i]
		local trueX = ((b - a) * x + b + a) / 2
		sum = sum + w * f(trueX)
	end
	return (b - a) * sum / 2
end

function NumericalMethods.integration.adaptiveQuadrature(
	f: (number) -> number,
	a: number,
	b: number,
	tolerance: number,
	guess: number
): number
	tolerance = tolerance or 10 ^ -13
	guess = guess or NumericalMethods.integration.fivePointGaussianQuadrature(f, a, b)
	if math.abs(b - a) < tolerance then
		return guess
	end
	local midPoint = (b + a) / 2
	local left = NumericalMethods.integration.fivePointGaussianQuadrature(f, a, midPoint)
	local right = NumericalMethods.integration.fivePointGaussianQuadrature(f, midPoint, b)
	local splitGuess = left + right
	if math.abs(guess - splitGuess) < tolerance then
		return splitGuess
	else
		return NumericalMethods.integration.adaptiveQuadrature(f, a, midPoint, tolerance, left)
			+ NumericalMethods.integration.adaptiveQuadrature(f, midPoint, b, tolerance, right)
	end
end

NumericalMethods.integration.integrate = NumericalMethods.integration.adaptiveQuadrature

NumericalMethods.integration.fivePointLaguerreGrid =
	{ 0.2635603197181409, 1.4134030591065208, 3.5964257710407237, 7.085810005858842, 12.640800844275782 }
NumericalMethods.integration.fivePointLaguerreWeights =
	{ 0.5217556105828091, 0.3986668110831633, 0.07594244968170724, 0.003611758679922023, 0.000023369972385776235 }

function NumericalMethods.integration.fivePointLaguerre(f: (number) -> number, a: number, negative: boolean): number
	local sum = 0
	for i = 1, 5, 1 do
		local x = NumericalMethods.integration.fivePointLaguerreGrid[i]
		local w = NumericalMethods.integration.fivePointLaguerreWeights[i]
		local trueX
		if negative then
			trueX = -a - x
		else
			trueX = x + a
		end
		sum = sum + w * f(trueX) * math.exp(x)
	end
	if negative then
		return -sum
	else
		return sum
	end
end

function NumericalMethods.integration.adaptiveLaguerre(
	f: (number) -> number,
	a: number,
	negative: boolean,
	tolerance: number,
	guess: number,
	stride: number,
	strideTolerance: number
): number
	tolerance = tolerance or 10 ^ -13
	stride = stride or 1
	guess = guess or NumericalMethods.integration.fivePointLaguerre(f, a, negative)
	strideTolerance = strideTolerance or 2 ^ 16
	if stride > strideTolerance then
		error(
			"Adaptive Laguerre Quadrature failed to converge! Either try a higher strideTolerance or the integral may not converge."
		)
	end
	if negative then
		local left = NumericalMethods.integration.adaptiveQuadrature(f, a - stride, a, tolerance)
		local right = NumericalMethods.integration.fivePointLaguerre(f, a - stride, negative)
		local splitGuess = left + right
		if math.abs(guess - splitGuess) < tolerance then
			return splitGuess
		else
			return left
				+ NumericalMethods.integration.adaptiveLaguerre(
					f,
					a - stride,
					negative,
					tolerance,
					right,
					2 * stride,
					strideTolerance
				)
		end
	else
		local left = NumericalMethods.integration.adaptiveQuadrature(f, a, a + stride, tolerance)
		local right = NumericalMethods.integration.fivePointLaguerre(f, a + stride)
		local splitGuess = left + right
		if math.abs(guess - splitGuess) < tolerance or math.abs(right) < tolerance then
			return splitGuess
		else
			return left
				+ NumericalMethods.integration.adaptiveLaguerre(
					f,
					a + stride,
					negative,
					tolerance,
					right,
					2 * stride,
					strideTolerance
				)
		end
	end
end

function NumericalMethods.integration.trapezoid(
	f: (number) -> number,
	a: number,
	b: number,
	tolerance: number,
	guess: number,
	table: { [number]: number }
): number
	local width = b - a
	if table then
		local leftVal = table[a]
		if not leftVal then
			leftVal = f(a)
			table[a] = leftVal
		end
		local rightVal = table[b]
		if not rightVal then
			rightVal = f(b)
			table[b] = rightVal
		end
		return width * (leftVal + rightVal) / 2
	else
		local leftVal, rightVal = f(a), f(b)
		return width * (leftVal + rightVal) / 2
	end
end

function NumericalMethods.integration.adaptiveTrapezoid(
	f: (number) -> number,
	a: number,
	b: number,
	tolerance: number,
	guess: number,
	table: { [number]: number }
): number
	tolerance = tolerance or 10 ^ -12
	table = table or {}
	local width = b - a
	if not guess then
		local leftVal = table[a]
		local rightVal = table[b]
		if leftVal then
			if rightVal then
				guess = width * (leftVal + rightVal) / 2
			else
				rightVal = f(b)
				table[b] = rightVal
				guess = width * (leftVal + rightVal) / 2
			end
		else
			leftVal = f(a)
			table[a] = leftVal
			if rightVal then
				guess = width * (leftVal + rightVal) / 2
			else
				rightVal = f(b)
				table[b] = rightVal
				guess = width * (leftVal + rightVal) / 2
			end
		end
	end
	local midPoint = (b + a) / 2
	local left = NumericalMethods.integration.trapezoid(f, a, midPoint, table)
	local right = NumericalMethods.integration.trapezoid(f, midPoint, b, table)
	local splitGuess = left + right
	if math.abs(guess - splitGuess) < tolerance or width < tolerance then
		return splitGuess
	else
		return NumericalMethods.integration.adaptiveTrapezoid(f, a, midPoint, tolerance, left, table)
			+ NumericalMethods.integration.adaptiveTrapezoid(f, midPoint, b, tolerance, right, table)
	end
end

return NumericalMethods
