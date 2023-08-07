--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require(script.Parent.Tools)
local MA = require(script.Parent.Matrices)

local EPSILON = 10 ^ -13

local NumericalMethods = {}

NumericalMethods.solvers = {}

--[=[
	Attempts to solve the problem f(x) = t for x between left and right bounds
	given by a and b respectively using Regula Falsi. This method only works if
	f(a) and f(b) have different signs.

	@param f --- A differentiable function
	@param t --- The right hand value in the equation f(x) = t
	@param a --- Left bound of the search interval
	@param b --- Right bound of the search interval
	@param tol --- An optional parameter to specify a tolerance for root finding.
--]=]
function NumericalMethods.solvers.regulaFalsi(
	f: (number) -> number,
	t: number,
	a: number,
	b: number,
	tol: number
): number
	local tolerance = tol or EPSILON
	local leftValue, rightValue, middleValue = f(a) - t, f(b) - t, 0
	local left, right, middle = a, b, (a + b) * 0.5
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

--[=[
	Attempts to solve the problem f(x) = t for x using Newton's method.
	This might take a long time to run if the function is bad.

	@param f --- A differentiable function
	@param fPrime --- The derivative of the function f
	@param t --- The right hand value in the equation f(x) = t
	@param x --- An optional initial guess at a solution. Highly recommended, but not needed.
	@param tol --- An optional parameter to specify a tolerance for root finding.
--]=]
function NumericalMethods.solvers.newtonsMethod(
	f: (number) -> number,
	fPrime: (number) -> number,
	t: number,
	x: number,
	tol: number
): number
	tol = tol or EPSILON
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

--[=[
	Attempts to solve the problem f(x) = t using a number of methods.

	@param f --- A differentiable function
	@param t --- The right hand value in the equation f(x) = t
	@param x --- An optional initial guess at a solution. Highly recommended, but not needed.
	@param tolerance --- An optional parameter to specify a tolerance for root finding.
--]=]
function NumericalMethods.solvers.solve(f: (number) -> number, t: number, x: number, tolerance: number): number
	tolerance = tolerance or EPSILON
	x = x or 0
	local a = x - 10 * tolerance
	local b = x + 10 * tolerance
	local leftValue = f(a) - t
	local rightValue = f(b) - t
	while leftValue >= tolerance do
		if leftValue * rightValue > 0 then
			local middleValue = 0
			local left, right, middle = a, b, (a + b) * 0.5
			if math.abs(leftValue) < tolerance then
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
		else
			b, a = (b * leftValue - a * rightValue) / (leftValue - rightValue), b
			rightValue, leftValue = f(b) - t, rightValue
		end
	end
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

--[=[
	Integrates a function using 5-point Gaussian quadrature. Legendre roots and
	quadrature weights sourced from https://pomax.github.io/bezierinfo/legendre-gauss.html
	with highest precision afforded by 64-bit floats.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
--]=]
function NumericalMethods.integration.fivePointGaussianQuadrature(f: (number) -> number, a: number, b: number): number
	if a == b then
		return 0
	end

	local sign = 1

	if b < a then
		b, a = a, b
		sign = -1
	end

	local halfInterval = (b - a) * 0.5
	local midpoint = (b + a) * 0.5

	local x2 = halfInterval * 0.5384693101056831
	local x3 = halfInterval * 0.906179845938664

	return sign
		* halfInterval
		* (
			0.5688888888888889 * f(midpoint)
			+ 0.47862867049936647 * (f(midpoint - x2) + f(midpoint + x2))
			+ 0.23692688505618908 * (f(midpoint - x3) + f(midpoint + x3))
		)
end

--[=[
	Integrates a function, f, from point a to point b using adaptive quadrature.
	Computes the integral by recursively subdividing the interval until guess is
	within tolerance of the subdivided guess.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
	@param tolerance --- An optional parameter to set a numerical tolerance
	@param guess --- An optional parameter to guess the integral. Mostly for internal use.
--]=]
function NumericalMethods.integration.adaptiveQuadrature(
	f: (number) -> number,
	a: number,
	b: number,
	tolerance: number,
	guess: number
): number
	tolerance = tolerance or EPSILON
	guess = guess or NumericalMethods.integration.fivePointGaussianQuadrature(f, a, b)
	if math.abs(b - a) < tolerance then
		return guess
	end
	local midPoint = (b + a) * 0.5
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

--[=[
	Integrates a function, f, from point a to point b using adaptive quadrature.
	Computes the integral by recursively subdividing the interval until guess is
	within tolerance of the subdivided guess.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
	@param tolerance --- An optional parameter to set a numerical tolerance
	@param guess --- An optional parameter to guess the integral. Mostly for internal use.
--]=]
NumericalMethods.integration.integrate = NumericalMethods.integration.adaptiveQuadrature

NumericalMethods.integration.fivePointLaguerreGrid =
	{ 0.2635603197181409, 1.4134030591065208, 3.5964257710407237, 7.085810005858842, 12.640800844275782 }
NumericalMethods.integration.fivePointLaguerreWeights =
	{ 0.5217556105828091, 0.3986668110831633, 0.07594244968170724, 0.003611758679922023, 0.000023369972385776235 }

--[=[
	Integrates a function, f, from point a to plus or minus infinity using
	five point Laguerre quadrature. Only works well on functions that decrease
	exponentially.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param negative --- A boolean value that determines the sign of infinity
	@param tolerance --- An optional parameter to set a numerical tolerance
--]=]
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

--[=[
	Integrates a function, f, from point a to plus of minus infinity using
	adaptive five point Laguerre quadrature. Only works well on functions of
	exponential decrease.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
	@param tolerance --- An optional parameter to set a numerical tolerance
	@param guess --- An optional parameter to guess the integral. Mostly for internal use.
--]=]
function NumericalMethods.integration.adaptiveLaguerre(
	f: (number) -> number,
	a: number,
	negative: boolean,
	tolerance: number,
	guess: number,
	stride: number,
	strideTolerance: number
): number
	tolerance = tolerance or EPSILON
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

--[=[
	Integrates a function, f, from point a to point b using the trapezoid rule. This
	should only be used as a last result for integration, or in the case where evaluating
	the function is very expensive.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
	@param table --- An optional parameter that allows you to pass in the function values. Mostly for internal use.
--]=]
function NumericalMethods.integration.trapezoid(
	f: (number) -> number,
	a: number,
	b: number,
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
		return width * (leftVal + rightVal) * 0.5
	else
		local leftVal, rightVal = f(a), f(b)
		return width * (leftVal + rightVal) * 0.5
	end
end

--[=[
	Integrates a function, f, from point a to point b using adaptive trapezid rule.
	Computes the integral by recursively subdividing the interval until guess is
	within tolerance of the subdivided guess. Only use this as a last result or if
	your function is very expoensive to evaluate.

	@param f --- An integrable function
	@param a --- The lower integral bound
	@param b --- The upper integral bound
	@param tolerance --- An optional parameter to set a numerical tolerance
	@param guess --- An optional parameter to guess the integral. Mostly for internal use.
	@param table --- An optional parameter that allows you to pass in the function values. Mostly for internal use.
--]=]
function NumericalMethods.integration.adaptiveTrapezoid(
	f: (number) -> number,
	a: number,
	b: number,
	tolerance: number,
	guess: number,
	table: { [number]: number }
): number
	tolerance = tolerance or EPSILON * 10
	table = table or {}
	local width = b - a
	if not guess then
		local leftVal = table[a]
		local rightVal = table[b]
		if leftVal then
			if rightVal then
				guess = width * (leftVal + rightVal) * 0.5
			else
				rightVal = f(b)
				table[b] = rightVal
				guess = width * (leftVal + rightVal) * 0.5
			end
		else
			leftVal = f(a)
			table[a] = leftVal
			if rightVal then
				guess = width * (leftVal + rightVal) * 0.5
			else
				rightVal = f(b)
				table[b] = rightVal
				guess = width * (leftVal + rightVal) * 0.5
			end
		end
	end
	local midPoint = (b + a) * 0.5
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
