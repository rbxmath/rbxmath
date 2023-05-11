--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require("../src/Tools")

local Scalars = {}

local Complex = {}

function Complex:new(a, b)
	if type(a) == "table" then
		setmetatable(a, self)
		self.__index = self
		return a
	elseif b == nil then
		local result = { a, 0 }
		setmetatable(result, self)
		self.__index = self
		return result
	else
		local result = { a, b }
		setmetatable(result, self)
		self.__index = self
		return result
	end
end

function Complex.__newindex(t, k, v)
	error("Complex numbers are immutable!")
end

function Complex.__eq(left, right)
	return left[1] == right[1] and left[2] == right[2]
end

function Complex.__add(left, right)
	return Complex:new(left[1] + right[1], left[2] + right[2])
end

function Complex.__sub(left, right)
	return Complex:new(left[1] - right[1], left[2] - right[2])
end

function Complex.__mul(left, right)
	return Complex:new(left[1] * right[1] - left[2] * right[2], left[1] * right[2] + left[2] * right[1])
end

function Complex.__unm(complex)
	return Complex:new(-complex[1], -complex[2])
end

function Complex:conjugate()
	return Complex:new(self[1], -self[2])
end

function Complex:norm()
	return math.sqrt(self[1] ^ 2 + self[2] ^ 2)
end

function Complex:scale(lambda)
	return Complex:new(lambda * self[1], lambda * self[2])
end

function Complex:inverse()
	return self:conjugate():scale(1 / self:norm() ^ 2)
end

function Complex:exp(theta)
	return self:new(math.cos(theta), math.sin(theta))
end

function Complex.__div(left, right)
	return left * right:inverse()
end

function Complex:pow(n)
	local product = Complex:new(1)
	for i = 1, n do
		product = product * self
	end
	return product
end

local _signString = function(x)
	if x >= 0 then
		return " + "
	else
		return " - "
	end
end

function Complex.__tostring(complex)
	local sign = function(x)
		if x >= 0 then
			return " + "
		else
			return " - "
		end
	end
	return tostring(complex[1]) .. _signString(complex[2]) .. "i " .. tostring(math.abs(complex[2]))
end

function Complex:sort(array, order)
	local comp = function(left, right)
		if order == "d" then
			if type(left) == "number" and type(right) == "number" then
				return math.abs(left) > math.abs(right)
			elseif getmetatable(left) == getmetatable(right) then
				return left:norm() > right:norm()
			elseif getmetatable(left) == Complex then
				return left:norm() > math.abs(right)
			else
				return math.abs(left) > right:norm()
			end
		else
			if type(left) == "number" and type(right) == "number" then
				return math.abs(left) < math.abs(right)
			elseif getmetatable(left) == getmetatable(right) then
				return left:norm() < right:norm()
			elseif getmetatable(left) == Complex then
				return left:norm() < math.abs(right)
			else
				return math.abs(left) < right:norm()
			end
		end
	end

	table.sort(array, comp)
end

local _rational = {}

local _rationalFromArray = function(array)
	if array[2] == 0 then
		error("Division by zero!", -1)
	end

	local result = { math.tointeger(array[1]), math.tointeger(array[2]) }

	setmetatable(result, _rational)

	return result
end

local _rationalFromTwoInputs = function(a, b)
	if b == 0 then
		error("Division by zero!", -1)
	end

	local result = { math.tointeger(a), math.tointeger(b) }

	setmetatable(result, _rational)

	return result
end

function _rational.__eq(left, right)
	return left[1] == right[1] and left[2] == right[2]
end

function _rational.__add(left, right)
	local denominator = left[2] * right[2]
	local numerator = left[1] * right[2] + right[1] * left[2]

	local divisor = Tools.integers.gcd(denominator, numerator)

	return _rationalFromTwoInputs(numerator / divisor, denominator / divisor)
end

function _rational.__sub(left, right)
	local denominator = left[2] * right[2]
	local numerator = left[1] * right[2] - right[1] * left[2]

	local divisor = Tools.integers.gcd(denominator, numerator)

	return _rationalFromTwoInputs(numerator / divisor, denominator / divisor)
end

function _rational.__mul(left, right)
	return _rationalFromTwoInputs(left[1] * right[1], left[2] * right[2])
end

function _rational.__div(left, right)
	return _rationalFromTwoInputs(left[1] * right[2], left[2] * right[1])
end

function _rational.__tostring(rational)
	return tostring(rational[1]) .. "/" .. tostring(rational[2])
end

local _rationalFromContinuedFraction = function(continuedFraction)
	local n = #continuedFraction
	local rational = _rationalFromTwoInputs(1, continuedFraction[n])

	for i = n - 1, 1, -1 do
		rational = _rationalFromTwoInputs(continuedFraction[i], 1) + _rationalFromTwoInputs(rational[2], rational[1])
	end

	return rational
end

local _rationalApproximation = function(float, tolerance)
	local n = -math.floor(math.log(tolerance, 10) + 0.5)

	if n <= 0 then
		return math.floor(float)
	end

	local continuedFraction = Tools.continuedFractions.compute(float, n)

	local i = 1

	while math.abs(float - Tools.continuedFractions.evaluate(continuedFraction)) > tolerance and i < math.maxinteger do
		continuedFraction = Tools.continuedFractions.compute(float, n + i)
		i = i + 1
	end

	return _rationalFromContinuedFraction(continuedFraction)
end

local _ZmodP = {}

local _ZmodPFromArray = function(array)
	local result = { array[1] % array[2], array[2] }

	setmetatable(result, _ZmodP)

	return result
end

local _ZmodPFromTwoInputs = function(a, b)
	local result = { a % b, b }

	setmetatable(result, _ZmodP)

	return result
end

function _ZmodP.__eq(left, right)
	return left[1] == right[1] and left[2] == right[2]
end

function _ZmodP.__add(left, right)
	if left[2] ~= right[2] then
		error("Modulus mismatch!", -1)
	end

	return _ZmodPFromTwoInputs(left[1] + right[1], left[2])
end

function _ZmodP.__sub(left, right)
	if left[2] ~= right[2] then
		error("Modulus mismatch!", -1)
	end

	return _ZmodPFromTwoInputs(left[2] + left[1] - right[1], left[2])
end

function _ZmodP.__mul(left, right)
	if left[2] ~= right[2] then
		error("Modulus mismatch!", -1)
	end

	return _ZmodPFromTwoInputs(left[1] * right[1], left[2])
end

function _ZmodP.__tostring(zinteger)
	return tostring(zinteger[1])
end

--- Credit to Wikipedia https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
local _ZmodPInverse = function(zinteger)
	local t = 0
	local r = zinteger[2]
	local newt = 1
	local newr = zinteger[1]

	local quotient
	while newr ~= 0 do
		quotient = r % newr
		t, newt = newt, t - quotient * newt
		r, newr = newr, r - quotient * newr
	end

	if r > 1 then
		error("Input is not invertible!", -1)
	end

	return _ZmodPFromTwoInputs(t + zinteger[2], zinteger[2])
end

function _ZmodP.__div(left, right)
	return left * _ZmodPInverse(right)
end

Scalars.Complex = Complex

Scalars.Rational = {}

function Scalars.Rational.new(a, b)
	return _rationalFromTwoInputs(a, b)
end

function Scalars.Rational.newFromArray(array)
	return _rationalFromArray(array)
end

function Scalars.Rational.toReal(rational)
	return rational[1] / rational[2]
end

function Scalars.Rational.rationalApproximation(float, tolerance)
	return _rationalApproximation(float, tolerance)
end

function Scalars.Rational.rationalFromContinuedFraction(continuedFraction)
	return _rationalFromContinuedFraction(continuedFraction)
end

function Scalars.Rational.inverse(rational)
	return _rationalFromTwoInputs(rational[2], rational[1])
end

Scalars.ZmodP = {}

function Scalars.ZmodP.new(a, b)
	return _ZmodPFromTwoInputs(a, b)
end

function Scalars.ZmodP.newFromArray(array)
	return _ZmodPFromArray(array)
end

function Scalars.ZmodP.inverse(zinteger)
	return _ZmodPInverse(zinteger)
end

return Scalars
