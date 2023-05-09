--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Complex = require("src/Scalars").Complex
local Tools = require("src/Tools")
local Primes = require("src/Primes")

local factorLists = {}
setmetatable(factorLists, { __mode = "v" })

local FastFourierTransform = {
	factorLists = factorLists,
}

-- FFT2 code taken from https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

function FastFourierTransform:FFT2(xList, realQ)
	local data = {}
	local logN = math.floor((math.log(#xList) / math.log(2) + 0.5))
	local n = math.pow(2, logN)
	if realQ then
		for i = n, 1, -1 do
			if xList[i] then
				data[i] = Complex:new(xList[i])
			else
				data[i] = Complex:new(0)
			end
		end
	else
		for i = n, 1, -1 do
			data[i] = xList[i] or 0
		end
	end
	local rootN = Complex:new(math.sqrt(n))
	for s = 1, logN do
		local m = math.pow(2, s)
		local omegaM = Complex:exp(-2 * math.pi / m)
		for k = 0, n - 1, m do
			local omega = Complex:new(1)
			for j = 1, m / 2 do
				local t = omega * data[k + j + m / 2]
				local u = data[k + j]
				data[k + j] = u + t
				data[k + j + m / 2] = u - t
				omega = omega * omegaM
			end
		end
	end
	return data
end

function FastFourierTransform:IFFT2(xList, realQ)
	local data = {}
	local logN = math.floor((math.log(#xList) / math.log(2) + 0.5))
	local n = math.pow(2, logN)
	if realQ then
		for i = n, 1, -1 do
			if xList[i] then
				data[i] = Complex:new(xList[i])
			else
				data[i] = Complex:new(0)
			end
		end
	else
		for i = n, 1, -1 do
			data[i] = xList[i] or 0
		end
	end
	local rootN = Complex:new(math.sqrt(n))
	for s = 1, logN do
		local m = math.pow(2, s)
		local omegaM = Complex:exp(-2 * math.pi / m)
		for k = 0, n - 1, m do
			local omega = Complex:new(1)
			for j = 1, m / 2 do
				local t = omega * data[k + j + m / 2]
				local u = data[k + j]
				data[k + j] = u + t
				data[k + j + m / 2] = u - t
				omega = omega * omegaM
			end
		end
	end
	n = Complex:new(n)
	for i = 1, #data do
		data[i] = data[i] / n
	end
	return data
end

function FastFourierTransform:FFT(xList, realQ)
	local data = {}
	local n = #xList
	local primeDecomp = Primes:decompose(n)
	if realQ then
		for i = n, 1, -1 do
			data[i] = Complex:new(xList[i])
		end
	end
	return FastFourierTransform:FFTS(data, n, 1, 1, primeDecomp)
end

function FastFourierTransform:FFTS(xList, n, startIndex, stride, primeDecomp)
	if n == 1 then
		return { xList[startIndex] }
	end
	local N1 = primeDecomp[#primeDecomp]
	local N2 = n / N1
	local newStride = stride * N1
	local u = {}
	for i = 0, N1 - 1 do
		u[i + 1] = FastFourierTransform:FFTS(
			xList,
			N2,
			startIndex + i * stride,
			newStride,
			Tools.list.sublist(primeDecomp, #primeDecomp - 1)
		)
	end
	local data = {}
	for k1 = 0, N1 - 1 do
		for k2 = 0, N2 - 1 do
			local iter = N2 * k1 + k2
			local sum = Complex:new(0)
			for n1 = 0, N1 - 1 do
				sum = sum + Complex:exp(-2 * math.pi * n1 * iter / n) * u[n1 + 1][k2 + 1]
			end
			data[iter + 1] = sum
		end
	end
	return data
end

function FastFourierTransform:FCT(xList)
	local newList = Tools.list.copy(xList)
	local n = #xList
	for i = n - 1, 2, -1 do
		newList[#newList + 1] = xList[i]
	end
	local length = #newList
	local FFT = FastFourierTransform:FFT(newList, true)
	local data = {}
	data[1] = FFT[1][1] / length
	data[n] = FFT[n][1] / length
	for i = 2, n - 1 do
		data[i] = (FFT[i][1] + FFT[length - i + 2][1]) / length
	end
	return data
end

return FastFourierTransform
