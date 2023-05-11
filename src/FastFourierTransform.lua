--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Complex = require("../src/Scalars").Complex
local Tools = require("../src/Tools")
local Primes = require("../src/Primes")

local factorLists = {}
setmetatable(factorLists, { __mode = "v" })

local FastFourierTransform = {
	factorLists = factorLists,
}

function FastFourierTransform:FFT(xList, realQ)
	local data = {}
	local n = #xList
	local primeDecomp = Primes:decompose(n)
	if realQ then
		for i = n, 1, -1 do
			data[i] = Complex:new(xList[i])
		end
	else
		data = xList
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
				sum = sum + Complex:exp(2 * math.pi * n1 * iter / n) * u[n1 + 1][k2 + 1]
			end
			data[iter + 1] = sum
		end
	end
	return data
end

function FastFourierTransform:IFFT(xList, realQ)
	local data = {}
	local n = #xList
	local primeDecomp = Primes:decompose(n)
	if realQ then
		for i = n, 1, -1 do
			data[i] = Complex:new(xList[i])
		end
	else
		data = xList
	end
	return FastFourierTransform:IFFTS(data, n, 1, 1, primeDecomp)
end

function FastFourierTransform:IFFTS(xList, n, startIndex, stride, primeDecomp)
	if n == 1 then
		return { xList[startIndex] }
	end
	local N1 = primeDecomp[#primeDecomp]
	local N2 = n / N1
	local newStride = stride * N1
	local u = {}
	for i = 0, N1 - 1 do
		u[i + 1] = FastFourierTransform:IFFTS(
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

function FastFourierTransform:FCT1(xList)
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
		data[i] = FFT[i][1] / length
	end
	return data
end

function FastFourierTransform:IFCT1(xList)
	local newList = Tools.list.copy(xList)
	local n = #xList
	for i = n - 1, 2, -1 do
		newList[#newList + 1] = xList[i]
	end
	local length = #newList
	local FFT = FastFourierTransform:IFFT(newList, true)
	local data = {}
	data[1] = FFT[1][1]
	data[n] = FFT[n][1]
	for i = 2, n - 1 do
		data[i] = FFT[i][1]
	end
	return data
end

function FastFourierTransform:FCT(xList)
	local newList = {}
	local n = #xList
	for i = 1, n do
		newList[2 * (i - 1) + 1] = 0
		newList[2 * i] = xList[i]
	end
	newList[2 * n + 1] = 0
	for i = 0, 2 * n - 1 do
		newList[4 * n - i] = newList[i + 2]
	end
	local FFT = FastFourierTransform:FFT(newList, true)
	local data = {}
	for i = 1, n do
		data[i] = FFT[i][1] / (2 * n)
	end
	return data
end

function FastFourierTransform:IFCT(xList)
	local newList = {}
	local n = #xList
	newList[1] = xList[1]
	newList[2 * n + 1] = -xList[1]
	for i = 2, n do
		newList[i] = xList[i]
		newList[n + i] = -xList[n - i + 2]
		newList[2 * n + i] = -xList[i]
		newList[3 * n + i] = xList[n - i + 2]
	end
	newList[n + 1] = 0
	newList[3 * n + 1] = 0
	local FFT = FastFourierTransform:IFFT(newList, true)
	local data = {}
	for i = 1, n do
		data[i] = FFT[2 * i][1] / 2
	end
	return data
end

return FastFourierTransform
