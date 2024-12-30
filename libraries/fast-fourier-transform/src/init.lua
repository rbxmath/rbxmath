--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Packages = script.Parent
local Complex = require(Packages.Scalars).Complex
local Tools = require(Packages.Tools)
local Primes = require(Packages.Primes)

local factorLists = {}
setmetatable(factorLists, { __mode = "v" })

local FastFourierTransform = {
   factorLists = factorLists,
}

function FastFourierTransform:RFFT(xList)
   local data = {}
   local n = #xList
   local primeDecomp = Primes:decompose(n)
   for i = n, 1, -1 do
      data[i] = Complex:new(xList[i])
   end
   return FastFourierTransform:RFFTS(data, n, 1, 1, primeDecomp)
end

function FastFourierTransform:RFFTS(xList, n, startIndex, stride, primeDecomp)
   if n == 1 then
      return { xList[startIndex] }
   end
   local N1 = primeDecomp[#primeDecomp]
   local N2 = n / N1
   local newStride = stride * N1
   local u = {}
   for i = 0, N1 - 1 do
      u[i + 1] = FastFourierTransform:RFFTS(
	 xList,
	 N2,
	 startIndex + i * stride,
	 newStride,
	 Tools.list.sublist(primeDecomp, #primeDecomp - 1)
      )
   end
   local data = {}
   -- Because n = N1 * N2, the following two nested for loops combined run for n iterations.
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

function FastFourierTransform:FFT(xList, realQ)
   local n = #xList
   local data = table.clone(xList)
   local primeDecomp = Primes:decompose(n)
   if realQ then
      for i = n, 1, -1 do
	 data[i] = Complex:new(xList[i])
      end
   else
      data = xList
   end
   return FastFourierTransform.FFTS(data, n, primeDecomp)
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
   return FastFourierTransform.IFFTS(data, n, primeDecomp)
end

local function iterateMultiBaseString(numberString, primeList, numPrimes, level)
   numPrimes = numPrimes or #primeList
   local currentIndex = 1
   local currentPrimeIndex = level
   local currentNumber = numberString[currentIndex] + 1
   local currentPrime = primeList[currentPrimeIndex]
   while currentNumber == currentPrime and currentPrimeIndex < numPrimes do
      numberString[currentIndex] = 0
      currentIndex += 1
      currentPrimeIndex += 1
      currentNumber = numberString[currentIndex] + 1
      currentPrime = primeList[currentPrimeIndex]
   end
   numberString[currentIndex] = currentNumber
   return numberString
end

local function computeStrideAndStart(numberString, primeList, numPrimes, level)
   numPrimes = numPrimes or #primeList
   local stride = 1
   local newStart = 1
   for i = numPrimes, level, -1 do
      newStart += stride * numberString[i - level + 1]
      stride *= primeList[i]
   end
   return stride, newStart
end

local function computeStride(primeList, level, numPrimes)
   numPrimes = numPrimes or #primeList
   local stride = 1
   for i = numPrimes, level, -1 do
      stride *= primeList[i]
   end
   return stride
end

function FastFourierTransform.DFT(list, n, twiddles, work)
   work = work or table.clone(list)
   n = n or #list
   if not twiddles then
      twiddles = table.create(#list, Complex:new(0))
      local rootOfUnity = 2 * math.pi / n
      for i = 1, n do
	 twiddles[i] = Complex:exp(rootOfUnity * (i - 1))
      end
   end
   for i = 0, n - 1 do
      local sum = Complex:new(0)
      for j = 0, n - 1 do
	 local index = (i * j) % n + 1
	 sum += twiddles[index] * list[j + 1]
      end
      work[i + 1] = sum
   end
   return work
end

function FastFourierTransform.FFTS(list, n, primeDecomp, twiddles, work)
   -- This FFT implementation is designed to avoid recursion
   -- We don't even need a stack to implement the FFT, but we
   -- do need a tree to keep track of where everything is. This
   -- tree will be tracked as a multi-base number. This is not
   -- memory safe!
   work = work or table.clone(list)
   if not twiddles then
      twiddles = table.create(#list, Complex:new(0))
      local rootOfUnity = 2 * math.pi / n
      for i = 1, n do
	 twiddles[i] = Complex:exp(rootOfUnity * (i - 1))
      end
   end
   if n == 1 then
      return list
   end
   local numOfPrimes = #primeDecomp
   -- The current working prime
   local N1 = 1
   -- The size of the current working vector
   local N2 = 1
   local startList = table.create(primeDecomp[numOfPrimes])
   for level = 1, numOfPrimes do
      N2 *= N1
      N1 = primeDecomp[level]
      local zeroIndexLevel = level - 1
      local pseudoN = N1 * N2
      local locNumOfPrimes = numOfPrimes - zeroIndexLevel
      local numberString = table.create(locNumOfPrimes, 0)
      local numOfLeaves = n / pseudoN
      local stride = 1
      for i = numOfPrimes, level, -1 do
	 stride *= primeDecomp[i]
      end
      for leaf = 1, numOfLeaves do
	 for i = 1, N1 do
	    stride = 1
	    local newStart = 1
	    for i = numOfPrimes, level, -1 do
	       newStart += stride * numberString[i - zeroIndexLevel]
	       stride *= primeDecomp[i]
	    end
	    startList[i] = newStart
	    local currentIndex = 1
	    local currentPrimeIndex = level
	    local currentNumber = numberString[currentIndex] + 1
	    local currentPrime = primeDecomp[currentPrimeIndex]
	    while currentNumber == currentPrime and currentPrimeIndex < numOfPrimes do
	       numberString[currentIndex] = 0
	       currentIndex += 1
	       currentPrimeIndex += 1
	       currentNumber = numberString[currentIndex] + 1
	       currentPrime = primeDecomp[currentPrimeIndex]
	    end
	    numberString[currentIndex] = currentNumber
	 end
	 for k1 = 0, N1 - 1 do
	    local start = startList[k1 + 1]
	    local along = 0
	    for k2 = 0, N2 - 1 do
	       local iter = N2 * k1 + k2
	       local sum = Complex:new(0)
	       local rootIndex = numOfLeaves * iter
	       for n1, tempStart in ipairs(startList) do
		  sum += twiddles[(rootIndex * (n1 - 1)) % n + 1] *
		     list[tempStart + along]
	       end
	       work[start + along] = sum
	       along += stride
	    end
	 end
      end
      for k, v in ipairs(work) do
	 list[k] = v
      end
   end
   return list
end

function FastFourierTransform.IFFTS(list, n, primeDecomp)
   -- This FFT implementation is designed to avoid recursion
   -- We don't even need a stack to implement the FFT, but we
   -- do need a tree to keep track of where everything is. This
   -- tree will be tracked as a multi-base number. This is not
   -- memory safe!
   local work = table.clone(list)
   local twiddles = table.create(#list, Complex:new(0))
   if n == 1 then
      return list
   end
   local rootOfUnity = -2 * math.pi / n
   for i = 1, n do
      twiddles[i] = Complex:exp(rootOfUnity * (i - 1))
   end
   local numOfPrimes = #primeDecomp
   -- The current working prime
   local N1 = 1
   -- The size of the current working vector
   local N2 = 1
   for level = 1, numOfPrimes do
      N2 *= N1
      N1 = primeDecomp[level]
      local zeroIndexLevel = level - 1
      local pseudoN = N1 * N2
      local locNumOfPrimes = numOfPrimes - zeroIndexLevel
      local numberString = table.create(locNumOfPrimes, 0)
      local numOfLeaves = n / pseudoN
      local stride = 1
      for i = numOfPrimes, level, -1 do
	 stride *= primeDecomp[i]
      end
      local startList = table.create(N1, 0)
      for leaf = 1, numOfLeaves do
	 for i = 1, N1 do
	    stride = 1
	    local newStart = 1
	    for i = numOfPrimes, level, -1 do
	       newStart += stride * numberString[i - zeroIndexLevel]
	       stride *= primeDecomp[i]
	    end
	    startList[i] = newStart
	    local currentIndex = 1
	    local currentPrimeIndex = level
	    local currentNumber = numberString[currentIndex] + 1
	    local currentPrime = primeDecomp[currentPrimeIndex]
	    while currentNumber == currentPrime and currentPrimeIndex < numOfPrimes do
	       numberString[currentIndex] = 0
	       currentIndex += 1
	       currentPrimeIndex += 1
	       currentNumber = numberString[currentIndex] + 1
	       currentPrime = primeDecomp[currentPrimeIndex]
	    end
	    numberString[currentIndex] = currentNumber
	 end
	 for k1 = 0, N1 - 1 do
	    local start = startList[k1 + 1]
	    local along = 0
	    for k2 = 0, N2 - 1 do
	       local iter = N2 * k1 + k2
	       local sum = Complex:new(0)
	       local rootIndex = numOfLeaves * iter
	       for n1, tempStart in ipairs(startList) do
		  sum += twiddles[(rootIndex * (n1 - 1)) % n + 1] *
		     list[tempStart + along]
	       end
	       work[start + along] = sum
	       along += stride
	    end
	 end
      end
      for k, v in ipairs(work) do
	 list[k] = v
      end
   end
   return list
end

function FastFourierTransform:FFTS_OLD(xList, n, startIndex, stride, primeDecomp)
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
   local const1 = 2 * math.pi / n
   -- Because n = N1 * N2, the following two nested for loops combined run for n iterations.
   for k1 = 0, N1 - 1 do
      for k2 = 0, N2 - 1 do
	 local iter = N2 * k1 + k2
	 local sum = Complex:new(0)
	 local const2 = const1 * iter
	 for n1 = 0, N1 - 1 do
	    sum = sum + Complex:exp(const2 * n1) * u[n1 + 1][k2 + 1]
	 end
	 data[iter + 1] = sum
      end
   end
   return data
end

function FastFourierTransform:IFFTS_OLD(xList, n, startIndex, stride, primeDecomp)
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
   local const1 = -2 * math.pi / n
   for k1 = 0, N1 - 1 do
      for k2 = 0, N2 - 1 do
	 local iter = N2 * k1 + k2
	 local sum = Complex:new(0)
	 local const2 = const1 * iter
	 for n1 = 0, N1 - 1 do
	    sum = sum + Complex:exp(const2 * n1) * u[n1 + 1][k2 + 1]
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
   data[1] = FFT[1][1] / (2 * length)
   data[n] = FFT[n][1] / (2 * length)
   for i = 2, n - 1 do
      data[i] = FFT[i][1] / length
   end
   return data
end

-- DOES NOT WORK CURRENTLY
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
