--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

--[[
   Welcome to RbxMath's Fast Fourier Transform package! The fast fourier transform
   is one of the most powerful tools in the arsenal of anyone looking for fast and
   stable numerical methods. Its uses are, unfortunately, hidden behind years of
   experience with linear algebra, mathematical analysis, or signal processing.
   Once you have uncovered its secrets, you will be overjoyed with the applications!
   
   Strictly speaking, the fast fourier transform is actually the fast discrete
   fourier transform. In particular, the discrete fourier transform (DFT) can be
   naively computed in O(n^2) operations. The fast fourier transform refers to a
   class of algorithms that implement this procedure in O(n log n) operations.

   This package will prioritize speed over usability. FastFourierTransform:FFT is
   the most user friendly function in the package, but it is not as efficient. The
   best workflow for applying the FFT is to first run FFTSetup. This will produce
   two arrays: twiddles and work. These will serve as the workspace for the FFT
   algorithm of size n. You can reuse these workspace vectors every time you apply
   an FFT of size n.
   
   FFTFormatInput will take a list of numbers (real or complex) and convert them 
   into a properly formated input. While the complex numbers you can find in Scalars
   are useful, they are slow. So, internally, all arrays of n complex numbers are
   stored as length 2n arrays of real numbers.
--]]

local Packages = script.Parent
local Complex = require(Packages.Scalars).Complex
local Tools = require(Packages.Tools)
local Primes = require(Packages.Primes)

local factorLists = {}
setmetatable(factorLists, { __mode = "v" })

local FastFourierTransform = {}

--[[
   +-----------+
   | Utilities |
   +-----------+
--]]

function FastFourierTransform.FFTSetup(n)
   local twiddles = table.create(2 * n, 0)
   local rootOfUnity = 2 * math.pi / n
   for i = 0, n - 1 do
      local twoI = 2 * i
      twiddles[twoI + 1] = math.cos(rootOfUnity * i)
      twiddles[twoI + 2] = math.sin(rootOfUnity * i)
   end
   return twiddle, table.clone(twiddles)
end

function FastFourierTransform.FFTFormatInput(list, scalarQ)
   scalarQ = scalarQ or true
   local n = #list
   local data
   if scalarQ then
      data = table.create(2 * n, 0)
      if type(list[1] == "number") then
	 for i = 0, n - 1 do
	    data[2 * i + 1] = list[i + 1]
	 end
      else
	 for i = 0, n - 1 do
	    data[2 * i + 1] = list[i + 1][1]
	    data[2 * i + 2] = list[i + 1][2]
	 end
      end
   end
   return data
end

--[[
   +-------------------------+
   | User-friendly Functions |
   +-------------------------+
--]]

function FastFourierTransform.FFT(xList, scalarQ)
   scalarQ = scalarQ or true
   local n = #xList
   local primeDecomp = Primes:decompose(n)
   local data
   if scalarQ then
      data = table.create(2 * n, 0)
      if type(xList[1] == "number") then
	 for i = 0, n - 1 do
	    data[2 * i + 1] = xList[i + 1]
	 end
      else
	 for i = 0, n - 1 do
	    data[2 * i + 1] = xList[i + 1][1]
	    data[2 * i + 2] = xList[i + 1][2]
	 end
      end
   else
      data = table.create(n, 0)
      n /= 2
   end
   FastFourierTransform.FFTS(data, n, primeDecomp)
   if scalarQ then
      local list = table.create(n, Complex:new(0))
      for i = 0, n - 1 do
	 list[i + 1] = Complex:new({data[2 * i + 1], data[2 * i + 2]})
      end
      return list
   else
      return data
   end
end

--[[
   +-------------------------------------------+
   | Implementations that need precomputations |
   +-------------------------------------------+
--]]

--[[
   This is the naive implementation of the DFT. It has an O(n^2) runtime.
--]]
function FastFourierTransform.DFTS(list, n, twiddles, work)
   local twoN = 2 * n
   work = work or table.clone(list)
   if not twiddles then
      twiddles = table.create(twoN, 0)
      local rootOfUnity = 2 * math.pi / n
      for i = 0, n - 1 do
	 local twoI = 2 * i
	 twiddles[twoI + 1] = math.cos(rootOfUnity * i)
	 twiddles[twoI + 2] = math.sin(rootOfUnity * i)
      end
   end
   for i = 0, n - 1 do
      local sumRe = 0
      local sumIm = 0
      for j = 0, n - 1 do
	 local twoIndex = 2 * ((i * j) % n)
	 local twiddleRe = twiddles[twoIndex + 1]
	 local twiddleIm = twiddles[twoIndex + 2]
	 local twoJ = 2 * j
	 local listRe = list[twoJ + 1]
	 local listIm = list[twoJ + 2]
	 sumRe += twiddlesRe * listRe - twiddlesIm * listIm
	 sumIm += twiddlesIm * listRe + twiddlesRe * listIm
      end
      local twoI = 2 * i
      work[twoI + 1] = sumRe
      work[twoI + 2] = sumIm
   end
   return work
end

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

function FastFourierTransform.FFTS(list, n, primeDecomp, twiddles, work)
   -- This FFT implementation is designed to avoid recursion
   -- We don't even need a stack to implement the FFT, but we
   -- do need a tree to keep track of where everything is. This
   -- tree will be tracked as a multi-base number. This is not
   -- memory safe!
   if n == 1 then
      return list
   end
   local twoN = 2 * n
   work = work or table.clone(list)
   if not twiddles then
      twiddles = table.create(twoN, 0)
      local rootOfUnity = 2 * math.pi / n
      for i = 0, n - 1 do
	 local twoI = 2 * i
	 twiddles[twoI + 1] = math.cos(rootOfUnity * i)
	 twiddles[twoI + 2] = math.sin(rootOfUnity * i)
      end
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
	    local along = -1
	    for k2 = 0, N2 - 1 do
	       local iter = N2 * k1 + k2
	       local sumRe = 0
	       local sumIm = 0
	       local rootIndex = numOfLeaves * iter
	       for n1, tempStart in ipairs(startList) do
		  local tempIndex = tempStart + along
		  local twiddleRe = twiddles[2 * ((rootIndex * (n1 - 1)) % n) + 1]
		  local twiddleIm = twiddles[2 * ((rootIndex * (n1 - 1)) % n) + 2]
		  local re = list[2 * tempIndex + 1]
		  local im = list[2 * tempIndex + 2]
		  sumRe += twiddleRe * re - twiddleIm * im
		  sumIm += re * twiddleIm + im * twiddleRe
	       end
	       work[2 * (start + along) + 1] = sumRe
	       work[2 * (start + along) + 2] = sumIm
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

function FastFourierTransform.FFTS_Complex(list, n, primeDecomp, twiddles, work)
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
