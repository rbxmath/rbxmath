--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local PrimeList = require(script.Parent.PrimeList)

local Primes = {
	List = PrimeList,
}

function Primes.BuildPrimeList(n: number): { number }
	for i = PrimeList[#PrimeList] + 1, n do
		local isPrime = true
		local sqrt = math.sqrt(i)

		for k, v in PrimeList do
			if k > sqrt then
				break
			end

			isPrime = (i % v ~= 0)
			if not isPrime then
				break
			end
		end

		if isPrime then
			table.insert(PrimeList, i)
		end
	end

	return PrimeList
end

function Primes.Decompose(n: number): { number }
	if math.sqrt(n) > PrimeList[#PrimeList] then
		Primes.BuildPrimeList(n)
	end

	local decomp = {}
	for k, v in PrimeList do
		if v > n then
			break
		end

		local temp = n
		while temp % v == 0 do
			temp /= v
			table.insert(decomp, v)
		end
	end

	return decomp
end

return Primes
