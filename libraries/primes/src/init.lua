--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local PrimeList = require(script.Parent.PrimeList)

local Primes = {
    primes = PrimeList
}

function Primes:buildPrimesList (n : number) : { [number] : number }
    local primesList = self.primes or {}
    local max = primesList[#primesList]
    local start
    if max then
        start = max + 1
    else
        start = 2
    end
    local primesList = {}
    for i = start, n do
        local isPrime = true
        local sqrt = math.sqrt(i)
        for k, v in ipairs(primesList) do
            if k > sqrt then
                break
            end
            isPrime = (i % v ~= 0)
            if not isPrime then
                break
            end
        end
        if isPrime then
            primesList[#primesList + 1] = i
        end
    end
    return primesList
end

function Primes:decompose (n : number) : { [number] : number }
    local sqrt = math.sqrt(n)
    local primeList = {}
    if sqrt > self.primes[#self.primes] then
        primeList = self:buildPrimesList(n)
    else
        primeList = self.primes
    end
    local decomp = {}
    for k, v in ipairs(primeList) do
        if v > n then
            break
        end
        local temp = n
        while temp % v == 0 do
            temp = temp / v
            decomp[#decomp + 1] = v
        end
    end
    return decomp
end

return Primes