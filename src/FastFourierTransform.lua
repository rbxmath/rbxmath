local Complex = require("src.Scalars").Complex
local Tools = require("src.Tools")
local Primes = require("src.Primes")

local factorLists = {}
setmetatable(factorLists, {__mode = "v"})

local FastFourierTransform = {
    factorLists = factorLists
}

-- Code taken from https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

function FastFourierTransform:FFT2 (xList, realQ)
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
            for j = 1, m/2 do
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

function FastFourierTransform:IFFT2 (xList, realQ)
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
            for j = 1, m/2 do
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

function FastFourierTransform:FFT (xList, realQ)
    local data = {}
    local n = #xList
    local primeDecomp = Primes:decompose(n)
    if realQ then
        for i = n, 1, -1 do
            data[i] = Complex:new(xList[i])
        end
    end
    local logN = #primeDecomp
    local m = 1
    for s = 1, logN do
        m = m * primeDecomp[s]
        local omegaM = Complex:exp(-2 * math.pi / m)
        for k = 0, n - 1, m do
            local omega = Complex:new(1)
            for j = 1, m / primeDecomp[s] do
                local u = {}
                for l = 0, primeDecomp[s] - 1 do
                    u[l + 1] = omega:pow(l) * data[k + j + l * m / primeDecomp[s]]
                end
                for l = 0, primeDecomp[s] - 1 do
                    local sum = Complex:new(0)
                    for ll = 0, primeDecomp[s] - 1 do
                        sum = sum + Complex:exp(-2 * math.pi * ll * j / primeDecomp[s])
                    end
                    data[k + j + l * m / primeDecomp[s]] = data[k + j + l * m / primeDecomp[s]] + sum
                end
                omega = omega * omegaM
            end
        end
    end
    return data
end

return FastFourierTransform