local FFT = require("../src/FastFourierTransform")
local Scalars = require("../src/Scalars")
local Tools = require("../src/Tools")
local Complex = Scalars.Complex
local Primes = require("../src/Primes")

print(#Primes.primes)
print(Tools.list.tostring(Primes:decompose(15)))
local data = { 1, 1, 1, 1 }
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 1, 1, 1, 1, 1, 1 }
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 1, -1, 1, -1 }
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 10., -5.41421, 2., -2.58579, 2., -2.58579, 2., -5.41421 }
print(Tools.list.tostring(FFT:IFFT(data, true)))

--[[local dataList = {}
for i = 1, 15 do
    local data = {}
    for j = 1, math.pow(2, i) do
        data[j] = 1
    end
    dataList[i] = data
end]]

--[[local tic
local save
for i = 1, 10 do
    tic = os.clock()
    save = FFT:FFT2(dataList[i], true)
    print(os.clock() - tic)
end]]
