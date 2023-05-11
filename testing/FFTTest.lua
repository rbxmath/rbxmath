local FFT = require("src.FastFourierTransform")
local Scalars = require("src.Scalars")
local Tools = require("src.Tools")
local Complex = Scalars.Complex
local Primes = require("src.Primes")

print(#Primes.primes)
print(Tools.list.tostring(Primes:decompose(15)))
local data = { 1, 1, 1, 1 }
print(Tools.list.tostring(FFT:FFT2(data, true)))
print(Tools.list.tostring(FFT:IFFT2(FFT:FFT2(data, true))))
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 1, 1, 1, 1, 1, 1 }
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 1, -1, 1, -1 }
print(Tools.list.tostring(FFT:FFT(data, true)))
data = { 10., -5.41421, 2., -2.58579, 2., -2.58579, 2., -5.41421 }
print(Tools.list.tostring(FFT:FFT(data, true)))

print("\n\n\n\n\n")

print("Atempting to recover Chebyshev series 0T_0 + 1T_1 + 2T_2 + 3T_3 + 4T_4 from values on grid...")
data = { 10., -5.414213562373094, 2., -2.585786437626906, 2. }
print("Values recovered:", Tools.list.tostring(FFT:FCT(data)))

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
