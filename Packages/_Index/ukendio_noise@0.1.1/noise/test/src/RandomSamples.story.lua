local ReplicatedStorage = game:GetService("ReplicatedStorage")
local StarterGui = game:GetService("StarterGui")

local Noise = require(ReplicatedStorage.Noise)
local Rand = Random.new()

local SAMPLES = 1e5 -- 1e6 hangs Studio and takes about 10 seconds

return function()
	local min = math.huge
	local max = -math.huge
	local mean = 0
	local samples = table.create(SAMPLES)
	local bins = table.create(10, 0)

	for i = 1, SAMPLES do
		local value = Noise.Perlin4D(
			Rand:NextNumber(-100, 100),
			Rand:NextNumber(-100, 100),
			Rand:NextNumber(-100, 100),
			Rand:NextNumber(-100, 100)
		)
		samples[i] = value

		min = math.min(min, value)
		max = math.max(max, value)
		mean += value

		local ok, message = pcall(function()
			bins[math.floor(value * 10) + 1] += 1
		end)
		if not ok then
			print("Outside of [0, 1]: ", value)
		end
	end

	mean /= SAMPLES

	local stdDeviation = 0
	for _, sample in samples do
		stdDeviation += (sample - mean) ^ 2
	end
	stdDeviation = math.sqrt(stdDeviation / SAMPLES)

	print("Samples: ", SAMPLES)
	print("Min: ", min)
	print("Max: ", max)
	print("Mean: ", mean)
	print("Std Deviation: ", stdDeviation)
	print("Bins: ", bins)
end
