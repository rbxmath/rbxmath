local ReplicatedStorage = game:GetService("ReplicatedStorage")
local StarterGui = game:GetService("StarterGui")

local Noise = require(ReplicatedStorage.Noise)
local Rand = Random.new()

local WIDTH = 64
local HEIGHT = 64
local PIXEL_SIZE = 4

return function(target)
	local tile = Noise.GetNoiseTile(
		WIDTH,
		HEIGHT,
		{ Rand:NextNumber(0, 100), Rand:NextNumber(0, 100), Rand:NextNumber(0, 100), Rand:NextNumber(0, 100) },
		Rand:NextNumber(0.4, 2)
	)

	local frame1 = Instance.new("Frame")
	frame1.Size = UDim2.fromOffset(WIDTH, HEIGHT)
	frame1.BackgroundTransparency = 1

	for x = 1, WIDTH do
		for y = 1, HEIGHT do
			local v = tile[x][y]

			local pixel = Instance.new("Frame", frame1)
			pixel.Size = UDim2.fromOffset(PIXEL_SIZE, PIXEL_SIZE)
			pixel.Position = UDim2.fromOffset(PIXEL_SIZE * (x - 1), PIXEL_SIZE * (y - 1))
			pixel.BorderSizePixel = 0
			pixel.BackgroundColor3 = Color3.new(v, v, v)
		end
	end

	frame1.Parent = target

	local frame2 = frame1:Clone()
	frame2.Position = UDim2.fromOffset(WIDTH * PIXEL_SIZE, 0)
	frame2.Parent = target

	local frame3 = frame1:Clone()
	frame3.Position = UDim2.fromOffset(0, HEIGHT * PIXEL_SIZE)
	frame3.Parent = target

	local frame4 = frame1:Clone()
	frame4.Position = UDim2.fromOffset(WIDTH * PIXEL_SIZE, HEIGHT * PIXEL_SIZE)
	frame4.Parent = target

	-- ScreenGui.Parent = target

	return function()
		-- ScreenGui:Destroy()
		frame1:Destroy()
		frame2:Destroy()
		frame3:Destroy()
		frame4:Destroy()
	end
end
