local ReplicatedStorage = game:GetService("ReplicatedStorage")
local RunService = game:GetService("RunService")

local Spline = require(ReplicatedStorage.Spline)
local gizmo = require(ReplicatedStorage.gizmo)

-- local N = 10
-- local rand = Random.new()
-- local points = table.create(N)
-- for i = 1, N do
-- 	points[i] = Vector3.new(rand:NextNumber(-10, 10), rand:NextNumber(-10, 10), rand:NextNumber(-10, 10))
-- end

-- local spline = Spline.CubicBezier.new(points)

-- RunService.Heartbeat:Connect(function()
-- 	local positions = spline:SolveMany(spline.SolvePosition, 100)
-- 	for _, pos in positions do
-- 		gizmo.drawPoint(pos)
-- 	end

-- 	-- for i = 1, 100 do
-- 	-- 	local t = (i - 1) / 100
-- 	-- 	gizmo.drawPoint(spline:SolvePosition(t))
-- 	-- end
-- end)

RunService.Heartbeat:Connect(function()
	local parts = workspace.Model:GetChildren()
	local curve = Spline.CubicHermite.new(parts[1].Position, parts[2].Position, parts[3].Position, parts[4].Position)

	for i = 0, 100 do
		local pos = curve:SolvePosition(i / 100)
		gizmo.drawPoint(pos)
	end
end)
