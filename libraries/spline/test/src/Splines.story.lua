local ReplicatedStorage = game:GetService("ReplicatedStorage")

local Spline = require(ReplicatedStorage.Packages.Spline)
local gizmo = require(ReplicatedStorage.gizmo)

return function(target)
	local N = 8

	local parts = table.create(N)

	for i = 1, N do
		local part = Instance.new("Part", workspace)
		part.Anchored = true
		part.CanCollide = false
		part.Size = Vector3.one * 0.5
		part.Transparency = 0.5
		part.BrickColor = BrickColor.new("Persimmon")
		part.CFrame = CFrame.new(0, 10, i)
		parts[i] = part
	end

	local heartbeat = game:GetService("RunService").Heartbeat:Connect(function()
		local points = {}
		local velocities = {}
		for i, part in parts do
			if i % 2 == 0 then
				table.insert(points, part.Position)
			else
				table.insert(velocities, part.Position)
			end
		end

		local spline = Spline.CubicHermite.Spline.new(points, velocities)
		spline:ToUnitSpeed()

		gizmo.style.color = Color3.new(1, 1, 1)
		for i = 0, 100 do
			gizmo.point.draw(spline:SolvePosition(i / 100))
		end
	end)

	return function()
		heartbeat:Disconnect()
		for _, part in parts do
			part:Destroy()
		end
	end
end
