--[[
This file is for use by Benchmarker (https://boatbomber.itch.io/benchmarker)

|WARNING| THIS RUNS IN YOUR REAL ENVIRONMENT. |WARNING|
--]]

local ReplicatedStorage = game:GetService("ReplicatedStorage")
local Spline = require(ReplicatedStorage.Packages.Spline)

local Rand = Random.new()
local POINTS = 100
local SOLVE_SAMPLES = 100
local N = 1e2

return {
	ParameterGenerator = function()
		local points = table.create(POINTS)

		for i = 1, POINTS do
			points[i] = vector.create(Rand:NextNumber(-100, 100), Rand:NextNumber(-100, 100), Rand:NextNumber(-100, 100))
		end

		return Spline.CubicBezier.new(points)
	end,

	Functions = {
		["SolvePosition Many"] = function(Profiler, spline)
			for _ = 1, N do
				spline:SolveMany(spline.SolvePosition, SOLVE_SAMPLES)
			end
		end,

		["SolvePosition"] = function(Profiler, spline)
			for _ = 1, N do
				for i = 1, SOLVE_SAMPLES do
					local t = (i - 1) / SOLVE_SAMPLES
					spline:SolvePosition(t)
				end
			end
		end,
	},
}
