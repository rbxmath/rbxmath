local Tools = require("src.Tools")
local Scalars = require("src.Scalars")
local MatrixAlgebra = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local SymbolicAlgebra = require("src.SymbolicAlgebra")
local ODE = require("src.ODE")

local test = Interpolation.Chebyshev.makeInterpolant(function (x) return math.exp(x) end, 0, 1, 10)
local testPrime = test:derivative()
local testInverse = test:inverse()

print("Interpolant Tests")

print("Solve test(x) = 1    ", Tools.list.tostring(test:solve(1)))
print("Solve test'(x) = 1   ", Tools.list.tostring(testPrime:solve(1)))
print("Evaluate test^(-1)(1.5)", testInverse:evaluate(1.5))
print("exp(0.405465)        ", math.exp(0.405465))
print("Compute test max     ", test:max())
print("exp(1)               ", math.exp(1))