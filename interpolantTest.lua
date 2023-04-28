local Tools = require("src.Tools")
local Scalars = require("src.Scalars")
local MatrixAlgebra = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local SymbolicAlgebra = require("src.SymbolicAlgebra")
local ODE = require("src.ODE")
local NM = require("src.NumericalMethods")

print("\nObject Tests\n")
local tic = os.clock()
local test = Interpolation.ChebyshevInterpolant:new(function (x) return math.exp(x) end, 0, 1, 10)
local testPrime = test:derivative()

print("Solve test(x) = 2    ", Tools.list.tostring(test:solve(2)))
print("Solve test'(x) = 2   ", Tools.list.tostring(testPrime:solve(2)))
print("Solve exp(x) = 2     ", Tools.solve.regulaFalsi(math.exp, 2, 0, 1))
print("exp(0.69314718)      ", math.exp(0.69314718))
print("test(0.69314718)     ", test:evaluate(0.69314718))
local testInverse = test:inverse()
print("Evaluate test^(-1)(2)", testInverse:evaluate(2))
print("exp(0.405465)        ", math.exp(0.405465))
print("Compute test max     ", test:max())
print("exp(1)               ", math.exp(1))
local toc = os.clock()
print("Total Time Taken:    ", toc - tic)

print("\nMonotone Tests\n")
tic = os.clock()
test.solveMethod = "Monotone"
testPrime = test:derivative()

print("Solve test(x) = 2    ", Tools.list.tostring(test:solve(2)))
print("Solve test'(x) = 2   ", Tools.list.tostring(testPrime:solve(2)))
print("Solve exp(x) = 2     ", Tools.solve.regulaFalsi(math.exp, 2, 0, 1))
print("exp(0.69314718)      ", math.exp(0.69314718))
print("test(0.69314718)     ", test:evaluate(0.69314718))
testInverse = test:inverse()
print("Evaluate test^(-1)(2)", testInverse:evaluate(2))
print("exp(0.405465)        ", math.exp(0.405465))
print("Compute test max     ", test:max())
print("exp(1)               ", math.exp(1))
toc = os.clock()
print("Total Time Taken:    ", toc - tic)

print("\nNumerical Method Interface Test\n")

tic = os.clock()
print("Solve test(x) = 2 with monotone method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
test.solveMethod = "RegulaFalsi"
tic = os.clock()
print("Solve test(x) = 2 with built in method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
tic = os.clock()
print("Solve test(x) = 2 with Newton's method:", NM.solvers.newtonsMethod(test.evaluationFunction, test:derivative().evaluationFunction, 2))
print("Time taken:                            ", os.clock() - tic, "\n")

test = Interpolation.ChebyshevInterpolant:new(function (x) return math.exp(x) end, 0, 1, 20)
test.solveMethod = "Monotone"

tic = os.clock()
print("Solve test(x) = 2 with monotone method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
test.solveMethod = "RegulaFalsi"
tic = os.clock()
print("Solve test(x) = 2 with built in method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
tic = os.clock()
print("Solve test(x) = 2 with Newton's method:", NM.solvers.newtonsMethod(test.evaluationFunction, test:derivative().evaluationFunction, 2))
print("Time taken:                            ", os.clock() - tic, "\n")

test = Interpolation.ChebyshevInterpolant:new(function (x) return math.exp(x) end, 0, 1, 40)
test.solveMethod = "Monotone"

tic = os.clock()
print("Solve test(x) = 2 with monotone method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
test.solveMethod = "RegulaFalsi"
tic = os.clock()
print("Solve test(x) = 2 with built in method:", Tools.list.tostring(test:solve(2)))
print("Time taken:                            ", os.clock() - tic)
tic = os.clock()
print("Solve test(x) = 2 with Newton's method:", NM.solvers.newtonsMethod(test.evaluationFunction, test:derivative().evaluationFunction, 2))
print("Time taken:                            ", os.clock() - tic, "\n")

print("Integration Test\n")
print("True value of integral of exp from 0 to 1:         ", math.exp(1) - 1)
print("Five point Gaussian quadrature for exp from 0 to 1:", NM.integration.fivePointGaussianQuadrature(function (x) return math.exp(x) end, 0, 1), "\n")

print("Full Worked Example of Arclength Parameterization\n")
tic = os.clock()
local n = 100
local grid = Interpolation.Chebyshev.grid(n)
local linearRescalingFunction = Interpolation.Chebyshev.linearRescalingFunction(0, 1)
local shiftedGrid = {}
for i = 1, #grid, 1 do
    shiftedGrid[i] = linearRescalingFunction(grid[i])
end
print("Chebyshev Grid:        ", Tools.list.tostring(grid))
print("Shifted Chebyshev Grid:", Tools.list.tostring(shiftedGrid))
-- We will assume that |gamma'(x)| is (exp(x) - 1) / (exp(1) - 2) and thus gamma is length 1
local gridValues = {0}
for i = 1, #shiftedGrid - 1, 1 do
    gridValues[i + 1] = gridValues[i] + NM.integration.fivePointGaussianQuadrature(function (x) return (math.exp(x) - 1) / (math.exp(1) - 2) end, shiftedGrid[i], shiftedGrid[i + 1])
end
print("Arc Length Function Grid Values:", Tools.list.tostring(gridValues))
local interpolant = Interpolation.ChebyshevInterpolant:new(gridValues, 0, 1, n)
interpolant.solveMethod = "Monotone"
interpolant = interpolant:inverse()
print("Time Taken:                             ", os.clock() - tic)
print("Arc Length Parameterization Grid Values:", Tools.list.tostring(interpolant.gridValues))
print("This should be close to 0.5:            ", NM.integration.fivePointGaussianQuadrature(function (x) return (math.exp(x) - 1) / (math.exp(1) - 2) end, 0, interpolant:evaluate(0.5)))
tic = os.clock()
for i = 1, 1000000, 1 do
    interpolant:evaluate(0.5)
end
print("Time Taken for 1,000,000 Evaluatations:", os.clock() - tic)