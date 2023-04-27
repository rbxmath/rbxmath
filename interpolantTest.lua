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
print("Five point Gaussian quadrature for exp from 0 to 1:", NM.integration.fivePointGaussianQuadrature(function (x) return math.exp(x) end, 0, 1))