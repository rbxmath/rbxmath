<<<<<<< HEAD
local Matrices = require("../src/Matrices")
local Matrix = Matrices.Matrix
local ComplexMatrix = Matrices.ComplexMatrix
local SparseMatrix = Matrices.SparseMatrix
local FFT = require("../src/FastFourierTransform")
local Interpolation = require("../src/Interpolation")
local Tools = require("../src/Tools")
local ODE = require("../src/ODE")
local NumericalMethods = require("../src/NumericalMethods")
local Sets = require("../src/Sets")
local Set = Sets.Set
local Grid = require("../src/Grid")

print("Matrix Tests:\n")
local zero = Matrix:new({ { 0, 0 }, { 0, 0 } })
local zer2 = Matrix:new({ { 0, 0 }, { 0, 0 } })
local ones = Matrix:new({ { 1, 1 }, { 1, 1 } })
local twos = Matrix:new({ { 2, 2 }, { 2, 2 } })
local iden = Matrix:new({ { 1, 0 }, { 0, 1 } })
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
	idet = iden * iden
end
print("  Matrix time for 1000\n  2x2 multiplications:      ", os.clock() - tic, "\n")
print("Complex Matrix Tests:\n")
zero = ComplexMatrix:new({ { 0, 0 }, { 0, 0 } }, true)
zer2 = ComplexMatrix:new({ { 0, 0 }, { 0, 0 } }, true)
ones = ComplexMatrix:new({ { 1, 1 }, { 1, 1 } }, true)
twos = ComplexMatrix:new({ { 2, 2 }, { 2, 2 } }, true)
iden = ComplexMatrix:new({ { 1, 0 }, { 0, 1 } }, true)
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
	idet = iden * iden
end
print("  Matrix time for ,000\n  2x2 multiplications:      ", os.clock() - tic, "\n")
print("Sparse Matrix Tests:\n")
zero = SparseMatrix:new({ { 0, 0 }, { 0, 0 } })
zer2 = SparseMatrix:new({ { 0, 0 }, { 0, 0 } })
ones = SparseMatrix:new({ { 1, 1 }, { 1, 1 } })
twos = SparseMatrix:new({ { 2, 2 }, { 2, 2 } })
iden = SparseMatrix:new({ { 1, 0 }, { 0, 1 } })
print(Tools.list.deeptostring(ones:arnoldiProcess(2)))
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
	idet = iden * iden
=======
local Tools = require("src.Tools")
local Scalars = require("src.Scalars")
local MatrixAlgebra = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local SymbolicAlgebra = require("src.SymbolicAlgebra")
local ODE = require("src.ODE")

print("Testing Tools")
print()

local testList1 = {1, 2, 3, 4, 5, 6}
print(Tools.list.binarySearch(4, testList1))
print(Tools.list.binarySearch(4.5, testList1))

local phi = (1 + math.sqrt(5)) / 2
print(phi)
local continuedFractionArrayOfPhi = Tools.continuedFractions.compute(phi, 10)
print(Tools.list.tostring(continuedFractionArrayOfPhi))
local phiApproximation = Tools.continuedFractions.evaluate(continuedFractionArrayOfPhi)
print(phiApproximation)
local phiRational = Scalars.Rational.rationalFromContinuedFraction(continuedFractionArrayOfPhi)
print(phiRational)
local phiRationalToReal = Scalars.Rational.toReal(phiRational)
print(phiRationalToReal)
print()

print("Testing Interpolation")
print()

local stepFunction = function (x)
    if x < 0 then
        return 0
    else
        return 1
    end
end
local thousandPointGrid = {}
for i = 0, 999, 1 do
    thousandPointGrid[i + 1] = -1 + 2 * i / 999
end
local mark1 = os.clock()
local degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(stepFunction, thousandPointGrid, 100)
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))
local mark2 = os.clock()
print("Approximating step function with degree 99 at 1000 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(stepFunction, thousandPointGrid))))
print()

local polynomial = function (x)
    return 1 + 10 * x + 100 * x^3 + 1000 * x^4
>>>>>>> main
end
mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(polynomial, thousandPointGrid, 100)
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))
print(Tools.list.tostring(Interpolation.Chebyshev.solve(polynomial, 0, 1000)))
mark2 = os.clock()
print("Approximating step function with degree 99 at 1000 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(stepFunction, thousandPointGrid))))
print()

<<<<<<< HEAD
print("ODE Test:\n")
local odes = ODE.SpectralMethods.chebyshevFirstOrder(function(x)
	return math.exp(x)
end, 1 / math.exp(1), 20)
print("First order method test:    ", math.abs(odes[#odes] - math.exp(1)) < 10 ^ -13)

print("\nNumerical Methods Test:\n")

local expF = function(x)
	return math.exp(x)
end
local nexF = function(x)
	return math.exp(-x)
end
local exqu = function(x)
	return 1 / (x ^ 4 + 1)
end
local test = math.abs(NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1, 10 ^ -15) - math.exp(1) + 1)
	< 10 ^ -15
local adat = math.abs(NumericalMethods.integration.fivePointGaussianQuadrature(expF, 0, 1) - math.exp(1) + 1)
	> math.abs(NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1) - math.exp(1) + 1)
local alia = NumericalMethods.integration.integrate(expF, 0, 1)
		- NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1)
	== 0
local lagu = math.abs(NumericalMethods.integration.fivePointLaguerre(nexF, 0) - 1) < 10 ^ -13
local alag = math.abs(NumericalMethods.integration.adaptiveLaguerre(nexF, 0) - 1)
	< math.abs(NumericalMethods.integration.fivePointLaguerre(nexF, 0) - 1)
local exqQ = math.abs(NumericalMethods.integration.adaptiveLaguerre(exqu, 0) - math.pi / (2 * math.sqrt(2))) < 10 ^ -11
local trap = math.abs(NumericalMethods.integration.adaptiveTrapezoid(expF, 0, 1) - math.exp(1) + 1) < 10 ^ -9
print("Adap. Quad. Error < 10^-15: ", test)
print("Adapt. Better Than 5PGQ:    ", adat)
print("Integrate Alias for A. Quad:", alia)
print("Lagu. Quad. Error < 10^-13: ", lagu)
print("Adap. Lag. Better Than Lag.:", alag)
print("Adap. Lag. Error < 10^-11:  ", exqQ)
print("Adap. Trap. Error < 10^-9: ", trap)
tic = os.clock()
for i = 1, 1000 do
	test = NumericalMethods.integration.integrate(expF, 0, 1)
=======
local exponential = function (x)
    return math.exp(x)
>>>>>>> main
end
mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(exponential, thousandPointGrid, 100)
mark2 = os.clock()
print("Approximating exponential function with degree 99 at 1000 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(exponential, thousandPointGrid))))
print()
local hundredPointGrid = {}
for i = 0, 99, 1 do
    hundredPointGrid[i + 1] = -1 + 2 * i / 99
end
mark1 = os.clock()
local degreeOnehundredInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(exponential, hundredPointGrid, 100)
mark2 = os.clock()
print("Approximating exponential function with degree 99 at 100 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAHundredPoints, Tools.list.map(exponential, hundredPointGrid))))
print()
mark1 = os.clock()
local degreeTenInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(exponential, hundredPointGrid, 10)
mark2 = os.clock()
print("Approximating exponential function with degree 9 at 100 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeTenInterpolationEvaluatedAtAHundredPoints, Tools.list.map(exponential, hundredPointGrid))))
print()

local splineArcLength = function (x)
    return math.sqrt(
        (4 * x^2 - 2 * x + 1)^2 +
        (x^2 - 12)^2 +
        (x^2 + x + 1)^2
    )
end
mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(splineArcLength, thousandPointGrid, 100)
mark2 = os.clock()
print("Approximating splineArcLength type function with degree 99 at 1000 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(splineArcLength, thousandPointGrid))))
print()
hundredPointGrid = {}
for i = 0, 99, 1 do
    hundredPointGrid[i + 1] = -1 + 2 * i / 99
end
mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(splineArcLength, hundredPointGrid, 100)
mark2 = os.clock()
print("Approximating splineArcLength type function with degree 99 at 100 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAHundredPoints, Tools.list.map(splineArcLength, hundredPointGrid))))
print()
mark1 = os.clock()
degreeTenInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(splineArcLength, hundredPointGrid, 10)
mark2 = os.clock()
print("Approximating splineArcLength type function with degree 9 at 100 points")
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeTenInterpolationEvaluatedAtAHundredPoints, Tools.list.map(splineArcLength, hundredPointGrid))))
print()
local errorVec, timeVec = Interpolation.Chebyshev.benchmark(splineArcLength, hundredPointGrid, 1, 25)
print("Benchmarking algorithm on splineArcLength function")
print("Error: " .. Tools.list.tostring(errorVec))
print("Time: " .. Tools.list.tostring(timeVec))
print()

<<<<<<< HEAD
test = Set:new({ 1, 1, 3, 2 })
local tes2 = test:copy():addTo(4)
local tes3 = test:copy():addTo(3)
print("Sets are sorted:                 ", test.data[1] == 1 and test.data[2] == 2 and test.data[3] == 3)
print(
	"Adding to sets works:            ",
	tes2.data[1] == 1
		and tes2.data[2] == 2
		and tes2.data[3] == 3
		and tes2.data[4] == 4
		and tes2.cardinality == 4
		and tes3.data[1] == 1
		and tes3.data[2] == 2
		and tes3.data[3] == 3
		and tes3.cardinality == 3
)
=======
local spline = function (x)
    return math.abs(x)^3
end
errorVec, timeVec = Interpolation.Chebyshev.benchmark(spline, hundredPointGrid, 1, 100)
print("Benchmarking algorithm on cubic spline")
print("Error: " .. Tools.list.tostring(errorVec))
print("Time: " .. Tools.list.tostring(timeVec))
print()

local randomVectorOfLength100 = {}
for i = 1, 100, 1 do
    randomVectorOfLength100[i] = 1 - 2 * math.random()
end
mark1 = os.clock()
local rootList = Interpolation.Chebyshev.solveFromData(randomVectorOfLength100, 0)
mark2 = os.clock()
print("Computing roots through random data")
print(Tools.list.tostring(randomVectorOfLength100))
print("Roots: " .. Tools.list.tostring(rootList))
print("Time: " .. tostring(mark2 - mark1))
print()

randomVectorOfLength100 = {}
randomVectorOfLength100[1] = -2
for i = 2, 100, 1 do
    randomVectorOfLength100[i] = randomVectorOfLength100[i - 1] + math.random()
end
mark1 = os.clock()
rootList = Interpolation.Chebyshev.solveFromData(randomVectorOfLength100, 0)
mark2 = os.clock()
print("Computing roots through random increasing data")
print(Tools.list.tostring(randomVectorOfLength100))
print("Roots: " .. Tools.list.tostring(rootList))
print("Time: " .. tostring(mark2 - mark1))
print()

print("Testing MatrixAlgebra")
print()

local ones = MatrixAlgebra.liSparseMatrix.new({{1, 1}, {1, 1}})
local twos = ones * ones
local rand = MatrixAlgebra.liSparseMatrix.random(2, 2, -1, 1, 0)
print(tostring(rand))
local rinv = MatrixAlgebra.liSparseMatrix.inverse(rand)
local oneC = MatrixAlgebra.liSparseMatrix.new({{1}, {1}})
local solv = MatrixAlgebra.liSparseMatrix.solve(rand, {1, 1})
local lurd = MatrixAlgebra.liSparseMatrix.lu(rand)
local l, u = lurd[1], lurd[2]
local iden = MatrixAlgebra.liSparseMatrix.new({{1, 0}, {0, 1}})
print(tostring(ones))
print(tostring(MatrixAlgebra.liSparseMatrix.scale(ones, 3)))
print(tostring(twos))
print(tostring(rand))
print(tostring(rinv))
print()
print(tostring(rinv * rand))
print(tostring(rinv * oneC))
print()
print(tostring(u))
print(tostring(solv))
print(tostring(rand * solv))
print()
print(iden)
print(iden + iden)
print(MatrixAlgebra.liSparseMatrix.scale(iden, {3, 2}))
print()

print("Testing sparseMatrix")
print()

ones = MatrixAlgebra.sparseMatrix.new({{1, 1}, {1, 1}})
twos = ones * ones
rand = MatrixAlgebra.sparseMatrix.random(2, 2, -1, 1, 0)
rinv = MatrixAlgebra.sparseMatrix.inverse(rand)
oneC = MatrixAlgebra.sparseMatrix.new({{1}, {1}})
solv = MatrixAlgebra.sparseMatrix.solve(rand, {1, 1})
lurd = MatrixAlgebra.sparseMatrix.lu(rand)
l, u = lurd[1], lurd[2]
iden = MatrixAlgebra.sparseMatrix.new({{1, 0}, {0, 1}})
print(tostring(ones))
print(tostring(twos))
print(tostring(rand))
print(tostring(rinv))
print()
print(tostring(rinv * rand))
print(tostring(rinv * oneC))
print()
print(tostring(u))
print(tostring(solv))
print(tostring(rand * solv))
print()
print(iden)
print(iden + iden)
print()

print("Testing matrix")
print()

ones = MatrixAlgebra.matrix.new({{1, 1}, {1, 1}})
twos = ones * ones
rand = MatrixAlgebra.matrix.random(2, 2, -1, 1, 0)
rinv = MatrixAlgebra.matrix.inverse(rand)
oneC = MatrixAlgebra.matrix.new({{1}, {1}})
solv = MatrixAlgebra.matrix.solve(rand, {1, 1})
lurd = MatrixAlgebra.matrix.lu(rand)
l, u = lurd[1], lurd[2]
iden = MatrixAlgebra.matrix.new({{1, 0}, {0, 1}})
print(tostring(ones))
print(tostring(twos))
print(tostring(rand))
print(tostring(rinv))
print()
print(tostring(rinv * rand))
print(tostring(rinv * oneC))
print()
print(tostring(u))
print(tostring(Tools.list.tostring(solv)))
print(tostring(Tools.list.tostring(MatrixAlgebra.matrix.apply(rand, solv))))
print()
print(iden)
print(iden + iden)
print()

print("Testing ODE")
print()

local diff3 = ODE.SpectralMethods.chebyshevDerivativeMatrix(2)
local diff3deg3 = ODE.SpectralMethods.chebyshevDerivativeMatrices(2, 3)
mark1 = os.clock()
local sine = ODE.SpectralMethods.chebyshevFirstOrder(math.cos, math.sin(-1), 32)
mark2 = os.clock()
print(diff3)
print(Tools.list.tostring(diff3deg3))
print()
print(Tools.list.tostring(sine))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
local splineInterp  = ODE.SpectralMethods.chebyshevFirstOrder(spline, 0, 32)
mark2 = os.clock()
print(Tools.list.tostring(splineInterp))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
local speed = function (x)
    return math.sqrt((x + 1)^2 + 1)
end
local length  = ODE.SpectralMethods.chebyshevFirstOrder(speed, 0, 32)
mark2 = os.clock()
print(Tools.list.tostring(length))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
local spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 8, 1}, function (x) return 0 end, {-1, 1}, {1, 0}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 32, 1}, function (x) return 0 end, {-1, 1}, {1, 0}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
local laplace  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({-1, 0, 0}, function (x) return 0 end, {-1, 1}, {1, 0}, 32)
mark2 = os.clock()
print(Tools.list.tostring(laplace))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 1, 100}, function (x) return x end, {-1, 1}, {0, 0}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderInitialValueProblem({1, 32, 1}, function (x) return 0 end, {-1, 1}, {0, 1}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderInitialValueProblem({1, 0, 10}, function (x) return 0 end, {-1, 1}, {0, 0.5}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderInitialValueProblem({1, function (x) return x end, 10}, function (x) return 0 end, {-1, 1}, {0, 0.5}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 0, function (x) return -x end}, function (x) return 0 end, {-10, 10}, {1, 0}, 32)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 0, function (x) return -x end}, function (x) return 0 end, {-30, 30}, {1, 0}, 128)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem({1, 0, function (x) return -x end}, function (x) return 0 end, {-30, 30}, {1, 0}, 256)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
>>>>>>> main
