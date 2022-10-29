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

local exponential = function (x)
    return math.exp(x)
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
local rinv = MatrixAlgebra.liSparseMatrix.inverse(rand)
local oneC = MatrixAlgebra.liSparseMatrix.new({{1}, {1}})
local solv = MatrixAlgebra.liSparseMatrix.solve(rand, {1, 1})
local lurd = MatrixAlgebra.liSparseMatrix.lu(rand)
local l, u = lurd[1], lurd[2]
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
local spring  = ODE.SpectralMethods.chebyshevSecondOrder({1, 8, 1}, function (x) return 0 end, 1, 0, 64)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
spring  = ODE.SpectralMethods.chebyshevSecondOrder({1, 32, 1}, function (x) return 0 end, 1, 0, 64)
mark2 = os.clock()
print(Tools.list.tostring(spring))
print("Time: " .. tostring(mark2 - mark1))
print()
mark1 = os.clock()
local laplace  = ODE.SpectralMethods.chebyshevSecondOrder({-1, 0, 0}, function (x) return 0 end, 1, 0, 64)
mark2 = os.clock()
print(Tools.list.tostring(laplace))
print("Time: " .. tostring(mark2 - mark1))
print()