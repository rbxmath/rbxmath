local Tools = require("src.Tools")
local Scalars = require("src.Scalars")
local MatrixAlgebra = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local SymbolicAlgebra = require("src.SymbolicAlgebra")

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

print(Tools.list.tostring(Interpolation.Chebyshev.grid(10)))

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
local mark2 = os.clock()
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(stepFunction, thousandPointGrid))))
print()

local exponential = function (x)
    return math.exp(x)
end
mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(exponential, thousandPointGrid, 100)
mark2 = os.clock()
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))
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
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAHundredPoints))
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAHundredPoints, Tools.list.map(exponential, hundredPointGrid))))
print()
mark1 = os.clock()
local degreeTenInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(exponential, hundredPointGrid, 10)
mark2 = os.clock()
print(Tools.list.tostring(degreeTenInterpolationEvaluatedAtAHundredPoints))
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
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))
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
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAHundredPoints))
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAHundredPoints, Tools.list.map(splineArcLength, hundredPointGrid))))
print()
mark1 = os.clock()
degreeTenInterpolationEvaluatedAtAHundredPoints = Interpolation.Chebyshev.evaluateAtPointList(splineArcLength, hundredPointGrid, 10)
mark2 = os.clock()
print(Tools.list.tostring(degreeTenInterpolationEvaluatedAtAHundredPoints))
print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeTenInterpolationEvaluatedAtAHundredPoints, Tools.list.map(splineArcLength, hundredPointGrid))))
print()
local errorVec, timeVec = Interpolation.Chebyshev.benchmark(splineArcLength, hundredPointGrid, 1, 100)
print("Error: " .. Tools.list.tostring(errorVec))
print("Time: " .. Tools.list.tostring(timeVec))