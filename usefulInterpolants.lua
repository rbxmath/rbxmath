local Tools = require("src.Tools")
local Scalars = require("src.Scalars")
local MatrixAlgebra = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local SymbolicAlgebra = require("src.SymbolicAlgebra")
local ODE = require("src.ODE")

local thousandPointGrid = {}
for i = 0, 9999, 1 do
    thousandPointGrid[i + 1] = -1 + 2 * i / 9999
end

local scaledCosine = function (x)
    return math.cos(math.pi * x)
end

local scaledSine = function (x)
    return math.sin(math.pi * x / 4 + math.pi / 4)
end

local mark1 = os.clock()
local degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(scaledCosine, thousandPointGrid, 21)
local mark2 = os.clock()

print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(scaledCosine, thousandPointGrid))))

mark1 = os.clock()
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointList(scaledSine, thousandPointGrid, 13)
mark2 = os.clock()

print("Time Taken: " .. tostring(mark2 - mark1))
print("Error: " .. tostring(Tools.list.error(degreeOnehundredInterpolationEvaluatedAtAThousandPoints, Tools.list.map(scaledSine, thousandPointGrid))))

print(Tools.list.tostring(Interpolation.Chebyshev.grid(12)))

local FastCosine = function (x)
    x = x % math.pi
    if math.pi/2 > x then
        return Interpolation.Chebyshev.evaluateOnData({0.0, 0.026758599128759, 0.10502933764983, 0.2280143241917, 0.38268343236509, 0.54979780960971, 0.70710678118655, 0.83529777238322, 0.92387953251129, 0.97365777764233, 0.99446912382076, 0.99964192457733, 1.0}, 2 * x / math.pi - 1, {-1.0, -0.96592582628907, -0.86602540378444, -0.70710678118655, -0.5, -0.25881904510252, 0, 0.25881904510252, 0.5, 0.70710678118655, 0.86602540378444, 0.96592582628907, 1.0})
    else
        return -1 * Interpolation.Chebyshev.evaluateOnData({0.0, 0.026758599128759, 0.10502933764983, 0.2280143241917, 0.38268343236509, 0.54979780960971, 0.70710678118655, 0.83529777238322, 0.92387953251129, 0.97365777764233, 0.99446912382076, 0.99964192457733, 1.0}, 2 * x / math.pi - 1, {-1.0, -0.96592582628907, -0.86602540378444, -0.70710678118655, -0.5, -0.25881904510252, 0, 0.25881904510252, 0.5, 0.70710678118655, 0.86602540378444, 0.96592582628907, 1.0})
    end
end

print(FastCosine(0) .. " " .. FastCosine(1))

mark1 = os.clock()
local x = 0
for i = 1, 10000, 1 do
    x = math.cos(math.random())
end
mark2 = os.clock()
print("Time Taken: " .. tostring(mark2 - mark1))
mark1 = os.clock()
x = 0
for i = 1, 10000, 1 do
    x = FastCosine(math.random())
end
mark2 = os.clock()
print("Time Taken: " .. tostring(mark2 - mark1))
mark1 = os.clock()
x = 0
for i = 1, 10000, 1 do
    x = math.asin(math.random())
end
mark2 = os.clock()
print("Time Taken: " .. tostring(mark2 - mark1))

local flist = {
    0.75487766624669276005,
    0.75014854970538310667,
    0.73563040222459897968,
    0.71025691749133585358,
    0.67198876270634578665,
    0.61739906059291265014,
    0.54130936396595314605,
    0.43789044527593884368,
    0.30632001947410781234,
    0.15634106108309008959,
    0,
    -0.15634106108309008959,
    -0.30632001947410781234,
    -0.43789044527593884368,
    -0.54130936396595314605,
    -0.61739906059291265014,
    -0.67198876270634578665,
    -0.71025691749133585358,
    -0.73563040222459897968,
    -0.75014854970538310667,
    -0.75487766624669276005
}
degreeOnehundredInterpolationEvaluatedAtAThousandPoints = Interpolation.Chebyshev.evaluateAtPointListOnData(flist, thousandPointGrid)
print(Tools.list.tostring(degreeOnehundredInterpolationEvaluatedAtAThousandPoints))