local Tools = require("src.Tools")

local Interpolation = {}

local _chebyshevGrid = function (n)
    local result = {}

    for i = 0, n, 1 do
        result[n + 1 - i] = math.cos(i * math.pi / n)
    end

    return result
end

local _barycentricInterpolationInChebyshevPointsAtAPoint = function (fList, x, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)

    local xDiffList = {}
    local diff
    for i, v in ipairs(chebyshevGrid) do
        diff = x - v
        if diff == 0 then
            return fList[i]
        end
        xDiffList[i] = diff
    end

    local numerator = 0.5 * fList[1] / (xDiffList[1])
    for i = 2, n-1, 1 do
        numerator = numerator + (-1)^(i - 1) * fList[i] / (xDiffList[i])
    end
    numerator = numerator + 0.5 * (-1)^(n - 1) * fList[n] / (xDiffList[n])

    local denominator = 0.5 / (xDiffList[1])
    for i = 2, n-1, 1 do
        denominator = denominator + (-1)^(i - 1) / (xDiffList[i])
    end
    denominator = denominator + 0.5 * (-1)^(n - 1) / (xDiffList[n])

    return numerator / denominator
end

local _linearRescalingFunction = function (a, b)
    return function (x)
        return (b - a) / 2 * (x - 1) + b
    end
end

Interpolation.Chebyshev = {}

Interpolation.Chebyshev.grid = function (n)
    return _chebyshevGrid(n)
end

Interpolation.Chebyshev.evaluate = function (f, x, n)
    local chebyshevGrid = _chebyshevGrid(n)
    local fList = {}
    for i = 1, n+1, 1 do
        fList[i] = f(chebyshevGrid[i])
    end
    return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x)
end

Interpolation.Chebyshev.rescaleAndEvaluate = function (f, a, b, x, n)
    local chebyshevGrid = _chebyshevGrid(n - 1)
    local rescalingFunction = _linearRescalingFunction(a, b)
    local fList = {}
    for i = 1, n, 1 do
        fList[i] = f(rescalingFunction(chebyshevGrid[i]))
    end
    return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x)
end

Interpolation.Chebyshev.evaluateOnData = function (fList, x)
    return _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x)
end

Interpolation.Chebyshev.evaluateAtPointList = function (f, xList, n)
    local chebyshevGrid = _chebyshevGrid(n - 1)
    local fList = {}
    for i = 1, n, 1 do
        fList[i] = f(chebyshevGrid[i])
    end

    local pList = {}
    for i = 1, #xList, 1 do
        pList[i] = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, xList[i], chebyshevGrid)
    end
    return pList
end

Interpolation.Chebyshev.evaluateAtPointListOnData = function (fList, xList)
    local chebyshevGrid = _chebyshevGrid(#fList - 1)
    local pList = {}
    for i = 1, #xList, 1 do
        pList[i] = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, xList[i], chebyshevGrid)
    end
    return pList
end

Interpolation.Chebyshev.benchmark = function (f, xList, deltaMax, maxIters)
    local delta = 0
    local iters = 0
    local mark, realAnswer, approximateAnswer
    local errorVector, timeVector = {}, {}

    while delta < deltaMax and iters < maxIters do
        iters = iters + 1
        realAnswer = Tools.list.map(f, xList)
        mark = os.clock()
        approximateAnswer = Interpolation.Chebyshev.evaluateAtPointList(f, xList, iters + 1)
        delta = os.clock() - mark
        errorVector[iters] = Tools.list.error(realAnswer, approximateAnswer)
        timeVector[iters] = delta
    end

    return errorVector, timeVector
end

return Interpolation