local Tools = require("src.Tools")
local MA = require("src.MatrixAlgebra")

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

local _barycentricInterpolationInChebyshevPointsAtPointList = function (fList, xList, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
    local result = {}
    for j = 1, #xList, 1 do
        local x = xList[j]
        local xDiffList = {}
        local diff
        local index = -1
        for i, v in ipairs(chebyshevGrid) do
            diff = x - v
            if diff == 0 then
                index = i
            end
            xDiffList[i] = diff
        end

        if index ~= -1 then
            result[j] = fList[index]
        else
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

            result[j] = numerator / denominator
        end
    end

    return result
end

local _barycentricInterpolationSolve = function (fList, t, tol, chebyshevGridPoints, gapTol)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n)
    local tolerance = tol or (10^-13)
    local gapTolerance = gapTol or (10^-13)

    local intervalList = {}
    local rootList = {}
    if math.abs(fList[1] - t) < tolerance then
        rootList[#rootList+1] = chebyshevGrid[1]
    end
    for i = 1, n - 1, 1 do
        if t <= fList[i + 1] - tolerance and t >= fList[i] + tolerance or t <= fList[i] - tolerance and t >= fList[i + 1] + tolerance then
            intervalList[#intervalList+1] = {chebyshevGrid[i], chebyshevGrid[i + 1]}
        end
        if math.abs(fList[i + 1] - t) < tolerance then
            rootList[#rootList+1] = chebyshevGrid[i + 1]
        end
    end

    while #intervalList ~= 0 do
        local newIntervalList = {}
        local x0, x1, x2, f0, f1, f2, temp
        for i, v in ipairs(intervalList) do
            x0, x1 = v[1], v[2]
            temp =  _barycentricInterpolationInChebyshevPointsAtPointList(fList, {x0, x1}, chebyshevGrid)
            f0, f1 = temp[1] - t, temp[2] - t
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            if x2 < x0 or x2 > x1 then
                x2 = (x0 + x1) / 2
            end
            f2 = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x2, chebyshevGrid) - t
            if math.abs(f2) < tolerance then
                rootList[#rootList+1] = x2
            else
                local gap = gapTolerance + (x1 - x0) + 1
                while math.abs(f2) >= tolerance and gap >= gapTolerance and gap ~= x1 - x0 do
                    gap = x1 - x0
                    if ((f0 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f0)) and ((f1 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f1)) then
                        newIntervalList[#newIntervalList+1] = {x2, x1}
                        x1 = x2
                        f1 = f2
                    elseif (f0 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f0) then
                        x1 = x2
                        f1 = f2
                    elseif (f1 < -tolerance and tolerance < f2) or (f2 < -tolerance and tolerance < f1) then
                        x0 = x2
                        f0 = f2
                    else
                        break
                    end
                    x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
                    if x2 < x0 or x2 > x1 then
                        x2 = (x0 + x1) / 2
                    end
                    f2 = _barycentricInterpolationInChebyshevPointsAtAPoint(fList, x2, chebyshevGrid) - t
                end
                if math.abs(f2) < tolerance or gap < gapTolerance then
                    rootList[#rootList+1] = x2
                end
            end
        end
        intervalList = newIntervalList
    end
    return rootList
end

local _chebyshevSpectralDifferentionMatrix = function (n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n)
    local d = {{}, {}}

    for i = 1, n+1, 1 do
        d[1][i] = 1
        d[2][i] = 0
    end

    for i = 1, n+1, 1 do
        for j = 1, n+1, 1 do
            if i ~= j then
                d[2][i] = d[2][i] + d[1][i] / (chebyshevGrid[i] - chebyshevGrid[j])
            end
        end
    end

    local D = {}

    for i = 1, n+1, 1 do
        D[i] = {}
    end

    local c = {}
    for i = 1, n+1, 1 do
        c[i] = 1
        for j = 1, n+1, 1 do
            if i ~= j then
                c[i] = c[i] * (chebyshevGrid[i] - chebyshevGrid[j])
            end
        end
    end

    for i = 1, n+1, 1 do
        for j = 1, n+1, 1 do
            if i == j then
                D[i][j] = d[2][i]
            else
                D[i][j] = c[i] / (c[j] * (chebyshevGrid[i] - chebyshevGrid[j]))
            end
        end
    end

    return MA.liSparseMatrix.new(D)
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

Interpolation.Chebyshev.evaluateAtPointList = function (f, xList, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
    local fList = {}
    for i = 1, n, 1 do
        fList[i] = f(chebyshevGrid[i])
    end

    return _barycentricInterpolationInChebyshevPointsAtPointList(fList, xList, chebyshevGrid)
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

Interpolation.Chebyshev.solve = function (f, t, n, method, tol, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
    local tolerance = tol or (10^-13)
    local fList = {}
    for i = 1, n, 1 do
        fList[i] = f(chebyshevGrid[i])
    end

    local methodToUse = method or "RegulaFalsi"

    if methodToUse == "RegulaFalsi" then
        return _barycentricInterpolationSolve(fList, t, tolerance, chebyshevGrid)
    end
end

Interpolation.Chebyshev.solveFromData = function (fList, t, method, tol, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or _chebyshevGrid(n - 1)
    local tolerance = tol or (10^-13)

    local methodToUse = method or "Secant"

    if methodToUse == "Secant" then
        return _barycentricInterpolationSolve(fList, t, tolerance, chebyshevGrid)
    end
end

Interpolation.Chebyshev.derivativeMatrix = function (n, chebyshexGridPoints)
    local chebyshevGrid = chebyshexGridPoints or _chebyshevGrid(n)
    return _chebyshevSpectralDifferentionMatrix(n, chebyshevGrid)
end

return Interpolation