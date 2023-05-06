local Tools = require("src/Tools")
local Matrices = require("src/Matrices")
local Matrix = Matrices.Matrix
local Interpolation = require("src/Interpolation")
local cheb = Interpolation.Chebyshev

local ODE = {}

local _linearRescalingFunction = function (a, b)
    return function (x)
        return (b - a) / 2 * (x - 1) + b
    end
end

local _chebyshevSpectralDifferentionMatrix = function (n, p, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local d = {}

    for i = 1, p + 1, 1 do
        d[i] = {}
    end

    for i = 1, n+1, 1 do
        d[1][i] = 1
        for j = 2, p + 1, 1 do
            d[j][i] = 0
        end
    end

    for i = 1, n+1, 1 do
        for j = 1, n+1, 1 do
            if i ~= j then
                for k = p + 1, 2, -1 do
                    d[k][i] = d[k][i] + (k - 1) * d[k - 1][i] / (chebyshevGrid[i] - chebyshevGrid[j])
                end
            end
        end
    end

    local D = {}

    for i = 1, p+1, 1 do
        D[i] = {}
        for j = 1, n+1, 1 do
            D[i][j] = {}
            for k = 1, n+1, 1 do
                if j == k then
                    D[i][j][k] = d[i][j]
                elseif i == 1 then
                    D[i][j][k] = 0
                end
            end
        end
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
            for k = 2, p+1, 1 do
                if i ~= j then
                    D[k][i][j] = (k - 1) * (d[k - 1][i] - D[k - 1][i][j]) / (chebyshevGrid[i] - chebyshevGrid[j])
                end
            end
        end
    end

    for i = 1, n+1, 1 do
        for j = 1, n+1, 1 do
            for k = 1, p+1, 1 do
                if i ~= j then
                    D[k][i][j] = D[k][i][j] * c[i] / c[j]
                end
            end
        end
    end

    local result = {}
    for i = 1, p+1, 1 do
        result[i] = Matrix:new(D[i])
    end

    return result
end

local _chebyshevFirstOrderSpectralMethod = function (fList, a, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)

    local D = cheb.derivativeMatrix(n - 1, chebyshevGrid)
    for i = 1, n, 1 do
        if i == 1 then
            D[1][i] = 1
        else
            D[1][i] = 0
        end
    end

    local vector = Tools.list.copy(fList)

    vector[1] = a

    return D:solve(vector)
end

local _chebyshevSecondtOrderBoundaryValueProblemSpectralMethod = function (coeffList, fList, a, b, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)
    a = a or 0
    b = b or 0

    local DList = _chebyshevSpectralDifferentionMatrix(n - 1, 2, chebyshevGrid)
    local D = DList[3]:scale( coeffList[1]) + DList[2]:scale(coeffList[2]) + DList[1]:scale(coeffList[3])

    for i = 1, n, 1 do
        if i == 1 then
            D[1][1] = 1
            D[n][1] = 0
        elseif i == n then
            D[1][i] = 0
            D[n][i] = 1
        else
            D[1][i] = 0
            D[n][i] = 0
        end
    end
    local vector = Tools.list.copy(fList)

    vector[1] = a
    vector[n] = b

    return D:solve(vector)
end

local _chebyshevSecondtOrderInitialValueProblemSpectralMethod = function (coeffList, fList, a, b, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)
    a = a or 0
    b = b or 0

    local DList = _chebyshevSpectralDifferentionMatrix(n - 1, 2, chebyshevGrid)
    local D = DList[3]:scale( coeffList[1]) + DList[2]:scale(coeffList[2]) + DList[1]:scale(coeffList[3])

    for i = 1, n, 1 do
        if i == 1 then
            D[1][i] = 1
            D[2][i] = DList[2][1][i]
        else
            D[1][i] = 0
            D[2][i] = DList[2][1][i]
        end
    end
    local vector = Tools.list.copy(fList)

    vector[1] = a
    vector[2] = b

    return D:solve(vector)
end

ODE.SpectralMethods = {}

function ODE.SpectralMethods.chebyshevDerivativeMatrix (n, chebyshevGridPoints)
    return cheb.derivativeMatrix(n, chebyshevGridPoints)
end

function ODE.SpectralMethods.chebyshevDerivativeMatrices (n, p, chebyshevGridPoints)
    return _chebyshevSpectralDifferentionMatrix(n, p, chebyshevGridPoints)
end

function ODE.SpectralMethods.chebyshevFirstOrder (f, a, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local fList = Tools.list.map(f, chebyshevGrid)
    return _chebyshevFirstOrderSpectralMethod(fList, a, chebyshevGrid)
end

function ODE.SpectralMethods.chebyshevSecondOrderBoundaryValueProblem (coeffList, f, boundary, boundaryValues, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local y, x = boundary[1], boundary[2]
    for i, v in ipairs(coeffList) do
        if type(v) == "function" then
            coeffList[i] = Tools.list.map(function (t) return ((x - y) / 2)^(i - #coeffList) * t end, Tools.list.map(v, Tools.list.map(_linearRescalingFunction(y, x), chebyshevGrid)))
        elseif type(v) == "number" then
            coeffList[i] = ((x - y) / 2)^(i - #coeffList) * coeffList[i]
        end
    end
    local a, b = boundaryValues[1], boundaryValues[2]
    local fList = Tools.list.map(f, Tools.list.map(_linearRescalingFunction(y, x), chebyshevGrid))
    return _chebyshevSecondtOrderBoundaryValueProblemSpectralMethod(coeffList, fList, a, b, chebyshevGrid)
end

function ODE.SpectralMethods.chebyshevSecondOrderInitialValueProblem (coeffList, f, boundary, boundaryValues, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local y, x = boundary[1], boundary[2]
    for i, v in ipairs(coeffList) do
        if type(v) == "function" then
            coeffList[i] = Tools.list.map(function (t) return ((x - y) / 2)^(i - #coeffList) * t end, Tools.list.map(v, Tools.list.map(_linearRescalingFunction(y, x), chebyshevGrid)))
        elseif type(v) == "number" then
            coeffList[i] = ((x - y) / 2)^(i - #coeffList) * coeffList[i]
        end
    end
    local a, b = boundaryValues[1], (x - y) * boundaryValues[2] / 2
    local fList = Tools.list.map(f, Tools.list.map(_linearRescalingFunction(y, x), chebyshevGrid))
    return _chebyshevSecondtOrderInitialValueProblemSpectralMethod(coeffList, fList, a, b, chebyshevGrid)
end

return ODE