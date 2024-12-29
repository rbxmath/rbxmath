--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require(script.Parent.Tools)
local Matrices = require(script.Parent.Matrices)
local Matrix = Matrices.Matrix
local Interpolation = require(script.Parent.Interpolation)
local cheb = Interpolation.Chebyshev

local ODE = {}

local _linearRescalingFunction = function (a, b)
    return function (x)
        return (b - a) / 2 * (x - 1) + b
    end
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
    vector = D:solve(vector)

    return vector
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
