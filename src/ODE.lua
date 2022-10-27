local Tools = require("src.Tools")
local MA = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local cheb = Interpolation.Chebyshev

local ODE = {}

local _chebyshevSpectralDifferentionMatrix = function (n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
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

local _chebyshevFirstOrderSpectralMethod = function (fList, a, b, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)

    local D = _chebyshevSpectralDifferentionMatrix(n - 1, chebyshevGrid)
    for i = 1, n, 1 do
        if i == 1 and i == n then
            D[i] = 1
            D[n] = 1
            break
        elseif i == 1 then
            D[1] = 1
            D[n * (n - 1) + 1] = 0
        elseif i == n then
            D[n * n] = 1
            D[n] = 0
        else
            D[i] = 0
            D[n * (n - 1) + i] = 0
        end
    end

    local vector = Tools.list.copy(fList)

    vector[1], vector[n] = a, b

    return MA.liSparseMatrix.solve(D, vector)
end

ODE.SpectralMethods = {}

ODE.SpectralMethods.chebyshevDerivativeMatrix = function (n, chebyshevGridPoints)
    return _chebyshevSpectralDifferentionMatrix(n, chebyshevGridPoints)
end

ODE.SpectralMethods.chebyshevFirstOrder = function (f, a, b, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local fList = Tools.list.map(f, chebyshevGrid)
    return _chebyshevFirstOrderSpectralMethod(fList, a, b, chebyshevGrid)
end

return ODE