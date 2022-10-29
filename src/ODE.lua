local Tools = require("src.Tools")
local MA = require("src.MatrixAlgebra")
local Interpolation = require("src.Interpolation")
local cheb = Interpolation.Chebyshev

local ODE = {}

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
                for k = 2, p + 1, 1 do
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
        result[i] = MA.liSparseMatrix.newNumeric(D[i], 10^-14)
    end

    return result
end

local _chebyshevFirstOrderSpectralMethod = function (fList, a, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)

    local D = cheb.derivativeMatrix(n - 1, chebyshevGrid)
    for i = 1, n, 1 do
        if i == 1 then
            D[i] = 1
        else
            D[i] = 0
        end
    end

    local vector = Tools.list.copy(fList)

    vector[1] = a

    return MA.liSparseMatrix.solve(D, vector)
end

local _chebyshevSecondtOrderSpectralMethod = function (coeffList, fList, a, b, chebyshevGridPoints)
    local n = #fList
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n - 1)
    a = a or 0
    b = b or 0

    local DList = _chebyshevSpectralDifferentionMatrix(n - 1, 2, chebyshevGrid)

    local D = MA.liSparseMatrix.scale(DList[3], coeffList[1]) + MA.liSparseMatrix.scale(DList[2], coeffList[2]) + MA.liSparseMatrix.scale(DList[1], coeffList[3])

    for i = 1, n, 1 do
        if i == 1 then
            D[i] = 1
            D[n * (n-1) + i] = 0
        elseif i == n then
            D[i] = 0
            D[n * (n-1) + i] = 1
        else
            D[i] = 0
            D[n * (n-1) + i] = 0
        end
    end

    local vector = Tools.list.copy(fList)

    vector[1] = a
    vector[n] = b

    return MA.liSparseMatrix.solve(D, vector)
end

ODE.SpectralMethods = {}

ODE.SpectralMethods.chebyshevDerivativeMatrix = function (n, chebyshevGridPoints)
    return cheb.derivativeMatrix(n, chebyshevGridPoints)
end

ODE.SpectralMethods.chebyshevDerivativeMatrices = function (n, p, chebyshevGridPoints)
    return _chebyshevSpectralDifferentionMatrix(n, p, chebyshevGridPoints)
end

ODE.SpectralMethods.chebyshevFirstOrder = function (f, a, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local fList = Tools.list.map(f, chebyshevGrid)
    return _chebyshevFirstOrderSpectralMethod(fList, a, chebyshevGrid)
end

ODE.SpectralMethods.chebyshevSecondOrder = function (coeffList, f, a, b, n, chebyshevGridPoints)
    local chebyshevGrid = chebyshevGridPoints or cheb.grid(n)
    local fList = Tools.list.map(f, chebyshevGrid)
    return _chebyshevSecondtOrderSpectralMethod(coeffList, fList, a, b, chebyshevGrid)
end

return ODE