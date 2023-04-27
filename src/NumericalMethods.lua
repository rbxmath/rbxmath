local Tools = require("src.Tools")
local MA = require("src.MatrixAlgebra")

local NumericalMethods = {}

NumericalMethods.solvers = {}

function NumericalMethods.solvers.regulaFalsi (f, t, a, b, tol)
    tolerance = tol or 10^(-13)
    local leftValue, rightValue, middleValue = f(a) - t, f(b) - t, 0
    local left, right, middle = a, b, (a + b) / 2
    if leftValue * rightValue > 0 then
        return nil
    elseif math.abs(leftValue) < tolerance then
        return left
    elseif math.abs(rightValue) < tolerance then
        return right
    end

    while math.abs(leftValue - rightValue) >= tolerance and math.abs(left - right) >= tolerance do
        middle = (right * leftValue - left * rightValue)/(leftValue - rightValue)
        middleValue = f(middle) - t
        if math.abs(middleValue) < tolerance or math.abs(left - middle) / math.abs(leftValue - rightValue) < tolerance or math.abs(right - middle) / math.abs(leftValue - rightValue) < tolerance then
            return middle
        elseif math.abs(leftValue) < tolerance then
            return left
        elseif math.abs(rightValue) < tolerance then
            return right
        elseif leftValue * middleValue > 0 then
            leftValue = middleValue
            left = middle
        else
            rightValue = middleValue
            right = middle
        end
    end
    
    return middle
end

function NumericalMethods.solvers.newtonsMethod (f, fPrime, t, x, tol)
    tol = tol or 10^(-13)
    x = x or 0
    local value = f(x) - t
    local valuePrime = fPrime(x)
    while math.abs(value) >= tol do
        x = x - value / valuePrime
        value = f(x) - t
        valuePrime = fPrime(x)
    end
    return x
end

NumericalMethods.integration = {}

NumericalMethods.integration.fivePointGaussianGrid = {-1 * math.sqrt(5 + 2 * math.sqrt(10 / 7))/3, -1 * math.sqrt(5 - 2 * math.sqrt(10 / 7))/3, 0, math.sqrt(5 - 2 * math.sqrt(10 / 7))/3, -1 * math.sqrt(5 + 2 * math.sqrt(10 / 7))/3}
NumericalMethods.integration.fivePointGaussianWeights = {(322 - 13 * math.sqrt(70)) / 900, (322 + 13 * math.sqrt(70)) / 900, 128 / 225, (322 + 13 * math.sqrt(70)) / 900, (322 - 13 * math.sqrt(70)) / 900}

function NumericalMethods.integration.fivePointGaussianQuadrature (f, a, b)
    local sum = 0
    for i = 1, 5, 1 do
        local x = NumericalMethods.integration.fivePointGaussianGrid[i]
        local w = NumericalMethods.integration.fivePointGaussianWeights[i]
        sum = sum + w * f(((b - a) * x + b + a) / 2)
    end
    return (b - a) * sum / 2
end

return NumericalMethods