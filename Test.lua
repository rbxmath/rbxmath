local Matrices         = require("src/Matrices")
local Matrix           = Matrices.Matrix
local ComplexMatrix    = Matrices.ComplexMatrix
local SparseMatrix     = Matrices.SparseMatrix
local FFT              = require("src/FastFourierTransform")
local Interpolation    = require("src/Interpolation")
local Tools            = require("src/Tools")
local ODE              = require("src/ODE")
local NumericalMethods = require("src/NumericalMethods")

print("Matrix Tests:\n")
local zero = Matrix:new({{0, 0}, {0, 0}})
local zer2 = Matrix:new({{0, 0}, {0, 0}})
local ones = Matrix:new({{1, 1}, {1, 1}})
local twos = Matrix:new({{2, 2}, {2, 2}})
local iden = Matrix:new({{1, 0}, {0, 1}})
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
    idet = iden * iden
end
print("  Matrix time for 1000\n  2x2 multiplications:      ", os.clock() - tic, "\n")
print("Complex Matrix Tests:\n")
zero = ComplexMatrix:new({{0, 0}, {0, 0}}, true)
zer2 = ComplexMatrix:new({{0, 0}, {0, 0}}, true)
ones = ComplexMatrix:new({{1, 1}, {1, 1}}, true)
twos = ComplexMatrix:new({{2, 2}, {2, 2}}, true)
iden = ComplexMatrix:new({{1, 0}, {0, 1}}, true)
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
    idet = iden * iden
end
print("  Matrix time for ,000\n  2x2 multiplications:      ", os.clock() - tic, "\n")
print("Sparse Matrix Tests:\n")
zero = SparseMatrix:new({{0, 0}, {0, 0}})
zer2 = SparseMatrix:new({{0, 0}, {0, 0}})
ones = SparseMatrix:new({{1, 1}, {1, 1}})
twos = SparseMatrix:new({{2, 2}, {2, 2}})
iden = SparseMatrix:new({{1, 0}, {0, 1}})
local idet = iden * iden
local onet = ones * ones
local onep = ones + ones
print("Matrix equality test:       ", zero == zer2)
print("Matrix multiplication test: ", twos == onet)
print("Matrix addition test:       ", twos == onep)
print("Matrix identity is identity:", iden == idet)
print("Matrix time test start:     ")
local tic = os.clock()
for i = 1, 1000 do
    idet = iden * iden
end
print("  Matrix time for 1000\n  2x2 multiplications:      ", os.clock() - tic, "\n")

print("ODE Test:\n")
local odes = ODE.SpectralMethods.chebyshevFirstOrder(function (x) return math.exp(x) end, 1 / math.exp(1), 20)
print("First order method test:    ", math.abs(odes[#odes] - math.exp(1)) < 10^(-13))

print("\nNumerical Methods Test:\n")

local expF = function (x) return math.exp(x) end
local nexF = function (x) return math.exp(-x) end
local exqu = function (x) return 1 / (x^4 + 1) end
local test = math.abs(NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1, 10^(-15)) - math.exp(1) + 1) < 10^(-15)
local adat = math.abs(NumericalMethods.integration.fivePointGaussianQuadrature(expF, 0, 1) - math.exp(1) + 1) > math.abs(NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1) - math.exp(1) + 1)
local alia = NumericalMethods.integration.integrate(expF, 0, 1) - NumericalMethods.integration.adaptiveQuadrature(expF, 0, 1) == 0
local lagu = math.abs(NumericalMethods.integration.fivePointLaguerre(nexF, 0) - 1) < 10^(-13)
local alag = math.abs(NumericalMethods.integration.adaptiveLaguerre(nexF, 0) - 1) < math.abs(NumericalMethods.integration.fivePointLaguerre(nexF, 0) - 1)
local exqQ = math.abs(NumericalMethods.integration.adaptiveLaguerre(exqu, 0) - math.pi / (2 * math.sqrt(2))) < 10^(-11)
local trap = math.abs(NumericalMethods.integration.adaptiveTrapezoid(expF, 0, 1) - math.exp(1) + 1) < 10^(-9)
print(math.abs(NumericalMethods.integration.adaptiveTrapezoid(expF, 0, 1) - math.exp(1) + 1))
print("Adap. Quad. Error < 10^-15: ", test)
print("Adapt. Better Than 5PGQ:    ", adat)
print("Integrate Alias for A. Quad:", alia)
print("Lagu. Quad. Error < 10^-13: ", lagu)
print("Adap. Lag. Better Than Lag.:", alag)
print("Adap. Lag. Error < 10^-11:  ", exqQ)
print("Adap. Trap. Error < 10^-9: ", trap)
tic = os.clock()
for i = 1, 1000 do
    test = NumericalMethods.integration.integrate(expF, 0, 1)
end
local toc = os.clock()
print("Time taken for 1000 integrations:", toc - tic)
print("                              HZ:", math.floor(1000 / (toc - tic)))