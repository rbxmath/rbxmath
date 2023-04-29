local Tools = require("src.Tools")
local Matrices = require("src.Matrices")
local Matrix = Matrices.Matrix
local MA = require("src.MatrixAlgebra")

print("Basic Tests:")

local ones = Matrix:new({{1, 1}, {1, 1}})
print(ones)
print(ones * ones)
local data = {}
n = 4
for i = 1, n, 1 do
    data[i] = {}
    for j = 1, n, 1 do
        data[i][j] = 1
    end
end
local Ones = Matrix:new(data)
data = {}
n = 16
for i = 1, n, 1 do
    data[i] = {}
    for j = 1, n, 1 do
        data[i][j] = 1
    end
end
local ONes = Matrix:new(data)
data = {}
n = 64
for i = 1, n, 1 do
    data[i] = {}
    for j = 1, n, 1 do
        data[i][j] = 1
    end
end
local ONEs = Matrix:new(data)
local MAONEs = MA.matrix.new(data)
data = {}
n = 256
for i = 1, n, 1 do
    data[i] = {}
    for j = 1, n, 1 do
        data[i][j] = 1
    end
end

print("Old Versus New Tests:")

local ONES = Matrix:new(data)
local MAONES = MA.matrix.new(data)
local tic = os.clock()
temp = Ones * Ones
time1 = (os.clock() - tic)
tic = os.clock()
temp = ONes * ONes
time2 = (os.clock() - tic)
tic = os.clock()
temp = ONEs * ONEs
time3 = (os.clock() - tic)
tic = os.clock()
temp = ONES * ONES
time4 = (os.clock() - tic)
print(time1, time2, time3, time4)
print(" ", time2 / time1, time3 / time2, time4 / time3)

print("Matrix Computation Tests:")

local iden = Matrix:identity(2)
print(iden)
print(iden:inverse())
print(Tools.list.tostring(iden:solve({1, 1})))
print(iden:transpose())
print(ones:determinant())
print(iden:determinant())
print(iden:permuted({2, 1}):determinant())

print("Eigenvalue Computation Tests:")

print(Ones:hessenbergForm())
tic = os.clock()
print(Tools.list.tostring(Ones:francisOne()))
print("2 x 2 Ones:", os.clock() - tic)
tic = os.clock()
local rand = Matrix:random(4, 4, 0, 1)
print(Tools.list.tostring(rand:francisOne()))
print("4 x 4 Random:", os.clock() - tic)
tic = os.clock()
rand = Matrix:random(8, 8, 0, 1)
print("Matrix:", rand)
tic = os.clock()
print(Tools.list.tostring(rand:eigenvalues()))
print("8 x 8 Random:", os.clock() - tic)
tic = os.clock()
for i = 1, 100 do
    rand = Matrix:random(8, 8, 0, 1)
    rand:eigenvalues()
end
print("Large Random Test:", (os.clock() - tic) / 100)
tic = os.clock()

print("Strassen Tests:")

local save
tic = os.clock()
save = ONEs * ONEs
print(os.clock() - tic)
tic = os.clock()
save = MAONEs * MAONEs
print(os.clock() - tic)
tic = os.clock()
save = ONES * ONES
print(os.clock() - tic)
tic = os.clock()
save = MAONES * MAONES
print(os.clock() - tic)