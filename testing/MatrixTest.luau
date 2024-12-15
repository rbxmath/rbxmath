local Tools = require("src.Tools")
local Matrices = require("src.Matrices")
local Matrix = Matrices.Matrix
local SparseMatrix = Matrices.SparseMatrix
local ComplexMatrix = Matrices.ComplexMatrix
local MA = require("src.MatrixAlgebra")
local Scalars = require("src.Scalars")
local Complex = Scalars.Complex

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
--[[tic = os.clock()
for i = 1, 100 do
    rand = Matrix:random(8, 8, 0, 1)
    rand:eigenvalues()
end
print("Large Random Test:", (os.clock() - tic) / 100)
tic = os.clock()]]

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

print("\nSparse Matrix Tests:")

local iden2
local idenP

tic = os.clock()
iden = SparseMatrix:identity(64)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)
--[[print(iden)
print(iden2)
print(idenP)]]

tic = os.clock()
iden = Matrix:identity(64)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)
--[[print(iden)
print(iden2)
print(idenP)]]

print("\nGrowth Test:")

tic = os.clock()
iden = SparseMatrix:identity(4)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(8)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(16)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(32)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(64)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(128)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(256)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = SparseMatrix:identity(512)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

print("\nSubmatrix Tests:")
local two = SparseMatrix:new({{2, 0}, {0, 2}})
iden:copy()
print(iden:submatrix(1, 4, 1, 4))
print(two)
print(iden:copy():setSubmatrix(1, 1, two):submatrix(1, 2, 1, 2))

print(Tools.admin.makeBanners("Testing", "The FitnessGram Pacer Test is a multistage aerobic capacity test that progressively gets more difficult as it continues. The 20 meter pacer test will begin in 30 seconds. Line up at the start. The running speed starts slowly, but gets faster each minute after you hear this signal. A single lap should be completed each time you hear this sound.  Remember to run in a straight line, and run as long as possible. The second time you fail to complete a lap before the sound, your test is over. The test will begin on the word start. On your mark, get ready, start."))
print(Tools.admin.makeBanners("Matrix Utilities", "This section contains many useful functions including those to copy and manipulate matrices. Functions of the form \"to_____\" or \"set_____\" will change the underlying matrix, while others will return a shallow copy."))
print(Tools.admin.makeBanners("Matrix Metamethods", "This section contains all of the metamethods for matrices. Addition and subtraction are relatively standard, but multiplication is an implementation of Strassen's method. The size of matrix at which Strassen multiplication will be used is set in Matrices.constant.STRASSENLIMIT."))
print(Tools.admin.makeBanners("Common Matrices", "This section contains methods for constructing many common matrices such as the identity and zero matrices."))
print(Tools.admin.makeBanners("Matrix Maps", "This section contains the common maps that send matrices to matrices. This includes functions like the transpose and inverse."))
print(Tools.admin.makeBanners("Linear Systems", "This section contains methods pertaining to solving systems of linear equations. This includes linear solve and the LU factorization."))
print(Tools.admin.makeBanners("Scalar Maps", "This section contains many common scalar maps including the trace and determinant."))
print(Tools.admin.makeBanners("Lua Math: Matrices", "This portion of the library focuses on the creation of many crucial objects in linear algebra. For most users, the Matrix object will be the most useful. Matrices can be multiplied, added, and subtracted with the standard symbols, and their contents can be accessed as you would access information from a table of tables. For the adventorous, there are the ComplexMatrix and SparseMatrix objects. Each of these has their own quirks, but they are all compatible with the standard operations (+, -, *)."))
print(Tools.admin.makeBanners("Eigenvalue Computations", "This section contains the math needed to compute the eigenvalues of a matrix. The primary tool for this is the Francis algorithm. This has been implemented for the Rayleigh shifts of degree 1, 2, 3, and 6. Much of the algorithms in this section come from the book Fundamentals of Matrix Computation by Watkins."))
print(Tools.admin.makeBanners("Matrix Constants", "This section contains some tools for computing constants related to matrix computation."))
print(Tools.admin.makeBanners("WARNING!!!", "While Matrix can be treated roughly like a table of tables in terms of data access, SparseMatrix uses a different data structure to store the values. As such, users should use SparseMatrix:set and SparseMatrix:get to get and set values."))
print(Tools.admin.makeBanners("Matrices", "Matrices are the bread-and-butter object of this library. The underlying data structure is a table of tables, so elements of a matrix can be accessed and set by two table accesses: M[i][j]. Matrices are mutable and many functions will overwrite the matrix itself for the sake of speed."))
print(Tools.admin.makeBanners("Sparse Matrices", "Sparse matrices are the go-to option if your problem involves a matrix with a lot of zero entries. Sparse matrices are significantly faster if the number of non-zero entries is low, but are, in general, slower to use than regular matrices. In terms of memory, sparse matrices are always cheaper to store than their traditional counterparts, so it may be worth the hit in CPU performance."))

print("\nComplex Matrices:\n")

tic = os.clock()
iden = SparseMatrix:identity(16)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = Matrix:identity(16)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = ComplexMatrix:identity(16)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = ComplexMatrix:identity(32)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = ComplexMatrix:identity(64)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

tic = os.clock()
iden = ComplexMatrix:identity(128)
iden2 = iden * iden
idenP = iden + iden
print(os.clock() - tic)

local ones = ComplexMatrix:new({{1,1},{1,1}},true)
print(ones * ones)
rand = ComplexMatrix:random(4)
print(ComplexMatrix:new({{Complex:new(0,1), Complex:new(0,1), Complex:new(0,1)}}))
print(ComplexMatrix:new({{Complex:new(0,1), Complex:new(0,1), Complex:new(0,1)}}):conjugateTranspose())
print(rand)
print(rand:hessenbergForm())