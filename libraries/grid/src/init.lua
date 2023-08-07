local Tools = require(script.Parent.Tools)
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object
local Scalars = require(script.Parent.Scalars)
local Complex = Scalars.Complex
local Vectors = require(script.Parent.Vectors)
local FFT = require(script.Parent.FastFourierTransform)
local Matrices = require(script.Parent.Matrices)
local SparseMatrix = Matrices.SparseMatrix

local Grid = {}

-- This will return an n by m grid where every element is a subgrid. Assumes that n and m divide the length and width!
function Grid.subdivide(grid: Tensor, n: number, m: number): Object
	local length = #grid
	local width = #grid[1]

	-- Initializing grid to be returned
	local returnGrid = {}
	for i = 1, n do
		returnGrid[i] = {}
		for j = 1, m do
			returnGrid[i][j] = {}
		end
	end

	local lengthStride = length / n
	local widthStride = width / m

	for i = 1, lengthStride do
		for j = 1, widthStride do
			for ii = 1, n do
				for jj = 1, m do
					if not returnGrid[ii][jj][i] then
						returnGrid[ii][jj][i] = {}
					end
					returnGrid[ii][jj][i][j] = grid[i + (ii - 1) * lengthStride][j + (jj - 1) * widthStride]
				end
			end
		end
	end

	return returnGrid
end

-- Takes coordinate of n by m subgrid and returns grid coordinate. subgridI and subgrid J determine which chunk you're in, and i and j are the coordinates in that grid.
function Grid.subgridCoordinateToGridCoordinate(
	gridLength: number,
	gridWidth: number,
	n: number,
	m: number,
	subgridI: number,
	subgridJ: number,
	i: number,
	j: number
): (number, number)
	local gridLengthStride = gridLength / n
	local gridWidthStride = gridWidth / m
	local totalLengthStride = gridLengthStride * (subgridI - 1)
	local totalWidthStride = gridWidthStride * (subgridJ - 1)
	return totalLengthStride + i, totalWidthStride + j
end

function Grid.compress(grid: Array<Vector>, tolerance: number): Object
	local compressedGrid = {}
	tolerance = tolerance or 10 ^ -13
	for k, v in ipairs(grid) do
		compressedGrid[k] = FFT:FCT(v)
	end
	local doubleCompressedGrid = {}
	for i = 1, #compressedGrid[1] do
		local data = {}
		for j = 1, #compressedGrid do
			data[j] = compressedGrid[j][i]
		end
		doubleCompressedGrid[i] = FFT:FCT(data)
	end
	return SparseMatrix:new(doubleCompressedGrid, #doubleCompressedGrid, #doubleCompressedGrid[1])
		:transpose()
		:numericalClean(tolerance)
end

function Grid.decompress(sparseGrid: Object): Array<Vector>
	local decompressedGrid = {}
	for i = 1, sparseGrid.length do
		decompressedGrid[i] = FFT:IFCT(sparseGrid:getColumn(i))
	end
	local transposedGrid = {}
	for i = 1, #decompressedGrid do
		for j = 1, #decompressedGrid[1] do
			if not transposedGrid[j] then
				transposedGrid[j] = {}
			end
			transposedGrid[j][i] = decompressedGrid[i][j]
		end
	end
	local doubleDecompressedGrid = {}
	for k, v in ipairs(transposedGrid) do
		doubleDecompressedGrid[k] = FFT:IFCT(v)
	end
	return doubleDecompressedGrid
end

return Grid
