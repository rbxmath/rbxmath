--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require(script.Parent.Tools)
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object
local Matrices = require(script.Parent.Matrices)
local SparseMatrix = Matrices.SparseMatrix
local Sets = require(script.Parent.Sets)
local Set = Sets.Set

local Graphs = {}

Graphs.Graph = {
	adjacencyMatrix = nil,
	numberOfVertices = 0,
	numberOfEdges = 0,
	labelTable = {},
	vertexSet = nil,
}

local Graph = Graphs.Graph

function Graph:new(edgeList: Object, vertexNames: Object): Object
	local o = {}
	setmetatable(o, self)
	self.__index = self
	vertexNames = vertexNames or {}
	local vertexSet = Set:new(vertexNames)
	for _, v in pairs(edgeList) do
		vertexSet:addTo(v[1]):addTo(v[2])
	end
	local labelTable = {}
	for k, v in ipairs(vertexSet) do
		labelTable[v] = k
	end
	local adjacencyMatrix = SparseMatrix:new({}, vertexSet.cardinality, vertexSet.cardinality)
	local numberOfVertices = #vertexSet
	local numberOfEdges = 0
	for _, v in pairs(edgeList) do
		local left = labelTable[v[1]]
		local right = labelTable[v[2]]
		if not adjacencyMatrix:get(left, right) then
			numberOfEdges += 1
			adjacencyMatrix:set(left, right, 1)
			adjacencyMatrix:set(right, left, 1)
		end
	end
	o.adjacencyMatrix = adjacencyMatrix
	o.labelTable = labelTable
	o.numberOfVertices = numberOfVertices
	o.numberOfEdges = numberOfEdges
	o.vertexSet = vertexSet
end

function Graph:edgeList()
	local adjacencyMatrix = self.adjacencyMatrix
	local length = adjacencyMatrix.length
	local width = adjacencyMatrix.width
	local edgeList = {}
	for k, v in pairs(adjacencyMatrix) do
		local c = k % width
		if c == 0 then
			c = width
		end
		local r = math.floor((k - c) / length)
		if r >= c then
			if v ~= 0 then
				edgeList[#edgeList + 1] = table.sort({ c, r })
			end
		end
	end
	return edgeList
end

return Graphs
