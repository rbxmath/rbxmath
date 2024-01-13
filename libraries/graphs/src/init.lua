--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Packages = script.Parent
local Tools = require(Packages.Tools)	
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
   type ScalarMap = Tools.ScalarMap
type Tensor = Tools.Tensor
type Object = Tools.Object
local Matrices = require(Packages.Matrices)
local SparseMatrix = Matrices.SparseMatrix
local Sets = require(Packages.Sets)
local Set = Sets.Set

local Graphs = {}

Graphs.Graph = {
   adjacencyMatrix = nil,
   labelTable = {},
   vertexSet = nil,
}

local Graph = Graphs.Graph

--[=[
   Creates a Graph object
   
   @param edgeList --- A list of edges between vertices in the graph. Edges should be of
                       the form {v1, v2, d, w} where d is an optional string (either 'd'
                       or 'u') to determine if the edge is directed or undirected and w
                       is an optional number for the edge weight.
--]=]
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
   local adjacencyMatrix = SparseMatrix:new({}, vertexSet.cardinality,
      vertexSet.cardinality)
   for _, v in pairs(edgeList) do
      local left = labelTable[v[1]]
      local right = labelTable[v[2]]
      local w = v[4] or 1
      if v[3] == 'd' then
	 if not adjacencyMatrix:get(left, right) then
	    adjacencyMatrix:set(left, right, w)
	 end
      else
	 if not adjacencyMatrix:get(left, right) then
	    adjacencyMatrix:set(left, right, w)
	    adjacencyMatrix:set(right, left, w)
	 end
      end
   end
   o.adjacencyMatrix = adjacencyMatrix
   o.labelTable = labelTable
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
      local val1, val2 = adjacencyMatrix:get(r, c), adjacencyMatrix:get(c, r)
      if val1 == val2 then
	 if r >= c and val1 ~= 0 then
	    edgeList[#edgeList + 1] = {r, c, 'u', val1}
	 end
      else
	 edgeList[#edgeList + 1] = {r, c, 'd', val1}
      end
   end
   return edgeList
end

function Graph:numberOfEdges()
   local adjacencyMatrix = self.adjacencyMatrix
   local length = adjacencyMatrix.length
   local width = adjacencyMatrix.width
   local numEdge, numDEdge = 0, 0
   for k, v in pairs(adjacencyMatrix) do
      local c = k % width
      if c == 0 then
	 c = width
      end
      local r = math.floor((k - c) / length)
      local val1, val2 = adjacencyMatrix:get(r, c), adjacencyMatrix:get(c, r)
      if val1 == val2 then
	 if r >= c and val1 ~= 0 then
	    numEdge += 1
	 end
      else
	 numDEdge += 1
      end
   end
   return numEdge, numDEdge
end

function Graph:vertexList()
   return self.vertexSet.data
end

function Graph:getIndex(vertex)
   return self.labelTable[vertex]
end

function Graph:isEdgeBetween(v1, v2)
   return self.adjacencyMatrix:get(self:getIndex(v1), self:getIndex(v2)) ~= 0
      and self.adjacencyMatrix:get(self:getIndex(v1), self:getIndex(v2)) ~= 0
end

function Graph:isEdgeFrom(v1, v2)
   return self.adjacencyMatrix:get(self:getIndex(v1), self:getIndex(v2)) ~= 0
end

function Graph:getEdgeBetween(v1, v2)
   return self.adjacencyMatrix:get(self:getIndex(v1), self:getIndex(v2)),
      self.adjacencyMatrix:get(self:getIndex(v2), self:getIndex(v1))

function Graph:getEdgeFrom(v1, v2)
   return self.adjacencyMatrix:get(self:getIndex(v1), self:getIndex(v2))
end

function Graph:setEdgeBetween(v1, v2, w)
   self.adjacencyMatrix:set(self:getIndex(v1), self:getIndex(v2), w)
end

function Graph:setEdgeFrom(v1, v2, w)
   self.adjacencyMatrix:set(self:getIndex(v1), self:getIndex(v2), w)
end

function Graph:getNeighbors(vertex, stateMatrix)
   local index = self:getIndex(vertex)
   local adjacencyMatrix = stateMatrix or self.adjacencyMatrix
   local neighborList = {}
   for i = 1, self.vertexSet.cardinality do
      local get = adjacencyMatrix:get(index, i)
      if get ~= 0 then
	 neighborList[#neighborList + 1] = get
      end
   end
   return neighborList, neighborIndexList
end

function Graph:isolateVertex(vertex, stateMatrix)
   local index = self:getIndex(vertex)
   -- Note the lack of a copy here!
   local adjacencyMatrix = stateMatrix or self.adjacencyMatrix
   for i = 1, self.vertexSet.cardinality do
      adjacencyMatrix:set(index, i, 0):set(i, index, 0)
   end
end

function Graph:DFS(stopFunction, startVertex, state)
   startVertex = startVertex or self.vertexSet[1]
   local stop = stopFunction(startVertex)
   if not stop then
      state = state or self.adjacencyMatrix:copy()
      local neighbors = self:getNeighbors(startVertex, state)
      self:isolateVertex(startVertex, state)
      local tail
      for k, v in ipairs(neighbors) do
	 tail = self:DFS(stopFunction, v, state)
	 if type(tail) == "table" then
	    break
	 end
      end
      if type(tail) == "table" then
	 return Tools.list.append({startVertex}, tail)
      else
	 return false
      end
   end
   return {startVertex}
end

function Graph:BFS(stopFunction, startVertex, state)
   startVertex = startVertex or self.vertexSet[1]
   local stop = stopFunction(startVertex)
   if not stop then
      state = state or self.adjacencyMatrix:copy()
      local neighbors = self:getNeighbors(startVertex, state)
      self:isolateVertex(startVertex, state)
      for k, v in ipairs(neighbors) do
	 stop = stopFunction(v)
	 if stop then
	    return {startVertex, v}
	 end
      end
   end
end

Graphs.Tree = {
   parent = nil,
   children = {},
   vertex = nil
}

local Tree = Graphs.Tree

--[=[
   This constructs a tree from a nested table. The expected format is {root, {children}}. 
   So, for instance, {1, {2, {3, 4}, 5, {6}}} would be the tree

     1
     |\
     2 5
     |\ \
     3 4 6

   @param nestedTable --- A table in the recursive form {root, {children}}.
--]=]
function Tree:new(nestedTable, parent)
   local o = {}
   setmetatable(o, self)
   self.__index = self
   o.parent = parent
   o.vertex = nestedTable[1]
   if #nestedTable == 1 then
      return o
   else
      local i = 2
      if type(nestedTable[2]) == "table" then
	 o.children = self:new(nestedTable[2], o)
	 if #nestedTable == 2 then
	    return o
	 end
	 i = 3
      end
      local parentsChildren = {o}
      while i <= #nestedTable do
	 if type(nestedTable[i + 1]) == "table" then
	    parentsChildren[#parentsChildren + 1] = self:new({nestedTable[i], nestedTable[i + 1]}, parent)
	    i += 1
	 else
	    parentsChildren[#parentsChildren + 1] = self:new({nestedTable[i]}, parent)
	 end
	 i += 1
      end
      return parentsChildren
   end
end

function Tree:addVertex(vertex, parent, children)
   
