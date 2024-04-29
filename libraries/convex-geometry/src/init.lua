--[[
   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Convex = {}

local SimplexMesh = {
   simplexDimension = 0,
   adjacencyTable = {},
   vertexTable = {},
   boundary = nil
}

function SimplexMesh:new(simD, adjT, verT, bnd)
   local o = {}
   setmetatable(o, self)
   self._index = self
   o.simplexDimension = simD
   o.adjacencyTable   = adjT
   o.vertexTable      = verT
   o.boundary         = bnd
   return o
end

function SimplexMesh:cycle(vertexList)
   local adjacencyTable = {}
   if #vertexList == 0 then
      return SimplexMesh:new(1, adjacencyTable, vertexTable, nil)
   elseif #vertexList == 1 then
      adjacencyTable[1] = {1, 1}
   else
      adjacencyTable[1] = {2, #vertexList}
      for i = 2, #vertexList - 1 do
	 adjacencyTable[vertexList[i]] = {vertexList[i + 1], vertexList[i - 1]}
      end
      adjacencyTable[vertexList[#vertexList]] = {vertexList[1], vertexList[#vertexList - 1]}
   end
   return SimplexMesh:new(1, adjacencyTable, vertexTable, nil)
end

function SimplexMesh:triangle(pointList)
   return SimplexMesh:new(2, {{}}, SimplexMesh:cycle({1, 2, 3}, pointList))
end

function SimplexMesh:DelaunayMerge(tri1, tri2)
   

--[[
function SimplexMesh:simplex(vertexList, vertexTable)
   -- We need only construct the boundary of the simplex as we have already input all of the information needed to make a simplex!
   for i = #vertexList, 1, -1 do
      local
      for j = 1, #vertexList do
	 if i == j then
	    continue
	 else
	    
--]]
