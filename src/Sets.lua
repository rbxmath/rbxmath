--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Tools = require("src/Tools")
type Vector = Tools.Vector
type Array<T> = Tools.Array<T>
type ScalarFunction = Tools.ScalarFunction
type Tensor = Tools.Tensor
type Object = Tools.Object

local Sets = {}

Sets.Set = {
	cardinality = 0,
	data = {},
}
local Set = Sets.Set

local function removeDuplicatesFromSortedList(array: Object): Object
	local result = { array[1] }

	for i = 2, #array do
		if array[i - 1] ~= array[i] then
			result[#result + 1] = array[i]
		end
	end

	return result
end

function Set:new(array: Object): Object
	local list = Tools.list.copy(array)
	table.sort(list)
	list = removeDuplicatesFromSortedList(list)
	local o = {}
	setmetatable(o, self)
	self.__index = self
	o.data = list
	o.cardinality = #list
	return o
end

function Set:copy(): Object
	local list = Tools.list.copy(self.data)
	local o = {}
	setmetatable(o, self)
	self.__index = self
	o.data = list
	o.cardinality = #list
	return o
end

function Set.__add(left, right)
	local result = Set:new({})

	local i = 0
	local j = 0
	local index = 1

	while i < left.cardinality or j < right.cardinality do
		if i == left.cardinality then
			result[index] = right[j + 1]
			j = j + 1
			index = index + 1
		elseif j == right.cardinality then
			result[index] = left[i + 1]
			i = i + 1
			index = index + 1
		elseif left[i + 1] < right[j + 1] then
			result[index] = left[i + 1]
			i = i + 1
			index = index + 1
		elseif left[i + 1] == right[j + 1] then
			result[index] = left[i + 1]
			i = i + 1
			j = j + 1
			index = index + 1
		else
			result[index] = right[j + 1]
			j = j + 1
			index = index + 1
		end
	end

	return result
end

function Set.__mul(left, right)
	local result = Set:new({})

	if left == nil or right == nil or left.cardinality == 0 or right.cardinality == 0 then
		return result
	end

	local i = 1
	local j = 1

	while i <= left.cardinality and j <= right.cardinality do
		while i <= left.cardinality and j <= right.cardinality and left.data[i] < right.data[j] do
			i = i + 1
		end
		while i <= left.cardinality and j <= right.cardinality and left.data[i] > right.data[j] do
			j = j + 1
		end
		if left.data[i] == right.data[j] then
			result[#result + 1] = left.data[i]
			i = i + 1
			j = j + 1
		end
	end

	return result
end

function Set.__eq(left, right)
	if left.cardinality ~= right.cardinality then
		return false
	end

	for i = 1, left.cardinality do
		if left.data[i] ~= right.data[i] then
			return false
		end
	end

	return true
end

function Set.__lt(left, right)
	if left.cardinality >= right.cardinality then
		return false
	end

	return left * right == left
end

function Set.__le(left, right)
	return left * right == left
end

function Set.__tostring(left)
	return Tools.list.tostring(left.data)
end

function Set:getIndex(item)
	return Tools.list.binarySearch(item, self.data)
end

function Set:addTo(item)
	local index = Tools.list.binarySearchBetween(item, self.data)
	if type(index) ~= "table" then
		local temp = self.data[index]
		self.data[index] = item
		for i = index + 1, self.cardinality + 1, 1 do
			self.data[i], temp = temp, self.data[i]
		end
		self.cardinality += 1
	end
	return self
end

return Sets
