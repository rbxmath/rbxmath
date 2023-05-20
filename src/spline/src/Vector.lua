local Vector = {}

Vector.__index = function(self, key)
	if key == "Magnitude" then
		local magnitude = 0
		for _, coordinate in self do
			magnitude += coordinate ^ 2
		end
		return math.sqrt(magnitude)
	elseif key == "Unit" then
		return self / self.Magnitude
	else
		return Vector[key]
	end
end

function Vector.new(coordinates: { number })
	return setmetatable(coordinates, Vector)
end

function Vector.isVector(value)
	return getmetatable(value) == Vector
end

function Vector:Dot(other): number
	assert(Vector.isVector(other), "Bad other")
	assert(#self == #other, "Cannot Dot vectors in different dimensions")

	local dot = 0

	for i, coord in self do
		dot += coord * other[i]
	end

	return dot
end

function Vector.__unm(self)
	local newCoordinates = table.create(#self)

	for i, coord in self do
		newCoordinates[i] = -coord
	end

	return Vector.new(newCoordinates)
end

function Vector.__add(a, b)
	assert(Vector.isVector(a) and Vector.isVector(b), "Attempted Vector subtraction with non-Vector")
	assert(#a == #b, "Cannot add Vectors in different dimensions")

	local newCoordinates = table.create(#a)
	for i, coord in a do
		newCoordinates[i] = coord + b[i]
	end

	return Vector.new(newCoordinates)
end

function Vector.__sub(a, b)
	assert(Vector.isVector(a) and Vector.isVector(b), "Attempted Vector addition with non-Vector")
	assert(#a == #b, "Cannot subtract Vectors in different dimensions")

	local newCoordinates = table.create(#a)
	for i, coord in a do
		newCoordinates[i] = coord - b[i]
	end

	return Vector.new(newCoordinates)
end

function Vector.__mul(a, b)
	if Vector.isVector(a) and type(b) == "number" then
		local newCoordinates = table.create(#a)

		for i, coord in a do
			newCoordinates[i] = coord * b
		end

		return Vector.new(newCoordinates)
	elseif type(a) == "number" and Vector.isVector(b) then
		local newCoordinates = table.create(#b)

		for i, coord in b do
			newCoordinates[i] = coord * a
		end

		return Vector.new(newCoordinates)
	else
		error("Attempted Vector multiplication with non-number")
	end
end

function Vector.__div(a, b)
	if Vector.isVector(a) and type(b) == "number" then
		local newCoordinates = table.create(#a)

		for i, coord in a do
			newCoordinates[i] = coord / b
		end

		return Vector.new(newCoordinates)
	else
		error("Bad Vector division")
	end
end

return Vector
