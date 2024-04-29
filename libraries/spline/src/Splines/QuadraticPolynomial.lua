local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticPolynomial = {}

QuadraticPolynomial.Interpolant = setmetatable({}, Position.Interpolant)
QuadraticPolynomial.Interpolant.__index = QuadraticPolynomial.Interpolant

function QuadraticPolynomial.Interpolant.new(a: Vector, b: Vector, c: Vector)
	local self = setmetatable(Position.Interpolant.new(), QuadraticPolynomial.Interpolant)

	self.a = a
	self.b = b
	self.c = c

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function QuadraticPolynomial.Interpolant:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * self.c)
end

function QuadraticPolynomial.Interpolant:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * 2 * self.c
end

function QuadraticPolynomial.Interpolant:SolveAcceleration(): Vector
	-- r''(t)
	return 2 * self.c
end

function QuadraticPolynomial.Interpolant:SolveJerk(): Vector
	-- r'''(t)
	return 0 * self.a
end

-- TODO: This is bad and does not work when the spline has unit speed.
function QuadraticPolynomial.Interpolant:SolveBoundingBox(): (Vector, Vector)
	local pointType = typeof(self.a)
	local boundCandidates = {}

	-- First-derivative coefficients
	local c1 = 2 * self.c
	local b1 = self.b

	local candidates = { self:SolvePosition(0) }

	-- The a and b here are to coincide with the usual linear notation of at + b
	local function addCandidate(a, b)
		local t = SplineUtils.SolveLine(a, b)

		if t ~= nil and t >= 0 and t <= 1 then
			table.insert(candidates, self:SolvePosition(t))
		end
	end

	if pointType == "Vector2" then
		addCandidate(c1.X, b1.X)
		addCandidate(c1.Y, b1.Y)
	elseif pointType == "Vector3" then
		addCandidate(c1.X, b1.X)
		addCandidate(c1.Y, b1.Y)
		addCandidate(c1.Z, b1.Z)
	elseif pointType == "table" then
		for i = 1, c1 do
			addCandidate(c1[i], b1[i])
		end
	else
		error("Bad point type")
	end

	local min = self:SolvePosition(1)
	local max = min

	for _, candidate in candidates do
		min = min:Min(candidate)
		max = max:Max(candidate)
	end

	return min, max
end

return QuadraticPolynomial
