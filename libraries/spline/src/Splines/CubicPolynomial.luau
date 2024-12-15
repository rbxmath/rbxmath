local Position = require(script.Parent.Position)
local SplineUtils = require(script.Parent.Parent.SplineUtils)
local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicPolynomial = {}

CubicPolynomial.Interpolant = setmetatable({}, Position.Interpolant)
CubicPolynomial.Interpolant.__index = CubicPolynomial.Interpolant

function CubicPolynomial.Interpolant.new(a: Vector, b: Vector, c: Vector, d: Vector)
	local self = setmetatable(Position.Spline.new(), CubicPolynomial.Interpolant)

	self.a = a
	self.b = b
	self.c = c
	self.d = d

	self.Codimension = SplineUtils.GetCodimensionFromPoint(a)
	self.Length = self:SolveLength()

	return self
end

function CubicPolynomial.Interpolant:SolvePosition(t: number): Point
	t = self:_accountForUnitSpeed(t)

	-- r(t) in Horner's form
	return self.a + t * (self.b + t * (self.c + t * self.d))
end

function CubicPolynomial.Interpolant:SolveVelocity(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r'(t) in Horner's form
	return self.b + t * (2 * self.c + t * 3 * self.d)
end

function CubicPolynomial.Interpolant:SolveAcceleration(t: number): Vector
	t = self:_accountForUnitSpeed(t)

	-- r''(t)
	return 2 * self.c + t * 6 * self.d
end

function CubicPolynomial.Interpolant:SolveJerk(): Vector
	-- r'''(t)
	return 6 * self.d
end

-- TODO: This is bad and does not work when the spline has unit speed.
function CubicPolynomial.Interpolant:SolveBoundingBox(): (Point, Point)
	local pointType = typeof(self.a)
	local boundCandidates = {}

	-- First-derivative coefficients
	local d1 = 3 * self.d
	local c1 = 2 * self.c
	local b1 = self.b

	local candidates = { self:SolvePosition(0) }

	-- The a, b, c here are to coincide with the usual quadratic notation of
	-- at^2 + bt^2 + c
	local function addCandidate(a, b, c)
		local t1, t2 = SplineUtils.SolveQuadratic(a, b, c)

		if t1 ~= nil then
			if t1 >= 0 and t1 <= 1 then
				table.insert(candidates, self:SolvePosition(t1))
			end

			if t2 ~= nil and t2 >= 0 and t2 <= 1 then
				table.insert(candidates, self:SolvePosition(t2))
			end
		end
	end

	if pointType == "Vector2" then
		addCandidate(d1.X, c1.X, b1.X)
		addCandidate(d1.Y, c1.Y, b1.Y)
	elseif pointType == "Vector3" then
		addCandidate(d1.X, c1.X, b1.X)
		addCandidate(d1.Y, c1.Y, b1.Y)
		addCandidate(d1.Z, c1.Z, b1.Z)
	elseif pointType == "table" then
		for i = 1, d1 do
			addCandidate(d1[i], c1[i], b1[i])
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

return CubicPolynomial
