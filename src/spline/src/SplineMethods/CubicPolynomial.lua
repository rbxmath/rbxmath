local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local CubicPolynomial = {}

function CubicPolynomial.SolvePosition(self, t: number): Point
	-- r(t) in Horner's form
	return self.a + t * (self.b + t * (self.c + t * self.d))
end

function CubicPolynomial.SolveVelocity(self, t: number): Vector
	-- r'(t) in Horner's form
	return self.b + t * (2 * self.c + t * 3 * self.d)
end

function CubicPolynomial.SolveAcceleration(self, t: number): Vector
	-- r''(t)
	return 2 * self.c + t * 6 * self.d
end

function CubicPolynomial.SolveJerk(self): Vector
	-- r'''(t)
	return 6 * self.d
end

return CubicPolynomial
