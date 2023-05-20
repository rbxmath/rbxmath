local Types = require(script.Parent.Parent.Types)

type Point = Types.Point
type Vector = Types.Vector

local QuadraticPolynomial = {}

function QuadraticPolynomial.SolvePosition(self, t: number): Point
	-- r(t) in Horner's form
	return self.a + t * (self.b + t * self.c)
end

function QuadraticPolynomial.SolveVelocity(self, t: number): Vector
	-- r'(t) in Horner's form
	return self.b + t * 2 * self.c
end

function QuadraticPolynomial.SolveAcceleration(self): Vector
	-- r''(t)
	return 2 * self.c
end

function QuadraticPolynomial.SolveJerk(self): Vector
	-- r'''(t)
	return 0 * self.a
end

return QuadraticPolynomial
