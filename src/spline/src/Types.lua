-- Point type
-- Determines ambient space of spline
export type Point = Vector2 | Vector3 | { number }
export type Vector = Vector2 | Vector3 | { number }

-- Basis curve types
-- Independent of ambient space
export type Linear = {
	-- a + bt
	a: { Point },
	b: { Point },
}

export type Quadratic = {
	-- a + bt + ct^2
	a: { Point },
	b: { Point },
	c: { Point },
}

export type Cubic = {
	-- a + bt + ct^2 + dt^3
	a: { Point },
	b: { Point },
	c: { Point },
	d: { Point },
}

export type Polynomial = { { number } }

export type Trigonometric = {}

export type BasisCurveType = Linear | Quadratic | Cubic | Polynomial | Trigonometric

-- Spline type
export type CatmullRom = Cubic & { alpha: number, tension: number }

export type SplineType = Linear | QuadraticBezier | CubicBezier | CatmullRom

-- Base spline
export type Spline = {
	-- Necessary data
	Type: SplineType,

	-- Convenient data for the user
	Points: { Point },

	-- Speedy data
	LUT: {},
	Domains: {},

	-- Global properties
	Closed: boolean,
	Regular: boolean,
	Convex: boolean,
	Length: number?, -- Assigned via CacheLength()
	WindingNumber: number?, -- Assigned via CacheWindingNumber()
	RotationIndex: number?, -- Assigned via CacheRotationIndex()
}

export type RotationalSpline = {
	Curves: { RotationalCurveType },
	Frames: { FrameType },
}

export type CFrameSpline = Spline & RotationalSpline

return nil
