--[[
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
]]

local Noise = {}

-- stylua: ignore
local perm = {
	[0] = 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
	140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247,
	120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177,
	33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71,
	134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133,
	230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161,
	1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130,
	116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250,
	124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227,
	47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44,
	154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98,
	108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34,
	242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14,
	239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121,
	50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243,
	141, 128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131,
	13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37,
	240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252,
	219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125,
	136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158,
	231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245,
	40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187,
	208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198,
	173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126,
	255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223,
	183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167,
	43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185,
	112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179,
	162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199,
	106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236,
	205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156,
	180,
}

local function grad(hash, x, y, z, t)
	local h = bit32.band(hash, 31)
	local u = h < 24 and x or y
	local v = h < 16 and y or z
	local w = h < 8 and z or t
	return ((bit32.band(h, 1) > 0) and -u or u)
		+ ((bit32.band(h, 2) > 0) and -v or v)
		+ ((bit32.band(h, 4) > 0) and -w or w)
end

local function fade(t)
	return t * t * t * (t * (t * 6 - 15) + 10)
end

local function lerp(t, a, b)
	return a + t * (b - a)
end

--[=[
	Generates 4D Perlin noise. Ported from Noise1234 by Stefan Gustavson
	(https://github.com/stegu/perlin-noise). Empirically, noise values are in
	the range [0, 1] with a mean of 0.5 and standard deviation of 0.129.

	@param x -- The x-coordinate to sample from
	@param y -- The y-coordinate to sample from
	@param z -- The z-coordinate to sample from
	@param w -- The w-coordinate to sample from
	@return -- A noise value
]=]
function Noise.Perlin4D(x: number, y: number, z: number, w: number): number
	local ix0 = math.floor(x) -- Integer part of x
	local iy0 = math.floor(y) -- Integer part of y
	local iz0 = math.floor(z) -- Integer part of z
	local iw0 = math.floor(w) -- Integer part of w
	local fx0 = x - ix0 -- Fractional part of x
	local fy0 = y - iy0 -- Fractional part of y
	local fz0 = z - iz0 -- Fractional part of z
	local fw0 = w - iw0 -- Fractional part of w
	local fx1 = fx0 - 1
	local fy1 = fy0 - 1
	local fz1 = fz0 - 1
	local fw1 = fw0 - 1
	local ix1 = bit32.band(ix0 + 1, 255)
	local iy1 = bit32.band(iy0 + 1, 255)
	local iz1 = bit32.band(iz0 + 1, 255)
	local iw1 = bit32.band(iw0 + 1, 255)
	ix0 = bit32.band(ix0, 255)
	iy0 = bit32.band(iy0, 255)
	iz0 = bit32.band(iz0, 255)
	iw0 = bit32.band(iw0, 255)

	local q = fade(fw0)
	local r = fade(fz0)
	local t = fade(fy0)
	local s = fade(fx0)

	local nxyz0 = grad(perm[ix0 + perm[iy0 + perm[iz0 + perm[iw0]]]], fx0, fy0, fz0, fw0)
	local nxyz1 = grad(perm[ix0 + perm[iy0 + perm[iz0 + perm[iw1]]]], fx0, fy0, fz0, fw1)
	local nxy0 = lerp(q, nxyz0, nxyz1)

	nxyz0 = grad(perm[ix0 + perm[iy0 + perm[iz1 + perm[iw0]]]], fx0, fy0, fz1, fw0)
	nxyz1 = grad(perm[ix0 + perm[iy0 + perm[iz1 + perm[iw1]]]], fx0, fy0, fz1, fw1)
	local nxy1 = lerp(q, nxyz0, nxyz1)

	local nx0 = lerp(r, nxy0, nxy1)

	nxyz0 = grad(perm[ix0 + perm[iy1 + perm[iz0 + perm[iw0]]]], fx0, fy1, fz0, fw0)
	nxyz1 = grad(perm[ix0 + perm[iy1 + perm[iz0 + perm[iw1]]]], fx0, fy1, fz0, fw1)
	nxy0 = lerp(q, nxyz0, nxyz1)

	nxyz0 = grad(perm[ix0 + perm[iy1 + perm[iz1 + perm[iw0]]]], fx0, fy1, fz1, fw0)
	nxyz1 = grad(perm[ix0 + perm[iy1 + perm[iz1 + perm[iw1]]]], fx0, fy1, fz1, fw1)
	nxy1 = lerp(q, nxyz0, nxyz1)

	local nx1 = lerp(r, nxy0, nxy1)

	local n0 = lerp(t, nx0, nx1)

	nxyz0 = grad(perm[ix1 + perm[iy0 + perm[iz0 + perm[iw0]]]], fx1, fy0, fz0, fw0)
	nxyz1 = grad(perm[ix1 + perm[iy0 + perm[iz0 + perm[iw1]]]], fx1, fy0, fz0, fw1)
	nxy0 = lerp(q, nxyz0, nxyz1)

	nxyz0 = grad(perm[ix1 + perm[iy0 + perm[iz1 + perm[iw0]]]], fx1, fy0, fz1, fw0)
	nxyz1 = grad(perm[ix1 + perm[iy0 + perm[iz1 + perm[iw1]]]], fx1, fy0, fz1, fw1)
	nxy1 = lerp(q, nxyz0, nxyz1)

	nx0 = lerp(r, nxy0, nxy1)

	nxyz0 = grad(perm[ix1 + perm[iy1 + perm[iz0 + perm[iw0]]]], fx1, fy1, fz0, fw0)
	nxyz1 = grad(perm[ix1 + perm[iy1 + perm[iz0 + perm[iw1]]]], fx1, fy1, fz0, fw1)
	nxy0 = lerp(q, nxyz0, nxyz1)

	nxyz0 = grad(perm[ix1 + perm[iy1 + perm[iz1 + perm[iw0]]]], fx1, fy1, fz1, fw0)
	nxyz1 = grad(perm[ix1 + perm[iy1 + perm[iz1 + perm[iw1]]]], fx1, fy1, fz1, fw1)
	nxy1 = lerp(q, nxyz0, nxyz1)

	nx1 = lerp(r, nxy0, nxy1)

	local n1 = lerp(t, nx0, nx1)

	return (0.87 * lerp(s, n0, n1)) / 2 + 0.5
end

--[=[
	Generates a tileable 2D array of Perlin noise by sampling on a Clifford
	torus in R^4.

	@param width -- The width of the tile
	@param height -- The height of the tile
	@param origin -- The 4D offset for sampling noise. Defaults to {0, 0, 0, 0}.
	@param scale -- The size of the noise. Higher scale gives smaller noise. Defaults to 1.
	@return -- A tileable 2D array of Perlin noise. The noise value at (x, y) is indexed by tile[x][y].
]=]
function Noise.GetNoiseTile(width: number, height: number, origin: { number }?, scale: number?): { { number } }
	origin = origin or { 0, 0, 0, 0 }
	scale = scale or 1

	local tile = table.create(width)

	for x = 1, width do
		local s = (x - 1) / width

		local column = table.create(height)
		tile[x] = column

		for y = 1, height do
			local t = (y - 1) / height

			-- Noise value on Clifford torus in 4D
			column[y] = Noise.Perlin4D(
				origin[1] + scale * math.cos(s * 2 * math.pi),
				origin[2] + scale * math.cos(t * 2 * math.pi),
				origin[3] + scale * math.sin(s * 2 * math.pi),
				origin[4] + scale * math.sin(t * 2 * math.pi)
			)
		end
	end

	return tile
end

return Noise
