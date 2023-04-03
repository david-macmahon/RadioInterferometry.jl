module TopoUVW

export hadec2uvw, azel2uvw, azza2uvw

import Rotations: RotY, RotZ

"""
    hadec2uvw(hob_rad::Real, dob_rad::Real, lon_rad::Real[, xyz])
    hadec2uvw(hob_rad::Real, dob_rad::Real, obslla[, xyz])

When `xyz` is not passed, return transformation (rotation and permutation)
matrix to convert coordinates from an ITRF aligned (X,Y,Z) frame to a (U,V,W)
frame where U is eastward, V is northward, and W points to the observed hour
angle `hob_rad` (west positive) and observed declination `dob_rad` as seen from
longitude `lon_rad` (all in radians).  An object with a `lon` field in degrees
(e.g. a `Geodesy.LLA` or `NamedTuple`) can be passed instead of `lon_rad`.

This works by rotating the XYZ frame anticlockwise about the Z (i.e. third)
axis by `lon_rad-hob_rad`, producing a (X',U,Z) frame, then rotating that frame
anticlockwise about the U (i.e. second) axis by `-dob_rad`, producing an
(W,U,V) frame which is then permuted to (U,V,W) where U is east, V is north,
and W is in the direction of the given hour angle and declination as seen from
the given longitude.

Left multiplying an ITRF-aligned (X,Y,Z) coordinate vector or matrix by the
returned rotation matrix will result in the corresponding (U,V,W) coordinates
for the given direction and longitude:

    uvw = xyz2uvw(ha, dec, lon) * xyz

When `xyz` is passed, this multiplication is performed and the transformed array
is returned rather than a transformation matrix.
"""
function hadec2uvw(hob_rad, dob_rad, lon_rad::Real)
    # Rotation matrix to transform XYZ frame to UVW frame in direction of
    # (hob, dob) at obslla.lon.
    xyz2uvw = [0 1 0
               0 0 1
               1 0 0] * RotY(dob_rad) * RotZ(hob_rad-lon_rad)

    # xyz to uvw
    xyz2uvw
end

function hadec2uvw(hob_rad, dob_rad, obslla)
    hadec2uvw(hob_rad, dob_rad, deg2rad(obslla.lon))
end

function hadec2uvw(hob_rad, dob_rad, lon_rad::Real, xyz)
    hadec2uvw(hob_rad, dob_rad, lon_rad) * xyz
end

function hadec2uvw(hob_rad, dob_rad, obslla, xyz)
    hadec2uvw(hob_rad, dob_rad, deg2rad(obslla.lon), xyz)
end

function hadec2uvw(hob_rad, dob_rad, lon_rad_or_obslla, antxyz, bls)
    antuvw = hadec2uvw(hob_rad, dob_rad, lon_rad_or_obslla, antxyz)
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function azel2uvw(aob_rad, eob_rad, lat_rad::Real, lon_rad::Real)
    hob, dob = ae2hd(aob_rad, eob_rad, lat_rad)
    hadec2uvw(hob, dob, lon_rad)
end

function azel2uvw(aob_rad, eob_rad, lat_rad::Real, lon_rad::Real, xyz)
    azel2uvw(aob_rad, eob_rad, lat_rad::Real, lon_rad::Real) * xyz
end

function azel2uvw(aob_rad, eob_rad, obslla)
    azel2uvw(aob_rad, eob_rad, deg2rad(obslla.lat), deg2rad(obslla.lon))
end

function azel2uvw(aob_rad, eob_rad, obslla, xyz)
    azel2uvw(aob_rad, eob_rad, deg2rad(obslla.lat), deg2rad(obslla.lon), xyz)
end

function azel2uvw(aob_rad, eob_rad, lat_rad::Real, lon_rad::Real, antxyz, bls)
    antuvw = azel2uvw(aob_rad, eob_rad, lat_rad, lon_rad, antxyz)
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function azel2uvw(aob_rad, eob_rad, obslla, antxyz, bls)
    azel2uvw(aob_rad, eob_rad, deg2rad(obslla.lat), deg2rad(obslla.lon), antxyz, bls)
end

function azza2uvw(aob_rad, zob_rad, lat_rad::Real, lon_rad::Real)
    axel2uvw(aob_rad, π/2-zob_rad, lat_rad, lon_rad)
end

function azza2uvw(aob_rad, zob_rad, lat_rad::Real, lon_rad::Real, xyz)
    axel2uvw(aob_rad, π/2-zob_rad, lat_rad, lon_rad, xyz)
end

function azza2uvw(aob_rad, zob_rad, obslla)
    axel2uvw(aob_rad, π/2-zob_rad, obslla)
end

function azza2uvw(aob_rad, zob_rad, obslla, xyz)
    axel2uvw(aob_rad, π/2-zob_rad, obslla, xyz)
end

function azza2uvw(aob_rad, zob_rad, lat_rad::Real, lon_rad::Real, antxyz, bls)
    axel2uvw(aob_rad, π/2-zob_rad, lat_rad, lon_rad, antxyz, bls)
end

function azza2uvw(aob_rad, zob_rad, obslla, antxyz, bls)
    axel2uvw(aob_rad, π/2-zob_rad, obslla, antxyz, bls)
end

#=
import ..RADec2Obs: radec2obs

function radec2uvw(ra, dec, jd, obslla, antxyz, bls;
                   dut1=nothing, xp=nothing, yp=nothing)
    # Convert nothings to default values
    dut1 = something(dut1, 0.0)
    xp = something(xp, 0.0)
    yp = something(yp, 0.0)

    obs = radec2obs(ra, dec, jd, obslla; dut1, xp, yp)
    hadec2uvw(obs.hob, obs.dob, obslla, antxyz, bls)
end
=#

end # module TopoUVW