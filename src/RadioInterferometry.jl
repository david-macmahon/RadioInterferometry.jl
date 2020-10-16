"""
A module for radio interferometry that uses ERFA for much of its underlying
functionality.
"""
module RadioInterferometry

export I3
export xyz2uvw

using ERFA
import ERFA.DD2R
import ERFA.DR2D
import ERFA.CMPS

using Geodesy

"""
`I3` is a 3x3 identity matrix.
"""
I3 = Float64[1 0 0
             0 1 0
             0 0 1]

# Add default I3 argument to ERFA.rx, ERFA.ry, ERFA.rz
ERFA.rx(phi::Real) = ERFA.rx(phi, I3)
ERFA.ry(phi::Real) = ERFA.ry(phi, I3)
ERFA.rz(phi::Real) = ERFA.rz(phi, I3)

# Some additional ways of expressing speed of light
"""
Speed of light in a vacuum (C) in meters/nanosecond
"""
CMPNS = CMPS / 1e9
"""
Speed of light in a vacuum (C) in nanoseconds/meter
"""
CNSPM = 1e9 / CMPNS
"""
Speed of light in a vacuum (C) in seconds/meter
"""
CSPM  = 1 / CMPS


"""
    xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert coordinates
in ITRF aligned (X,Y,Z) frame to a (U,V,W) projection with W pointing to the
hour angle `ha_rad` (west positive) and declination `dec_rad` as seen from
longitude `lon_rad` (all in radians).

This works by rotating the XYZ frame anticlockwise about the Z (i.e. third)
axis by `lon_rad-ha_rad`, producing a (X',U,Z) frame, then rotating that frame
anticlockwise about the U (i.e. second) axis by `-dec_rad`, producing an
(W,U,V) frame which is then permuted to (U,V,W) where U is east, V is north,
and W is in the direction of projection.  Left multiplying an ITRF-aligned
(X,Y,Z) coordinate vector by the returned rotation matrix will result in the
(U,V,W) projection of the (X,Y,Z) coordinates in the specified direction at the
specified longitude:

    uvw = xyz2uvw(lon, ha, dec) * xyz

!!! warning
    A final rotation around W for the paralactic angle is likely to be needed
    to align the U and V axes with RA/Dec directions.
"""
function xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  [0 1 0
   0 0 1
   1 0 0] * ERFA.ry(-dec_rad,ERFA.rz(lon_rad-h_rad))
end

"""
    xyz2uvw(xyz::Union{Array{<:Real,1},Array{<:Real,2},ECEF},
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}

    xyz2uvw(x::Real, y::Real, z::Real,
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}

Return the coordinataes of point(s) `xyz` (or [x,y,z]) projected in the
direction of hour angle `ha_rad` (west positive) and `dec_rad` from longitude
`lon_rad`, all in radians, and permuted to a (U,V,W) frame where U is east, V
is north, and W is in the direction of projection.
"""
function xyz2uvw(xyz::Union{Array{<:Real,1},Array{<:Real,2},ECEF},
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * xyz
end

function xyz2uvw(x::Real, y::Real, z::Real,
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * [x, y, z]
end

end # module
