"""
A module for radio interferometry that uses ERFA for much of its underlying
functionality.
"""
module RadioInterferometry

export I3
export dms2d
export hms2h
export xyz2uvw
export enu2uvw

using ERFA
import ERFA.DD2R
import ERFA.DR2D
import ERFA.CMPS

import Geodesy: ECEF, ENU

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
    dms2d(dms::AbstractString)::Float64

Parse `dms` from "dd:mm:ss.s" format to Float64 degrees.
"""
function dms2d(dms::AbstractString)::Float64
  d, m, s = map(x->parse(Float64,x), split(dms * ":0:0", ":"))
  sign = 1
  if '-' in dms
    sign = -1
    d = -d
  end
  d += m/60 + s/3600
  sign * d
end

dms2d(dms::Real)::Float64 = Float64(dms)
hms2h = dms2d


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
and W is in the direction of projection.

Left multiplying an ITRF-aligned (X,Y,Z) coordinate vector by the returned
rotation matrix will result in the (U,V,W) projection of the (X,Y,Z)
coordinates in the specified direction at the specified longitude:

    uvw = xyz2uvw(ha, dec, lon) * xyz
"""
function xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  [0 1 0
   0 0 1
   1 0 0] * ERFA.ry(-dec_rad,ERFA.rz(lon_rad-ha_rad))
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

"""
    enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert coordinates
in topocentric (E,N,Up) frame to a (U,V,W) projection with W pointing to the
hour angle `ha_rad` (west positive) and declination `dec_rad` as seen from
latitude `lat_rad` (all in radians).

This works by rotating the ENU frame anticlockwise about the E (i.e. first)
axis by `lat_rad`, producing a (E,Z,X') frame, then rotating that frame
anticlockwise about the Z (i.e. second) axis by `-ha_rad`, producing a
(U,Z,X") frame, then rotating anticlockwise about the U (i.e. first) axis by
`-dec_rad`, producing the (U,V,W) frame where U is east, V is north, and W is
in the direction of projection.

Left multiplying a topocentric (E,N,Up) coordinate vector by the returned
rotation matrix will result in the (U,V,W) projection of the (E,N,U)
coordinates in the specified direction at the specified latitude:

    uvw = enu2uvw(ha, dec, lat) * enu
"""
function enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,2}
   ERFA.rx(-dec_rad,ERFA.ry(-ha_rad,ERFA.rx(lat_rad)))
end

"""
    enu2uvw(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
            ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}

    enu2uvw(e::Real, n::Real, u::Real,
            ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}

Return the coordinataes of point(s) `enu` (or [e,n,u]) projected in the
direction of hour angle `ha_rad` (west positive) and `dec_rad` from latitude
`lat_rad`, all in radians, and permuted to a (U,V,W) frame where U is east, V
is north, and W is in the direction of projection.
"""
function enu2uvw(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
                 ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}
  enu2uvw(ha_rad, dec_rad, lat_rad) * enu
end

function enu2uvw(e::Real, n::Real, u::Real,
                 ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}
  enu2uvw(ha_rad, dec_rad, lat_rad) * [e, n, u]
end

end # module
