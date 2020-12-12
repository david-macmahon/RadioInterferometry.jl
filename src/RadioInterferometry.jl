"""
A module for radio interferometry that uses ERFA for much of its underlying
functionality.
"""
module RadioInterferometry

export dms2d
export @dms_str
export hms2h
export @hms_str

export d2dms
export d2dmsstr

export h2hms
export h2hmsstr

export xyz2uvw
export enu2uvw
export enu2xyz

using ERFA
import ERFA.DD2R
import ERFA.DR2D
import ERFA.CMPS

import Geodesy: ECEF, ENU

"""
`I3` is a 3x3 identity matrix.
"""
const I3 = Float64[1 0 0
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
    dms2d(d::Real)::Float64

Parse `dms::AbstractString` from "dd:mm:ss.s" format to `Float64` degrees or
convert `d::Real` to `Float64`.
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

dms2d(d::Real)::Float64 = Float64(d)

"""
    @dms_str("dd:mm:ss.s")::Float64
    dms"dd:mm:ss.s"

Convert sexagesimal degrees:minutes:seconds string to decimal degrees.
"""
macro dms_str(s)
  dms2d(s)
end

"""
    hms2h(hms::AbstractString)::Float64
    hms2h(h::Real)::Float64

Parse `hms::AbstractString` from "hh:mm:ss.s" format to Float64 hours or
convert `h::Real` to Float64.
"""
hms2h = dms2d

"""
    @hms_str("hh:mm:ss.s")
    hms"hh:mm:ss.s"

Convert sexagesimal hours:minutes:seconds string to decimal hours.
"""
macro hms_str(s)
  hms2h(s)
end

"""
    d2dms(d::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Rational{Int64}}

Convert degrees `d` into `(sign, degrees, arcminutes, arcseconds, fraction)`
using `ERFA.a2af()`.

`sign` is `-1` if `d < 0`, otherwise `+1`.  `ndp` is the number of decimal
places of precision.  For `ndp > 0`, the denominator of `fraction` will be
`10^ndp`.  For `ndp <= 0`, `fraction` will be `0//1` (i.e. zero), and the
returned resolution of `degrees`, `arcminutes`, `arcseconds` will be limited
accordingly.  For example, `ndp = -3` will return results rounded to the
nearest multiple of 10 arcminutes.  See `ERFA.a2af` for more details.
"""
function d2dms(d::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}
  ndp <= 9 || @warn "fraction may be inaccurate when ndp ($ndp) > 9"
  sign, deg, min, sec, frac = ERFA.a2af(ndp, d*ERFA.DD2R)
  (Int32(sign == '-' ? -1 : +1), deg, min, sec, ndp > 0 ? frac//10^ndp : 0//1)
end

"""
    d2dmsstr(d::Real, ndp::Integer=3; <kwargs>)::String

Return a sexagesimal string representing `d` degrees, rounded to the precision
specified by `ndp`.

# Keyword arguments
- `posind::AbstractString="+"`: string to indicate nonnegative values.
- `degwidth::Integer=0`: minimum width of the degree field.
- `degpad::Union{AbstractChar,AbstractString}='0'`: padding for degree field.

See also: [`d2dms`](@ref)
"""
function d2dmsstr(d::Real, ndp::Integer=3;
                  posind::AbstractString="+",
                  degwidth::Integer=0,
                  degpad::Union{AbstractChar,AbstractString}='0'
                 )::String
  sign, deg, min, sec, frac = d2dms(d,ndp)
  signstr = sign < 0 ? "-" : posind
  degstr = lpad(deg, degwidth, degpad)
  minstr = lpad(min, 2, '0')
  secstr = lpad(sec, 2, '0')
  fracstr = ndp > 0 ? ".$(lpad(Int(round(frac*10^ndp)), ndp, '0'))" : ""
  "$(signstr)$(degstr):$(minstr):$(secstr)$(fracstr)"
end

"""
    h2hms(h::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Rational{Int64}}

Convert hours `h` into `(sign, hours, minutes, seconds, fraction)`
using `ERFA.a2tf()`.

`sign` is `-1` if `d < 0`, otherwise `+1`.  `ndp` is the number of decimal
places of precision.  For `ndp > 0`, the denominator of `fraction` will be
`10^ndp`.  For `ndp <= 0`, `fraction` will be `0//1` (i.e. zero), and the
returned resolution of `degrees`, `minutes`, `seconds` will be limited
accordingly.  For example, `ndp = -3` will return results rounded to the
nearest multiple of 10 minutes.  See `ERFA.a2tf` for more details.
"""
function h2hms(h::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}
  ndp <= 9 || @warn "fraction may be inaccurate when ndp ($ndp) > 9"
  sign, deg, min, sec, frac = ERFA.a2tf(ndp, h*15*ERFA.DD2R)
  (Int32(sign == '-' ? -1 : +1), deg, min, sec, ndp > 0 ? frac//10^ndp : 0//1)
end

"""
    h2hmsstr(h::Real, ndp::Integer=3; <kwargs>)::String

Return a sexagesimal string representing `h` hours, rounded to the precision
specified by `ndp`.

# Keyword arguments
- `posind::AbstractString=""`: string to indicate nonnegative values.
- `hourwidth::Integer=0`: minimum width of the hour field.
- `hourpad::Union{AbstractChar,AbstractString}='0'`: padding for hour field.

See also: [`hhdms`](@ref)
"""
function h2hmsstr(h::Real, ndp::Integer=3;
                  posind::AbstractString="",
                  hourwidth::Integer=0,
                  hourpad::Union{AbstractChar,AbstractString}='0'
                 )::String
  sign, hour, min, sec, frac = h2hms(h,ndp)
  signstr = sign < 0 ? "-" : posind
  hourstr = lpad(hour, hourwidth, hourpad)
  minstr = lpad(min, 2, '0')
  secstr = lpad(sec, 2, '0')
  fracstr = ndp > 0 ? ".$(lpad(Int(round(frac*10^ndp)), ndp, '0'))" : ""
  "$(signstr)$(hourstr):$(minstr):$(secstr)$(fracstr)"
end

"""
    xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert coordinates
from an ITRF aligned (X,Y,Z) frame to a (U,V,W) frame where U is eastward, V is
northward, and W points to the hour angle `ha_rad` (west positive) and
declination `dec_rad` as seen from longitude `lon_rad` (all in radians).

This works by rotating the XYZ frame anticlockwise about the Z (i.e. third)
axis by `lon_rad-ha_rad`, producing a (X',U,Z) frame, then rotating that frame
anticlockwise about the U (i.e. second) axis by `-dec_rad`, producing an
(W,U,V) frame which is then permuted to (U,V,W) where U is east, V is north,
and W is in the direction of the given hour angle and declination as seen from
the given longitude.

Left multiplying an ITRF-aligned (X,Y,Z) coordinate vector or matrix by the
returned rotation matrix will result in the corresponding (U,V,W) coordinates
for the given direction and longitude:

    uvw = xyz2uvw(ha, dec, lon) * xyz
"""
function xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  [0 1 0
   0 0 1
   1 0 0] * ERFA.ry(-dec_rad,ERFA.rz(lon_rad-ha_rad))
end

"""
    xyz2uvw(xyz::Union{Array{<:Real,1},Array{<:Real,2},ECEF},
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0
           )::Union{Array{Float64,1},Array{Float64,2}}

    xyz2uvw(x::Real, y::Real, z::Real,
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}

Transform point(s) `xyz` (or [x,y,z]) from an (X,Y,Z) ITRF aligned frame to a
(U,V,W) frame where U is eastward, V is northward, and W points to the hour
angle `ha_rad` (west positive) and declination `dec_rad` as seen from longitude
`lon_rad` (all in radians):

    uvw = xyz2uvw(xyz, ha, dec, lon)
"""
function xyz2uvw(xyz::Union{Array{<:Real,1},Array{<:Real,2},ECEF},
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0
                )::Union{Array{Float64,1},Array{Float64,2}}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * xyz
end

function xyz2uvw(x::Real, y::Real, z::Real,
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * [x, y, z]
end

"""
    enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,2}

Return rotation matrix to convert coordinates from a topocentric
(East,North,Up) frame to a (U,V,W) frame where U is eastward, V is northward,
and W points to the hour angle `ha_rad` (west positive) and declination
`dec_rad` as seen from latitude `lat_rad` (all in radians).

This works by rotating the ENU frame anticlockwise about the East (i.e. first)
axis by `lat_rad`, producing a (East,Z,X') frame, then rotating that frame
anticlockwise about the Z (i.e. second) axis by `-ha_rad`, producing a
(U,Z,X") frame, then rotating anticlockwise about the U (i.e. first) axis by
`-dec_rad`, producing the (U,V,W) frame where U is east, V is north, and W is
in the direction of projection.

Left multiplying a topocentric (East,North,Up) coordinate vector or matrix by
the returned rotation matrix will result in the corresponding (U,V,W)
coordinates for the given direction and latitude:

    uvw = enu2uvw(ha, dec, lat) * enu
"""
function enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,2}
   ERFA.rx(-dec_rad,ERFA.ry(-ha_rad,ERFA.rx(lat_rad)))
end

"""
    enu2uvw(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
            ha_rad::Real, dec_rad::Real, lat_rad::Real=0
           )::Union{Array{Float64,1},Array{Float64,2}}

    enu2uvw(e::Real, n::Real, u::Real,
            ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}

Transform point(s) `enu` (or `[e,n,u]`) from a topocentric (East,North,Up)
frame to a (U,V,W) frame where U is eastward, V is northward, and W points to
the hour angle `ha_rad` (west positive) and declination `dec_rad` as seen from
latitude `lat_rad` (all in radians):

    uvw = enu2uvw(enu, ha, dec, lat)
"""
function enu2uvw(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
                 ha_rad::Real, dec_rad::Real, lat_rad::Real=0
                )::Union{Array{Float64,1},Array{Float64,2}}
  enu2uvw(ha_rad, dec_rad, lat_rad) * enu
end

function enu2uvw(e::Real, n::Real, u::Real,
                 ha_rad::Real, dec_rad::Real, lat_rad::Real=0)::Array{Float64,1}
  enu2uvw(ha_rad, dec_rad, lat_rad) * [e, n, u]
end

"""
    enu2xyz(lat_rad::Real=0, lon_rad::Real=0)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert coordinates
from a topocentric (East,North,Up) frame to a topocentric ITRF aligned (X,Y,Z)
frame for topocentric origin at geodetic latitude `lat_rad` and longitude
`lon_rad` (both in radians).

This works by rotating the ENU frame anticlockwise about the East (i.e. first)
axis by `lat_rad`, producing a (East,Z,X') frame, then rotating that frame
anticlockwise about the Z (i.e. second) axis by `-lon_rad`, producing a
(Y,Z,X) frame which is then permuted to (X,Y,Z).

Left multiplying a topocentric (East,North,Up) coordinate vector or matrix by
the returned transformation matrix will result in the topocentric (X,Y,Z)
coordinate(s) for the given latitude and longitude:

    xyz = enu2xyz(lat, lon) * enu
"""
function enu2xyz(lat_rad::Real=0, lon_rad::Real=0)::Array{Float64,2}
   [0 0 1
    1 0 0
    0 1 0] * ERFA.ry(-lon_rad,ERFA.rx(lat_rad))
end

"""
    enu2xyz(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
            lat_rad::Real=0, lon_rad::Real=0
           )::Union{Array{Float64,1},Array{Float64,2}}

    enu2uvw(e::Real, n::Real, u::Real,
            lat_rad::Real=0, lon_rad::Real=0)::Array{Float64,1}

Transform point(s) `enu` or `[e,n,u]` from a topocentric (East,North,Up) frame
to a topocentric ITRF aligned (X,Y,Z) frame for topocentric origin at geodetic
latitude `lat_rad` and longitude `lon_rad` (both in radians):

    xyz = enu2xyz(enu, lat, lon)
"""
function enu2xyz(enu::Union{Array{<:Real,1},Array{<:Real,2},ENU},
                 lat_rad::Real=0, lon_rad::Real=0
                )::Union{Array{Float64,1},Array{Float64,2}}
  enu2xyz(lat_rad, lon_rad) * enu
end

function enu2xyz(e::Real, n::Real, u::Real,
                 lat_rad::Real=0, lon_rad::Real=0)::Array{Float64,1}
  enu2xyz(lat_rad, lon_rad) * [e, n, u]
end

end # module
