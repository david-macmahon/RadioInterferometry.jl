"""
A module for radio interferometry that uses ERFA and Rotations for much of its
underlying functionality.
"""
module RadioInterferometry

export dms2deg
export dms2rad
export dms2ha
export @dms_str
export hms2ha
export hms2deg
export hms2rad
export @hms_str

export deg2dms
export deg2dmsstr
export ha2dmsstr
export rad2dmsstr

export ha2hms
export ha2hmsstr
export deg2hmsstr
export rad2hmsstr

export deg2ha
export ha2deg
export rad2ha
export ha2rad

# Export RADec2Obs function
export radec2obs

export xyz2uvw
export xyz2enu
export enu2uvw
export enu2xyz

# TODO Allow these to reuse a user supplied Array for the rotation matrix
export xyz2uvw!
export xyz2enu!
export enu2uvw!
export enu2xyz!

using ERFA
import ERFA.DD2R
import ERFA.DR2D
import ERFA.CMPS

# Add some Julia (i.e. composable) implementations of ERFA utility functions.
include("ae2hd.jl")
include("hd2ae.jl")
include("hd2pa.jl")

# Helper function `radec2obs`
include("radec2obs.jl")
import .RADec2Obs: radec2obs

import LinearAlgebra: mul!
import Rotations: RotX, RotY, RotZ

# For Julia < 1.9.0
if !isdefined(Base, :get_extension)
using Requires
end
@static if !isdefined(Base, :get_extension)
# Add Symbolics functionality if/when it is imported.
function __init__()
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        include("../ext/RadioInterferometrySymbolicsExt.jl")
    end
end
end

# A method is added to `signclean` if/when the Symbolics package is loaded
function signclean end

"""
    rx(phi::Real)::AbstractArray{<:Real,2}
    ry(phi::Real)::AbstractArray{<:Real,2}
    rz(phi::Real)::AbstractArray{<:Real,2}

Return a 3x3 rotation matrix that will rotate a right handed reference frame by
`phi` radians anticlockwise around the X, Y, or Z axis.

# Example
```jldoctest
julia> RadioInterferometry.rx(π/2) * [0, 1, 0] ≈ [0, 0, -1]
true

julia> RadioInterferometry.ry(π/2) * [1, 0, 0] ≈ [0, 0, 1]
true

julia> RadioInterferometry.rz(π/2) * [1, 0, 0] ≈ [0, -1, 0]
true
```
"""
function rx(phi::Real)::AbstractArray{<:Real,2}
    RotX(-phi)
end,
function ry(phi::Real)::AbstractArray{<:Real,2}
    RotY(-phi)
end,
function rz(phi::Real)::AbstractArray{<:Real,2}
    RotZ(-phi)
end

# Some additional ways of expressing speed of light
"""
Speed of light in a vacuum (C) in meters/nanosecond
"""
CMPNS = CMPS / 1e9
"""
Speed of light in a vacuum (C) in nanoseconds/meter
"""
CNSPM = 1 / CMPNS
"""
Speed of light in a vacuum (C) in seconds/meter
"""
CSPM  = 1 / CMPS

"""
    dms2deg(dms::AbstractString)::Float64
    dms2rad(dms::AbstractString)::Float64
    dms2ha(dms::AbstractString)::Float64
    dms2deg(d::Real)::Float64
    dms2rad(d::Real)::Float64
    dms2ha(d::Real)::Float64

Parse `dms::AbstractString` from "dd:mm:ss.s" format to `Float64`
degrees/radians/hour angle or convert `d::Real` to `Float64`
degrees/radians/hour angle.
"""
function dms2deg end,
function dms2ha end,
function dms2rad end

function dms2deg(dms::AbstractString)::Float64
    d, m, s = map(x->parse(Float64,x), split(dms * ":0:0", ":"))
    sign = 1
    if '-' in dms
        sign = -1
        d = -d
    end
    d += m/60 + s/3600
    sign * d
end

@deprecate dms2d dms2deg

dms2rad(dms::AbstractString) = deg2rad(dms2deg(dms))
dms2ha(dms::AbstractString) = dms2deg(dms)/15
dms2deg(d::Real)::Float64 = Float64(d)
dms2rad(d::Real)::Float64 = Float64(deg2rad(d))
dms2ha(d::Real)::Float64 = Float64(d)/15

"""
    @dms_str("dd:mm:ss.s", units="deg")::Float64
    dms"dd:mm:ss.s"deg (or dms"dd:mm:ss.s")
    dms"dd:mm:ss.s"ha
    dms"dd:mm:ss.s"rad

Convert sexagesimal degrees:minutes:seconds string to decimal degrees,
optionally scaled to `units`. `units` can be `"deg"` for degrees (default),
`"ha"` for hour angle, or `"rad"` for radians. Note that passing `"ha"` for
units converts to degrees first, then scales to hour angle.
"""
macro dms_str(s, u="deg")
    @assert u in ("deg", "ha", "rad")
    d = dms2deg(s)

    u == "deg" ? dms2deg(s) :
    u == "rad" ? dms2rad(s) : dms2ha(d)
end

"""
    hms2ha(hms::AbstractString)::Float64
    hms2ha(h::Real)::Float64
    hms2deg(hms::AbstractString)::Float64
    hms2deg(h::Real)::Float64
    hms2rad(hms::AbstractString)::Float64
    hms2rad(h::Real)::Float64

Parse `hms::AbstractString` from "hh:mm:ss.s" format to Float64 hour
angle/degrees/radians or convert `h::Real` to Float64 hour
angle/degrees/radians.
"""
function hms2ha end,
function hms2deg end,
function hms2rad end

hms2ha(h) = dms2deg(h)
hms2deg(h) = 15*hms2ha(h)
hms2rad(h) = deg2rad(hms2deg(h))

@deprecate hms2h hms2ha


"""
    @hms_str("hh:mm:ss.s", units="deg")::Float64
    hms"hh:mm:ss.s"ha (or hms"hh:mm:ss.s")
    hms"hh:mm:ss.s"deg
    hms"hh:mm:ss.s"rad

Convert sexagesimal hours:minutes:seconds string to decimal hours, optionally
scaled to `units`. `units` can be `"ha"` for hour angle (default), `"deg"`
for degrees, or `"rad"` for radians. Note that passing `"deg"` for units
converts to hour angle first, then scales to degrees.
"""
macro hms_str(s, u="ha")
    @assert u in ("deg", "ha", "rad")

    u == "ha"  ? hms2ha(s)  :
    u == "rad" ? hms2rad(s) : hms2deg(s)
end

"""
    deg2dms(d::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}

Convert degrees `d` into `(sign, degrees, arcminutes, arcseconds, fraction)`
using `ERFA.a2af()`.

`sign` is `-1` if `d < 0`, otherwise `+1`.  `ndp` is the number of decimal
places of precision.  For `ndp > 0`, the denominator of `fraction` will be
`10^ndp`.  For `ndp <= 0`, `fraction` will be `0//1` (i.e. zero), and the
returned resolution of `degrees`, `arcminutes`, `arcseconds` will be limited
accordingly.  For example, `ndp = -3` will return results rounded to the
nearest multiple of 10 arcminutes.  See `ERFA.a2af` for more details.
"""
function deg2dms(d::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}
    ndp <= 9 || @warn "fraction may be inaccurate when ndp ($ndp) > 9"
    sign, deg, min, sec, frac = ERFA.a2af(ndp, d*ERFA.DD2R)
    (Int32(sign == '-' ? -1 : +1), deg, min, sec, ndp > 0 ? frac//10^ndp : 0//1)
end

@deprecate d2dms deg2dms

"""
    deg2dmsstr(d::Real, ndp::Integer=3; <kwargs>)::String

Return a sexagesimal string representing `d` degrees, rounded to the precision
specified by `ndp`.

# Keyword arguments
- `posind::AbstractString="+"`: string to indicate nonnegative values.
- `degwidth::Integer=0`: minimum width of the degree field.
- `degpad::Union{AbstractChar,AbstractString}='0'`: padding for degree field.

See also: [`deg2dms`](@ref)
"""
function deg2dmsstr(d::Real, ndp::Integer=3;
                    posind::AbstractString="+",
                    degwidth::Integer=0,
                    degpad::Union{AbstractChar,AbstractString}='0'
                   )::String
    sign, deg, min, sec, frac = deg2dms(d,ndp)
    signstr = sign < 0 ? "-" : posind
    degstr = lpad(deg, degwidth, degpad)
    minstr = lpad(min, 2, '0')
    secstr = lpad(sec, 2, '0')
    fracstr = ndp > 0 ? ".$(lpad(Int(round(frac*10^ndp)), ndp, '0'))" : ""
    "$(signstr)$(degstr):$(minstr):$(secstr)$(fracstr)"
end

@deprecate d2dmsstr deg2dmsstr

"""
    ha2dmsstr(h::Real, ndp::Integer=3; kwargs...)::String

Return `deg2dmsstr(ha2deg(h), ndp; kwargs...)`
"""
function ha2dmsstr(h::Real, ndp::Integer=3; kwargs...)
    deg2dmsstr(ha2deg(h), ndp; kwargs...)
end

"""
    rad2dmsstr(r::Real, ndp::Integer=3; kwargs...)::String

Return `deg2dmsstr(rad2deg(r), ndp; kwargs...)`
"""
function rad2dmsstr(r::Real, ndp::Integer=3; kwargs...)
    deg2dmsstr(rad2deg(r), ndp; kwargs...)
end

"""
    deg2dmsstr(;digits=3, kwargs...)
    deg2hmsstr(;digits=3, kwargs...)
    ha2dmsstr(;digits=3, kwargs...)
    ha2hmsstr(;digits=3, kwargs...)
    rad2dmsstr(;digits=3, kwargs...)
    rad2hmsstr(;digits=3, kwargs...)

Return 1-input function that passes its argument as the first parameter and
`digits` as the second parameter to the 2-arg method of the same function.
This can be useful with function chaining/piping.  `kwargs` are passed
through too.

# Examples
```jldoctest
julia> 3661//3600 |> deg2dmsstr
"+1:01:01.000"

julia> 34.567 |> deg2dmsstr
"+34:34:01.200"

julia> -0.261799315 |> rad2hmsstr
"-0:59:59.999"

julia> -0.261799315 |> rad2hmsstr(digits=2)
"-1:00:00.00"

julia> 3+sqrt(2)/10 |> rad2dmsstr(digits=6)
"+179:59:24.667385"
```
"""
function deg2dmsstr(;digits=3, kwargs...)
    x->deg2dmsstr(x, digits; kwargs...)
end,
function deg2hmsstr(;digits=3, kwargs...)
    x->deg2hmsstr(x, digits; kwargs...)
end,
function ha2dmsstr( ;digits=3, kwargs...)
    x->ha2dmsstr( x, digits; kwargs...)
end,
function ha2hmsstr( ;digits=3, kwargs...)
    x->ha2hmsstr( x, digits; kwargs...)
end,
function rad2dmsstr(;digits=3, kwargs...)
    x->rad2dmsstr(x, digits; kwargs...)
end,
function rad2hmsstr(;digits=3, kwargs...)
    x->rad2hmsstr(x, digits; kwargs...)
end

"""
    ha2hms(h::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}

Convert hours `h` into `(sign, hours, minutes, seconds, fraction)`
using `ERFA.a2tf()`.

`sign` is `-1` if `d < 0`, otherwise `+1`.  `ndp` is the number of decimal
places of precision.  For `ndp > 0`, the denominator of `fraction` will be
`10^ndp`.  For `ndp <= 0`, `fraction` will be `0//1` (i.e. zero), and the
returned resolution of `degrees`, `minutes`, `seconds` will be limited
accordingly.  For example, `ndp = -3` will return results rounded to the
nearest multiple of 10 minutes.  See `ERFA.a2tf` for more details.
"""
function ha2hms(h::Real, ndp::Integer=3)::Tuple{Int32, Int32, Int32, Int32, Rational{Int64}}
    ndp <= 9 || @warn "fraction may be inaccurate when ndp ($ndp) > 9"
    sign, deg, min, sec, frac = ERFA.a2tf(ndp, h*15*ERFA.DD2R)
    (Int32(sign == '-' ? -1 : +1), deg, min, sec, ndp > 0 ? frac//10^ndp : 0//1)
end

@deprecate h2hms ha2hms

"""
    ha2hmsstr(h::Real, ndp::Integer=3; <kwargs>)::String

Return a sexagesimal string representing `h` hours, rounded to the precision
specified by `ndp`.

# Keyword arguments
- `posind::AbstractString=""`: string to indicate nonnegative values.
- `hourwidth::Integer=0`: minimum width of the hour field.
- `hourpad::Union{AbstractChar,AbstractString}='0'`: padding for hour field.

See also: [`hhdms`](@ref)
"""
function ha2hmsstr(h::Real, ndp::Integer=3;
                   posind::AbstractString="",
                   hourwidth::Integer=0,
                   hourpad::Union{AbstractChar,AbstractString}='0'
                  )::String
    sign, hour, min, sec, frac = ha2hms(h,ndp)
    signstr = sign < 0 ? "-" : posind
    hourstr = lpad(hour, hourwidth, hourpad)
    minstr = lpad(min, 2, '0')
    secstr = lpad(sec, 2, '0')
    fracstr = ndp > 0 ? ".$(lpad(Int(round(frac*10^ndp)), ndp, '0'))" : ""
    "$(signstr)$(hourstr):$(minstr):$(secstr)$(fracstr)"
end

@deprecate h2hmsstr ha2hmsstr

"""
    deg2hmsstr(d::Real, ndp::Integer=3; kwargs...)::String

Return `ha2hmsstr(deg2ha(d), ndp; kwargs...)`
"""
function deg2hmsstr(d::Real, ndp::Integer=3; kwargs...)
    ha2hmsstr(deg2ha(d), ndp; kwargs...)
end

"""
    rad2hmsstr(r::Real, ndp::Integer=3; kwargs...)::String

Return `ha2hmsstr(rad2ha(d), ndp; kwargs...)`
"""
function rad2hmsstr(r::Real, ndp::Integer=3; kwargs...)
    ha2hmsstr(rad2ha(r), ndp; kwargs...)
end

"""
    deg2ha(d)
Convert `d` from degrees to hour angle.
"""
deg2ha(d) = d/15

"""
    ha2deg(h)
Convert `h` from hour angle to degrees.
"""
ha2deg(h) = h*15

"""
    rad2ha(r)
Convert `r` from radians to hour angle.
"""
rad2ha = deg2ha ∘ rad2deg

"""
    ha2rad(h)
Convert `h` from hour angle to radians.
"""
ha2rad = deg2rad ∘ ha2deg

"""
    xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}

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
function xyz2uvw(ha_rad::Real, dec_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}
    [0 1 0
     0 0 1
     1 0 0] * ry(-dec_rad) * rz(lon_rad-ha_rad)
end

"""
    xyz2uvw(xyz::AbstractArray{<:Real},
            ha_rad::Real, dec_rad::Real, lon_rad::Real
           )::AbstractArray{<:Real}

Transform point(s) `xyz` from an (X,Y,Z) ITRF aligned frame to a (U,V,W) frame
where U is eastward, V is northward, and W points to the hour angle `ha_rad`
(west positive) and declination `dec_rad` as seen from longitude `lon_rad` (all
in radians):

    uvw = xyz2uvw(xyz, ha, dec, lon)
"""
function xyz2uvw(xyz::AbstractArray{<:Real},
                 ha_rad::Real, dec_rad::Real, lon_rad::Real
                )::AbstractArray{<:Real}
    xyz2uvw(ha_rad, dec_rad, lon_rad) * xyz
end

"""
    xyz2uvw!(uvw::AbstractArray{T},
             xyz::AbstractArray{<:Real},
             ha_rad::Real, dec_rad::Real, lon_rad::Real
            )::AbstractArray{T} where {T<:Real}

Transform point(s) `xyz` from an (X,Y,Z) ITRF aligned frame to a (U,V,W) frame
where U is eastward, V is northward, and W points to the hour angle `ha_rad`
(west positive) and declination `dec_rad` as seen from longitude `lon_rad` (all
in radians) and store result in in `uvw`:

    xyz2uvw!(uvw, xyz, ha, dec, lat)
"""
function xyz2uvw!(uvw::AbstractArray{T},
                  xyz::AbstractArray{<:Real},
                  ha_rad::Real, dec_rad::Real, lon_rad::Real
                 )::AbstractArray{T} where {T<:Real}
    mul!(uvw, xyz2uvw(ha_rad, dec_rad, lon_rad), xyz)
end

"""
    xyz2enu(lat_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}

Return transformation (rotation and permutation) matrix to convert
coordinates from a topocentric ITRF aligned (X,Y,Z) frame to a topocentric
(East,North,Up) frame for topocentric origin at geodetic latitude `lat_rad`
and longitude `lon_rad` (both in radians).

This works by rotating the XYZ frame anticlockwise about the Z (i.e. third)
axis by `lon_rad`, producing a (X',East,Z) frame, then rotating that frame
anticlockwise about the E (i.e. second) axis by `-lat_rad`, producing a
(U,E,N) frame which is then permuted to (E,N,U).

Left multiplying a topocentric (X,Y,Z) coordinate vector or matrix by the
returned transformation matrix will result in the topocentric (East,North,Up)
coordinate(s) for the given latitude and longitude:

    enu = xyz2enu(lat, lon) * xyz
"""
function xyz2enu(lat_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}
    [0 1 0
     0 0 1
     1 0 0] * ry(-lat_rad) * rz(lon_rad)
end

"""
    xyz2enu(xyz::AbstractArray{<:Real},
            lat_rad::Real, lon_rad::Real
           )::AbstractArray{<:Real}

Transform point(s) `xyz` from a topocentric ITRF aligned (X,Y,Z) frame to a
topocentric (East,North,Up) frame for topocentric origin at geodetic latitude
`lat_rad` and longitude `lon_rad` (both in radians):

    enu = xyz2enu(xyz, lat, lon)
"""
function xyz2enu(xyz::AbstractArray{<:Real},
                 lat_rad::Real, lon_rad::Real
                )::AbstractArray{<:Real}
    xyz2enu(lat_rad, lon_rad) * xyz
end

"""
    xyz2enu!(enu::AbstractArray{T},
             xyz::AbstractArray{<:Real},
             lat_rad::Real, lon_rad::Real
            )::AbstractArray{T} where {T<:Real}

Transform point(s) `xyz` from a topocentric ITRF aligned (X,Y,Z) frame to a
topocentric (East,North,Up) frame for topocentric origin at geodetic latitude
`lat_rad` and longitude `lon_rad` (both in radians) and store result in
`enu`.

    xyz2enu!(enu, xyz, lat, lon)
"""
function xyz2enu!(enu::AbstractArray{T},
                  xyz::AbstractArray{<:Real},
                  lat_rad::Real, lon_rad::Real
                 )::AbstractArray{T} where {T<:Real}
    mul!(enu, xyz2enu(lat_rad, lon_rad), xyz)
end

"""
    enu2uvw(az_rad::Real, el_rad::Real)::AbstractArray{<:Real,2}

Return rotation matrix to convert coordinates from a topocentric
(East,North,Up) frame to a (U,V,W) frame where U is eastward, V is northward,
and W points to azimuth `az_rad` (east positive) and elevation `el_rad` (both
in radians).

This works by rotating the ENU frame anticlockwise about the Up (i.e. third)
axis by `-az_rad`, producing a (W',U,Up) frame, then rotating that frame
anticlockwise about the U (i.e. second) axis by `-el_rad`, producing the
(W,U,V) frame, which is then permuted to (U,V,W).

Left multiplying a topocentric (East,North,Up) coordinate vector or matrix by
the returned rotation matrix will result in the corresponding (U,V,W)
coordinates for the given direction:

# Examples
```jldoctest
# East-only vector projected due east at the horizon is all `w`.
julia> enu2uvw(π/2, 0) * [1, 0, 0] ≈ [0, 0, 1]
true
```
```jldoctest
# North-only vector projected due south at the horizon is all `-w`.
julia> enu2uvw(π, 0) * [0, 2, 0] ≈ [0, 0, -2]
true
```
"""
function enu2uvw(az_rad::Real, el_rad::Real)::AbstractArray{<:Real,2}
    [0 1 0
     0 0 1
     1 0 0] * ry(-el_rad) * rz(π/2-az_rad)
end

"""
    enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real)::AbstractArray{<:Real,2}

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
# Examples
```jldoctest
# East-only vector projected due east at the horizon is all `w`.
julia> enu2uvw(-π/2, 0, 0) * [1, 0, 0] ≈ [0, 0, 1]
true
```
```jldoctest
# North-only vector projected due south at the horizon is all `-w`.
# This latitude is south of the equator and the hour angle is 12 hours.
julia> enu2uvw(π, -π/4, -π/4) * [0, 2, 0] ≈ [0, 0, -2]
true
```
"""
function enu2uvw(ha_rad::Real, dec_rad::Real, lat_rad::Real)::AbstractArray{<:Real,2}
    rx(-dec_rad) * ry(-ha_rad) * rx(lat_rad)
end

"""
    enu2uvw(enu::AbstractArray{<:Real},
            ha_rad::Real, dec_rad::Real, lat_rad::Real
           )::AbstractArray{<:Real}

Transform point(s) `enu` from a topocentric (East,North,Up) frame to a (U,V,W)
frame where U is eastward, V is northward, and W points to the hour angle
`ha_rad` (west positive) and declination `dec_rad` as seen from latitude
`lat_rad` (all in radians):

    uvw = enu2uvw(enu, ha, dec, lat)
"""
function enu2uvw(enu::AbstractArray{<:Real},
                 ha_rad::Real, dec_rad::Real, lat_rad::Real
                )::AbstractArray{<:Real}
    enu2uvw(ha_rad, dec_rad, lat_rad) * enu
end

"""
    enu2uvw!(uvw::AbstractArray{T}
             enu::AbstractArray{<:Real},
             ha_rad::Real, dec_rad::Real, lat_rad::Real
            )::AbstractArray{T} where {T<:Real}

Transform point(s) `enu` from a topocentric (East,North,Up) frame to a (U,V,W)
frame where U is eastward, V is northward, and W points to the hour angle
`ha_rad` (west positive) and declination `dec_rad` as seen from latitude
`lat_rad` (all in radians) and store result in `uvw`:

    enu2uvw!(uvw, enu, ha, dec, lat)
"""
function enu2uvw!(uvw::AbstractArray{T},
                  enu::AbstractArray{<:Real},
                  ha_rad::Real, dec_rad::Real, lat_rad::Real
                 )::AbstractArray{T} where {T<:Real}
    mul!(uvw, enu2uvw(ha_rad, dec_rad, lat_rad), enu)
end

"""
    enu2xyz(lat_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}

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
function enu2xyz(lat_rad::Real, lon_rad::Real)::AbstractArray{<:Real,2}
    [0 0 1
     1 0 0
     0 1 0] * ry(-lon_rad) * rx(lat_rad)
end

"""
    enu2xyz(enu::AbstractArray{<:Real},
            lat_rad::Real, lon_rad::Real
           )::AbstractArray{<:Real}

Transform point(s) `enu` from a topocentric (East,North,Up) frame to a
topocentric ITRF aligned (X,Y,Z) frame for topocentric origin at geodetic
latitude `lat_rad` and longitude `lon_rad` (both in radians):

    xyz = enu2xyz(enu, lat, lon)
"""
function enu2xyz(enu::AbstractArray{<:Real},
                 lat_rad::Real, lon_rad::Real
                )::AbstractArray{<:Real}
    enu2xyz(lat_rad, lon_rad) * enu
end

"""
    enu2xyz!(xyz::AbstractArray{T},
             enu::AbstractArray{<:Real},
             lat_rad::Real, lon_rad::Real
            )::AbstractArray{T} where {T<:Real}

Transform point(s) `enu` from a topocentric (East,North,Up) frame to a
topocentric ITRF aligned (X,Y,Z) frame for topocentric origin at geodetic
latitude `lat_rad` and longitude `lon_rad` (both in radians) and store result
in `xyz`.

    enu2xyz!(xyz, enu, lat, lon)
"""
function enu2xyz!(xyz::AbstractArray{T},
                  enu::AbstractArray{<:Real},
                  lat_rad::Real, lon_rad::Real
                 )::AbstractArray{T} where {T<:Real}
    mul!(xyz, enu2xyz(lat_rad, lon_rad), enu)
end

end # module
