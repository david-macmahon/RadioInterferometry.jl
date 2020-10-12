"""
A module for radio interferometry that uses ERFA for much of its underlying
functionality.
"""
module RadioInterferometry

export I3
export ECEF_from_LLA
export LLA_from_ECEF
export enu2xyz
export xyz2enu
export xyz2uvw

using ERFA
import ERFA.DD2R
import ERFA.DR2D
import ERFA.WGS84

"""
`I3` is a 3x3 identity matrix.
"""
I3 = Float64[1 0 0
             0 1 0
             0 0 1]

"""
    xyz = ECEF_from_LLA(lat_deg::Real, lon_deg::Real, alt_m::Real,
                        ellipsoid::ERFA.Ellipsoid=ERFA.WGS84)::Array{Float64,1}

    xyz = ECEF_from_LLA(lla::Array{<:Real,1},
                        ellipsoid::ERFA.Ellipsoid=ERFA.WGS84)::Array{Float64,1}

Return Earth-Centered-Earth-Fixed geocentric ITRF coordinates for the given
latitude, longitude, and altitude using the specified reference ellipsoid.
`lat_deg` and `lon_deg` are the latitude and longitude, resp., in degress,
`alt_m` is the height in meters above the reference ellipsoid specified by
`ellipsoid`, which must be one of ERFA.WGS84, ERFA.GRS80, or ERFA.WGS72.
The `lla` array is `[lat_deg, lon_deg, alt_m]`.

The returned coordinates are in meters.  See `ERFA.gd2gc` for more details.
"""
function ECEF_from_LLA(lat_deg::Real, lon_deg::Real, alt_m::Real,
                       ellipsoid::ERFA.Ellipsoid=WGS84)::Array{Float64,1}
  ERFA.gd2gc(ellipsoid, lon_deg*DD2R, lat_deg*DD2R, alt_m)
end

function ECEF_from_LLA(lla::Array{<:Real,1},
                       ellipsoid::ERFA.Ellipsoid=WGS84)::Array{Float64,1}
  ECEF_from_LLA(lla..., ellipsoid)
end

"""
    lla = LLA_from_ECEF(x_m::Real, y_m::Real, z_m::Real,
                        ellipsoid::ERFA.Ellipsoid=ERFA.WGS84)::Array{Float64,1}

    lla = LLA_from_ECEF(xyz_m::Array{<:Real,1},
                        ellipsoid::ERFA.Ellipsoid=WGS84)::Array{Float64,1}

Return latitude, longitude, and altitude for the given
Earth-Centered-Earth-Fixed geocentric ITRF coordinates using the specified
reference ellipsoid.  `x_m`, `y_m`, and `z_m` are in meters.

The returned latitude and lingitude are in degrees.  The returned altidude is
in meters above the reference ellipsoid specified by `ellipsoid`, which must be
one of ERFA.WGS84 (the default), ERFA.GRS80, or ERFA.WGS72.  See `ERFA.gc2gd`
for more details.
"""
function LLA_from_ECEF(x_m::Real, y_m::Real, z_m::Real,
                       ellipsoid::ERFA.Ellipsoid=WGS84)::Array{Float64,1}
  ERFA.gc2gd(ellipsoid, x_m, y_m, z_m) .* [DR2D, DR2D, 1]
end

function LLA_from_ECEF(xyz_m::Array{<:Real,1},
                       ellipsoid::ERFA.Ellipsoid=WGS84)::Array{Float64,1}
  LLA_from_ECEF(xyz_m..., ellipsoid)
end

# Add default I3 argument to ERFA.rx, ERFA.ry, ERFA.rz
ERFA.rx(phi::Real) = ERFA.rx(phi, I3)
ERFA.ry(phi::Real) = ERFA.ry(phi, I3)
ERFA.rz(phi::Real) = ERFA.rz(phi, I3)

"""
    enu2xyz(lat_rad::Real, lon_rad::Real)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert topocentric
coordinates in east-north-up (E,N,U) frame with origin at latitude `lat_rad`
longitude `lon_rad` (both in radians) to topocentric coordinates in
ITRF-aligned (X,Y,Z) frame.  This works by rotating the ENU frame anticlockwise
about the E (i.e.  first) axis by `lat_rad`, producing a (E,Z,U') frame, then
rotating that frame anticlockwise about the Z (i.e. second) axis by `-lon_rad`,
producing a (Y,Z,X) frame which is then permuted to (X,Y,Z).  Left multiplying
an (E,N,U) coordinate vector by the returned rotation matrix will result in the
ITRF-aligned (X,Y,Z) coordinates:

    xyz = enu2xyz(lat, lon) * enu
"""
function enu2xyz(lat_rad::Real, lon_rad::Real)::Array{Float64,2}
  [0 0 1
   1 0 0
   0 1 0] * ERFA.ry(-lon_rad,ERFA.rx(lat_rad))
end

"""
    enu2xyz(enu::Array{<:Real,1},
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}

    enu2xyz(enu::Array{<:Real,2},
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,2}

    enu2xyz(e::Real, n::Real, u::Real,
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}

Return the coordinataes of point(s) `enu` transformed to the XYZ ITRF aligned
frame.
"""
function enu2xyz(enu::Array{<:Real,1},
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  enu2xyz(lat_rad, lon_rad) * enu
end

function enu2xyz(enu::Array{<:Real,2},
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  enu2xyz(lat_rad, lon_rad) * enu
end

function enu2xyz(e::Real, n::Real, u::Real,
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  enu2xyz(lat_rad, lon_rad) * [e, n, u]
end

"""
    xyz2enu(lat_rad::Real, lon_rad::Real)::Array{Float64,2}

Return transformation (rotation and permutation) matrix to convert topocentric
coordinates in ITRF aligned (X,Y,Z) frame to topocentric east-north-up (E,N,U)
frame with origin at latitude `lat_rad` longitude `lon_rad` (both in radians).
This works by rotating the XYZ frame anticlockwise about the Z (i.e. third)
axis by `lon_rad`, producing a (X',E,Z) frame, then rotating that frame
anticlockwise about the E (i.e. second) axis by `-lat_rad`, producing an
(U,E,N) frame which is then permuted to (E,N,U).  Left multiplying an
ITRF-aligned (X,Y,Z) coordinate vector by the returned rotation matrix will
result in the topocentric (E,N,U) coordinates relative to latitude `lat_rad`
and longitude `lon_rad`:

    enu = xyz2enu(lat, lon) * xyz
"""
function xyz2enu(lat_rad::Real, lon_rad::Real)::Array{Float64,2}
  [0 1 0
   0 0 1
   1 0 0] * ERFA.ry(-lat_rad,ERFA.rz(lon_rad))
end

"""
    xyz2enu(xyz::Array{<:Real,1},
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}

    xyz2enu(xyz::Array{<:Real,2},
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,2}

    xyz2enu(x::Real, y::Real, z::Real,
            lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}

Return the coordinataes of point(s) `xyz` (or [x,y,z]) transformed to the local
ENU frame at latitude `lat_rad` and longitude `lon_rad`, both in radians.
"""
function xyz2enu(xyz::Array{<:Real,1},
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2enu(lat_rad, lon_rad) * xyz
end

function xyz2enu(xyz::Array{<:Real,2},
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  xyz2enu(lat_rad, lon_rad) * xyz
end

function xyz2enu(x::Real, y::Real, z::Real,
                 lat_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2enu(lat_rad, lon_rad) * [x, y, z]
end

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
    xyz2uvw(xyz::Array{<:Real,1},
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}

    xyz2uvw(xyz::Array{<:Real,2},
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}

    xyz2uvw(x::Real, y::Real, z::Real,
            ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}

Return the coordinataes of point(s) `xyz` (or [x,y,z]) projected in the
direction of hour angle `ha_rad` (west positive) and `dec_rad` from longitude
`lon_rad`, all in radians, and permuted to a (U,V,W) frame where U is east, V
is north, and W is in the direction of projection.
"""
function xyz2uvw(xyz::Array{<:Real,1},
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * xyz
end

function xyz2uvw(xyz::Array{<:Real,2},
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,2}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * xyz
end

function xyz2uvw(x::Real, y::Real, z::Real,
                 ha_rad::Real, dec_rad::Real, lon_rad::Real=0)::Array{Float64,1}
  xyz2uvw(ha_rad, dec_rad, lon_rad) * [x, y, z]
end

end # module
