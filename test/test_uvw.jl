include("../src/RadioInterferometry.jl")
using .RadioInterferometry
using ERFA
using Geodesy
using EarthOrientation

import ERFA.rx, ERFA.ry, ERFA.rz

function test_uvw()

# Test parameters and expected results adapted from pyuvdata's "test_utils.py"
# script, specifically the "test_phasing_funcs()" function.  From comments in
# that function:
#
#   > these tests are based on a notebook where I tested against the mwa_tools
#   > phasing code

### Test Parameters

# RA, Dec, and epoch (CIRS?)
ra_deg = 12.1 * 15
dec_deg = -42.3
mjd = 55780.1

# Array center in ITRF XYZ (m)
array_center_xyz = [-2559454.08, 5095372.14, -2849057.18]

# Antenna location in topocentric east/north/up (m)
antenna_location_enu = [-101.94, 156.41, 1.24]

### Expected results

u_exp = -97.122828
v_exp = 50.388281
w_exp = -151.27976

### Intermediate transforms/calculations/data (NOT from pyuvdata)

# Create ECEF and LLA objects from array_center_xyz (using WGS84 ellipsoid)
ref_ecef = ECEF(array_center_xyz)
ref_lla = LLA(ref_ecef, wgs84)
@show ref_ecef
@show ref_lla


# Create ENU and ECEF objects from ant_enu (relative to ref_ecef, using WGS84)
ant_enu = ENU(antenna_location_enu)
ant_ecef = ECEF(ant_enu, ref_ecef, wgs84)

# DUT1 at given MJD
dut1 = getΔUT1(ERFA.DJM0 + mjd)

# Polar motion at given MJD
xp, yp = polarmotion(ERFA.DJM0 + mjd) .* ERFA.DMAS2R

# Compute various "observed" coordinates when viewing ra/dec (CIRS?) from array
# center at given MJD. aob is azimuth, zof is zenitt distance, hob is hour
# angle, dob is declination, rob is right ascension.
#aob, zob, hob, dob, rob = ERFA.atio13(ra_deg*ERFA.DD2R, dec_deg*ERFA.DD2R,
aob, zob, hob, dob, rob = ERFA.atco13(ra_deg*ERFA.DD2R, dec_deg*ERFA.DD2R, 0, 0, 0, 0,
                                      ERFA.DJM0, mjd, dut1,
                                      ref_lla.lon*ERFA.DD2R,
                                      0*ref_lla.lat*ERFA.DD2R,
                                      ref_lla.alt,
                                      xp, yp,
                                      0, 0, 0, 0) # No refraction

az_deg = aob * ERFA.DR2D
el_deg = 90 - zob * ERFA.DR2D

@show az_deg
@show el_deg

# Baseline vector in ECEF frame
bl_ecef = ant_ecef - ref_ecef

# Rotate baseline ECEF XYZ frame by `rob` around the Z axis to get X′Y′Z′ frame
# Rotate baseline X'Y'Z' frame by `dob` around the Y axis to get WUV frame
# Then permute from WUV to UVW
ecef2uvw_halon = [0 1 0
                  0 0 1
                  1 0 0] * ry(-dob) * rz(-hob+ref_lla.lon*ERFA.DD2R)

ecef2uvw_rob   = [0 1 0
                  0 0 1
                  1 0 0] * ry(-dob) * rz(+rob)

uvw_exp = [u_exp, v_exp, w_exp]
uvw_halon = ecef2uvw_halon * bl_ecef
uvw_rob = ecef2uvw_rob * bl_ecef

println("using observed hour angle and declination")
uvw_diff(uvw_halon, uvw_exp)
println()
println("using observed right ascension and declination")
uvw_diff(uvw_rob, uvw_exp)
println()

[uvw_halon, uvw_rob]
end

function uvw_diff(uvw_calc, uvw_exp; show=true)
  uvw_diff = uvw_calc - uvw_exp
  @show uvw_calc
  @show uvw_exp
  @show uvw_diff
  @show hypot(uvw_diff...)
  uvw_diff
end

test_uvw()
