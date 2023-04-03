module GCRSUVW

export radec2uvw

import Rotations: RotY, RotZ
import ERFA: DAS2R, c2t06a, apco13, ldsun, ab, s2c, c2s, utctai, taitt

function radec2uvw(ra, dec, jdutc, obslla;
                   jdutc2=0.0, dut1=0.0, xp=0.0, yp=0.0)
    # Get TT
    tai1, tai2 = utctai(jdutc, jdutc2)
    tt1, tt2 = taitt(tai1, tai2)
    
    # Convert lat/lon to radians
    lat = deg2rad(obslla.lat)
    lon = deg2rad(obslla.lon)
    alt = obslla.alt

    # Convert ra,dec from ICRS to GCRS.  For completeness, this involves space
    # motion, parallax, light deflection, and aberration, but here we just do
    # light deflection due to the sun and aberration.  We use `apco13` rather
    # than `apcg13` because the former includes observatory position/velocity in
    # the astrometry parameters used to calculate light deflection and
    # aberration whereas the latter does not.
    astrom = apco13(jdutc, jdutc2, dut1, lon, lat, alt, xp, yp, 0, 0, 0, 0)[1]
    srcvec = s2c(ra, dec)
    ldvec = ldsun(srcvec, astrom.eh, astrom.em)
    gcvec = ab(ldvec, astrom.v, astrom.em, astrom.bm1)
    ra_gcrs, dec_gcrs = c2s(gcvec)

    # Rotation matrix to transform GCRS frame to UVW frame in direction of
    # (ra_gcrs, dec_gcrs).
    gcrs2uvw = [0 1 0
                0 0 1
                1 0 0] * RotY(dec_gcrs) * RotZ(-ra_gcrs)

    # Rotation matrix to transform ITRF (or XYZ) to GCRS (oriented) is the
    # transpose of ERFA's celestial-to-terrestrial ("c2t") matrix.
    xyz2gcrs = c2t06a(tt1, tt2, jdutc, jdutc2+dut1/86400, xp*DAS2R, yp*DAS2R)'

    # xyz to gcrs to uvw
    gcrs2uvw * xyz2gcrs
end

function radec2uvw(ra, dec, jdutc, obslla, xyz;
                   jdutc2=0.0, dut1=0.0, xp=0.0, yp=0.0)
    radec2uvw(ra, dec, jdutc, obslla; jdutc2, dut1, xp, yp) * xyz
end

function radec2uvw(ra, dec, jdutc, obslla, antxyz, bls;
                   jdutc2=0.0, dut1=0.0, xp=0.0, yp=0.0)
    antuvw = radec2uvw(ra, dec, jdutc, obslla, antxyz; jdutc2, dut1, xp, yp)
    # Compute baselines
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

end # module GCRSUVW