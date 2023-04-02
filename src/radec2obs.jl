# Helper function to simplify ERFA usage

module RADec2Obs

export radec2obs

import ERFA: DAS2R, atco13

"""
    radec2obs(α, δ, jdutc, obslla;
              jdutc2=0.0,
              dut1=0.0,
              xp=0.0,
              yp=0.0,
              tc::Union{Nothing,Real}=nothing,
              rh=0.5,
              wl=0.21e6
             ) -> (; aob, zob, hob, dob, rob, eo)

Convert ICRS right ascension `α` and declination `δ` to observed hour angle
`hob` and observed declination `δob` as seen from the location specified by
`obslla` at the Julian date `jdutc+jdutc2`.  Earth orientation parameters
`dut1`, `xp`, and `yp` may be specified if desired, otherwise the default to
0.0.

If `tc` (temperature in Celcius) is given (and not `nothing`), then refraction
will be included.  Pressure will be estimated from `altitude` and `tc`.
Relative humidity defaults to 50%.  Wavelength defaults to 21cm.

`obslla` must have `lat`, `lon` and `alt` fields.  It can be a NamedTuple or a
`Geodesy.LLA` object.  The `lat` and `lon` fields should be given in degrees.
They may be numeric values or Strings in decimal ("dd.ddd") or sexagesimal
("dd:mm:ss.sss") format.  The `alt` field is height above the WGS84 ellipsoid in
meters.

This is a convenience wrapper around `ERFA.atco13`, but note that some of the
input units differ.

!!! note
    Be sure to use the right units for the various parameters!

| Parameter | Description            | Units   |
|-----------|:-----------------------|:--------|
|         α | Right ascention (ICRS) | radians |
|         δ | Declination (ICRS)     | radians |
|     jdutc | Julian Date (UTC)      | days    |
| longitude | Longitude              | degrees |
|  latitude | Geodetic latitude      | degrees |
|  altitude | Altitude               | meters  |
|    jdutc2 | Julian Date (UTC)      | days    |
|      dut1 | UTC-UT1                | seconds |
|        xp | Polar motion X         | arcsec  |
|        yp | Polar motion Y         | arcsec  |
|        tc | Ambient temperature    | Celcius |
|        rh | Realtive humidity      | %/100.0 |
|        wl | Wavelength             | µmeters |

The returned NamedTuple contains these fields (as returned by `ERFA.atco13`):

| Name | Description                          | Units   |
|------|:-------------------------------------|:--------|
| aob  | Observed azimuth (N=0,E=90)          | radians |
| zob  | Observed zenith distance             | radians |
| hob  | Observed hour angle                  | radians |
| dob  | Observed declination                 | radians |
| rob  | Observed right ascension (CIO-based) | radians |
| eo   | Equation of the origins (ERA-GST)    | radians |
"""
function radec2obs(α, δ, jdutc, obslla;
                   jdutc2=0.0,
                   dut1=0.0,
                   xp=0.0,
                   yp=0.0,
                   tc::Union{Nothing,Real}=nothing,
                   _ignored_kwargs...
                  )
    if tc === nothing
        phpa = tc = rh = wl = 0.0
    else
        # phpa formula from ERFA.atco13 documentation
        phpa = 1013.25 * exp(-altitude/(29.3*(273.15+tc)))
        rh = 0.5
        wl = 0.21e6
    end

    atco13(α, δ,
        0, 0, 0, 0,
        jdutc, jdutc2, dut1,
        deg2rad(obslla.lon), deg2rad(obslla.lat), obslla.alt,
        xp*DAS2R, yp*DAS2R,
        phpa, tc, rh, wl
    ) |> NamedTuple{(:aob, :zob, :hob, :dob, :rob, :eo)}
end

end # module RADec2Obs