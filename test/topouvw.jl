let # Use let block to avoid polluting variable namespace

ra  = 3.2686162621088966
dec = 0.03582094364245917
jd  = 2.4598932634780095e6

obslla = (
    lat =   34.07861111111111,
    lon = -107.61777777777777,
    alt = 2124.0
)

antxyz = [
     75.3561   -1407.4597
    489.4072     -77.4865
    721.5445    -735.1909
]

uvw_expected = [-617.4059528467405, -1505.1348066165679, 1412.6458483084875]

o = radec2obs(ra, dec, jd, obslla)

xyz_uvw1 = hadec2uvw(o.hob, o.dob, deg2rad(obslla.lon))
xyz_uvw2 = hadec2uvw(o.hob, o.dob, obslla)

antuvw1 = xyz_uvw1 * antxyz
antuvw2 = hadec2uvw(o.hob, o.dob, deg2rad(obslla.lon), antxyz)
antuvw3 = hadec2uvw(o.hob, o.dob, obslla, antxyz)

uvw1 = antuvw1[:,2] - antuvw1[:,1]
uvw2 = hadec2uvw(o.hob, o.dob, deg2rad(obslla.lon), antxyz, [(1,2)])
uvw3 = hadec2uvw(o.hob, o.dob, obslla, antxyz, [(1,2)])

@testset "hadec2uvw" begin
    @test xyz_uvw1 == xyz_uvw2
    @test antuvw1 == antuvw2
    @test antuvw1 == antuvw3
    @test uvw1 == uvw_expected
    @test uvw2 == uvw_expected
    @test uvw3 == uvw_expected
end

end # let