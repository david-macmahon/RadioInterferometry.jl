let 
    
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

uvw_expected = [-617.7848812531224, -1504.9793145403532, 1412.645848308487]

xyz_uvw = radec2uvw(ra, dec, jd, obslla)

antuvw1 = xyz_uvw * antxyz
antuvw2 = radec2uvw(ra, dec, jd, obslla, antxyz)

uvw1 = antuvw1[:,2] - antuvw1[:,1]
uvw2 = radec2uvw(ra, dec, jd, obslla, antxyz, [(1,2)])

@testset "radec2uvw" begin
    @test antuvw1 == antuvw2
    @test uvw1 == uvw_expected
    @test uvw2 == uvw_expected
end

end # let