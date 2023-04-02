using Symbolics
using RadioInterferometry

@testset "extension" begin
    @test length(methods(RadioInterferometry.signclean)) == 1
end