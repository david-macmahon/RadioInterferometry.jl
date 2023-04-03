using Test
using Documenter
using RadioInterferometry

DocMeta.setdocmeta!(RadioInterferometry,
                    :DocTestSetup,
                    :(using RadioInterferometry);
                    recursive=true)

randlat() = (rand() - 0.5) * π
randlon() = (rand() - 0.5) * 2π
randha() = (rand() - 0.5) * 6
randdec = randlat

# Run doc tests
doctest(RadioInterferometry, manual=false)

@testset "string macros" begin
    @testset "dms_str tests" begin
        @test dms"1:30:00" == 1.5
        @test dms"-00:30:00" == -0.5
        @test dms"90:00:00" == 90.0
        @test dms"90:00:00"ha == 6.0
        @test dms"90:00:00"rad == π/2
    end

    @testset "hms_str tests" begin
        @test hms"1:30:00" == 1.5
        @test hms"-00:30:00" == -0.5
        @test hms"06:00:00" == 6.0
        @test hms"06:00:00"deg == 90.0
        @test hms"06:00:00"rad == π/2
    end
end

include("rotations.jl")

@testset "xyz2enu tests" begin
    # Atlantic ocean
    @test xyz2enu(0, 0) == [0 1 0
                            0 0 1
                            1 0 0]
    # North Pole
    @test xyz2enu(π/2, 0) ≈ [0 1 0
                            -1 0 0
                             0 0 1]
    # Indian Ocean
    @test xyz2enu(0, π/2) ≈ [-1 0 0
                              0 0 1
                              0 1 0]
end

include("extension.jl")
include("topouvw.jl")
include("gcrsuvw.jl")

nothing
