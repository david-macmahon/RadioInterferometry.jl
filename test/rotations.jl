using RadioInterferometry
using Symbolics

import RadioInterferometry: rx, ry, rz

@variables ϕ θ ψ

r123 = rz(ψ)ry(θ)rx(ϕ)

r123_erfa =  [
     cos(ψ)cos(θ)   cos(ψ)sin(θ)sin(ϕ)+sin(ψ)cos(ϕ)  -cos(ψ)sin(θ)cos(ϕ)+sin(ψ)sin(ϕ)
    -sin(ψ)cos(θ)  -sin(ψ)sin(θ)sin(ϕ)+cos(ψ)cos(ϕ)   sin(ψ)sin(θ)cos(ϕ)+cos(ψ)sin(ϕ)
           sin(θ)                     -cos(θ)sin(ϕ)                      cos(θ)cos(ϕ)
]

ϕval = 2π*rand()
θval = 2π*rand()
ψval = 2π*rand()

v = 2*rand(3, 1000) .- 1

r123_num      = rz(ψval)ry(θval)rx(ϕval)
r123_val      = substitute.(r123     , [[ϕ=>ϕval, θ=>θval, ψ=>ψval]])
r123_erfa_val = substitute.(r123_erfa, [[ϕ=>ϕval, θ=>θval, ψ=>ψval]])

@testset "rotations" begin
    @test all(r123_num .≈ r123_val)
    # This sometimes passes and sometimes fails with:
    # `MethodError: no method matching decompose(::Num)`
    @test all(r123_val .≈ r123_erfa_val) skip=true
    @test all((r123_num * v) .≈ (r123_val * v))
    # This sometimes passes and sometimes fails with:
    # `MethodError: no method matching decompose(::Num)`
    @test all((r123_val * v) .≈ (r123_erfa_val * v)) skip=true
end
