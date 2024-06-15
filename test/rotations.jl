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

# There has to be an easier way to do this, right???  We convert to Float64s
# because using `@test all(r123_val .≈ r123_erfa_val)` sometimes throws a weird
# exception (for reasons unknown) even though other times it works fine.  If we
# could understand/fix that latter problem then we wouldn't need to convert to
# Float64s.
r123_float = Symbolics.value.(parent(r123_val))
r123_erfa_float = Symbolics.value.(parent(r123_erfa_val))

@testset "rotations" begin
    @test all(r123_num .≈ r123_val)
    @test all(r123_float .≈ r123_erfa_float)
    @test all((r123_num * v) .≈ (r123_val * v))
    @test all((r123_float * v) .≈ (r123_erfa_float * v))
end
