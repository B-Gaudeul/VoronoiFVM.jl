module test_timestepnorm

using Test
using VoronoiFVM

function test()
    c = VoronoiFVM.SolverControl()
    @test VoronoiFVM.timestep_norm(nothing, [3.0, 4.0], [0.0, 0.0], :L1) ≈ 7.0
    @test VoronoiFVM.timestep_norm(nothing, [3.0, 4.0], [0.0, 0.0], :L2) ≈ 5.0
    @test VoronoiFVM.timestep_norm(nothing, [3.0, 4.0], [0.0, 0.0], :Linfty) ≈ 4.0

    c.delta_norm = :L2
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈ 5.0
    c.delta_norm = 3.0
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈ (3.0^3 + 4.0^3)^(1 / 3)
    c.delta_norm = 1.0
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈ 7.0
    c.delta_norm = Inf
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈ 4.0
    c.delta_norm = (system, u, v, t, Δt) -> VoronoiFVM.timestep_norm(system, u, v, :L1)
    c.delta = nothing
    c2 = VoronoiFVM.SolverControl(delta_norm = :L1, delta = nothing)
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈
          VoronoiFVM.timestep_delta(c2, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1)
    c.delta_norm = (system, u, v, t, Δt) -> VoronoiFVM.norm(system, u - v, 3.0)
    c2.delta_norm = 3.0
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) ≈
          VoronoiFVM.timestep_delta(c2, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1)
    c.delta_norm = (system, u, v, t, Δt) -> 11.0
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) == 11.0
    c.delta = (system, u, v, t, Δt) -> 42.0
    @test VoronoiFVM.timestep_delta(c, nothing, [3.0, 4.0], [0.0, 0.0], 0.0, 0.1) == 42.0

    @test_throws ArgumentError VoronoiFVM.timestep_norm(nothing, [1.0], [0.0], :Wasserstein)
    c.delta = nothing
    c.delta_norm = 0.0
    @test_throws ArgumentError VoronoiFVM.timestep_delta(c, nothing, [1.0], [0.0], 0.0, 0.1)
    return true
end

end
