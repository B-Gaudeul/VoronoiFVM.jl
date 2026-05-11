#=

# 450: Timestep delta callbacks
 ([source code](@__SOURCE_URL__))

Demonstrate that callback-based `delta_norm` can reproduce symbolic and numeric norm choices.
=#

module Example450_TimestepDeltaCallbacks

using VoronoiFVM
using Test

function main()
    u = [3.0, 4.0]
    v = [0.0, 0.0]

    c_sym = VoronoiFVM.SolverControl(delta_norm = :L1, delta = nothing)
    c_cb_sym = VoronoiFVM.SolverControl(
        delta_norm = (system, u, v, t, Δt) -> VoronoiFVM.timestep_norm(system, u, v, :L1),
        delta = nothing,
    )

    c_num = VoronoiFVM.SolverControl(delta_norm = 3.0, delta = nothing)
    c_cb_num = VoronoiFVM.SolverControl(
        delta_norm = (system, u, v, t, Δt) -> VoronoiFVM.norm(system, u - v, 3.0),
        delta = nothing,
    )

    d_sym = VoronoiFVM.timestep_delta(c_sym, nothing, u, v, 0.0, 0.1)
    d_cb_sym = VoronoiFVM.timestep_delta(c_cb_sym, nothing, u, v, 0.0, 0.1)

    d_num = VoronoiFVM.timestep_delta(c_num, nothing, u, v, 0.0, 0.1)
    d_cb_num = VoronoiFVM.timestep_delta(c_cb_num, nothing, u, v, 0.0, 0.1)

    return (d_sym = d_sym, d_cb_sym = d_cb_sym, d_num = d_num, d_cb_num = d_cb_num)
end

function runtests()
    r = main()
    @test r.d_sym ≈ r.d_cb_sym
    @test r.d_num ≈ r.d_cb_num
    return nothing
end

end
