# same as example 150 but with two species and a different system. This time the exact solution is not known, thus we compare with finite-difference approximation of the impedance.
# PDE system on x∈(0,L):
#   C∂ₜu₁ - ∂ₓ(D∂ₓ(u₁u₂)) + (Ru₁u₂ - u₂) = 0
#   C∂ₜu₂ - ∂ₓ(D∂ₓu₂)      + Ru₁u₂        = 0
# consistent with VoronoiFVM convention ∂ₜs(u)+∇⋅j(u)+r(u)=0 and diffusion flux j=-D∇(⋅).
# Boundary conditions used here:
#   Dirichlet: u₂(0,t)=1, u₁(0,t)=excitation(t), u₁(L,t)=0
#   Neumann (natural, zero flux): D∂ₓu₂(L,t)=0

module Example152_Impedance_Multispecies

using VoronoiFVM
using ExtendableGrids: geomspace, simplexgrid, num_nodes
using GridVisualize
using OrdinaryDiffEqSDIRK
using Printf

function main(;
        nref = 0,
        Plotter = nothing,
        verbose = false,
        unknown_storage = :sparse,
        assembly = :edgewise,
        L = 1.0, R = 1.0, D = 1.0, C = 1.0,
        ω0 = 1.0e-3, ω1 = 5.0e1,
        ω_incfactor = 1.1,
        N_periods = :auto,
        Ndt = 200,
        fdtest::Bool = true #false
    )

    # Create array which is refined close to 0
    h0 = 0.005 / 2.0^nref
    h1 = 0.1 / 2.0^nref

    X = geomspace(0, L, h0, h1)

    # Create discretitzation grid
    grid = simplexgrid(X)

    # Create and fill data
    data = (R = R, D = D, C = C)

    # Declare constitutive functions
    flux = function (f, u, edge, data)
        f[1] = data.D * (u[1, 1] * u[2, 1] - u[1, 2] * u[2, 2])
        return f[2] = data.D * (u[2, 1] - u[2, 2])
    end

    storage = function (f, u, node, data)
        f[1] = data.C * u[1]
        return f[2] = data.C * u[2]
    end

    reaction = function (f, u, node, data)
        f[1] = data.R * u[1] * u[2] - u[2]
        return f[2] = data.R * u[2] * u[1]
    end

    excited_bc = 1
    excited_bcval = 1.0
    excited_spec = 1
    meas_bc = 2

    bc = function (f, u, node, data)
        p = parameters(u)
        boundary_dirichlet!(f, u, node; species = 2, region = 1, value = 1.0)
        boundary_dirichlet!(f, u, node; species = 1, region = excited_bc, value = p[1])
        return boundary_dirichlet!(f, u, node; species = 1, region = meas_bc, value = 0.0)
    end

    sys = VoronoiFVM.System(
        grid; unknown_storage = unknown_storage,
        data = data,
        flux = flux,
        storage = storage,
        reaction = reaction,
        bcondition = bc,
        nparams = 1,
        assembly = assembly
    )

    enable_species!(sys, 1, [1])
    enable_species!(sys, 2, [1])
    factory = TestFunctionFactory(sys)
    measurement_testfunction = testfunction(factory, [excited_bc], [meas_bc])

    steadystate = solve(sys; inival = 1.0, params = [1.0])

    function meas_stdy(meas, U)
        if !(typeof(U) <: AbstractMatrix)
            u = reshape(U, sys)
        else
            u = U
        end
        meas[1] = -VoronoiFVM.integrate_stdy(sys, measurement_testfunction, u, params = [1.0])[excited_spec]
        return nothing
    end

    function meas_tran(meas, U)
        if !(typeof(U) <: AbstractMatrix)
            u = reshape(U, sys)
        else
            u = U
        end
        meas[1] = -VoronoiFVM.integrate_tran(sys, measurement_testfunction, u, params = [1.0])[excited_spec]
        return nothing
    end

    dmeas_stdy = measurement_derivative(sys, meas_stdy, steadystate)
    dmeas_tran = measurement_derivative(sys, meas_tran, steadystate)
    meas_tran_ref = zeros(1)
    meas_stdy_ref = zeros(1)
    meas_cos = zeros(1)
    meas_sin = zeros(1)
    meas_tran(meas_tran_ref, steadystate)
    meas_stdy(meas_stdy_ref, steadystate)

    # Create Impeadancs system from steady state
    isys = VoronoiFVM.ImpedanceSystem(sys, steadystate)

    # Prepare recording of impedance results
    allomega = zeros(0)

    # for calculated data
    allI0 = zeros(Complex{Float64}, 0)
    allIL = zeros(Complex{Float64}, 0)

    # for exact data
    allIx0 = zeros(Complex{Float64}, 0)
    allIxL = zeros(Complex{Float64}, 0)

    ω = ω0

    UZ = unknowns(isys)

    outflux_ref = zeros(2)
    outflux_cos = zeros(2)
    outflux_sin = zeros(2)
    nnodes = num_nodes(grid)
    lastedge = (nnodes - 1):nnodes
    @views flux(outflux_ref, steadystate[:, lastedge], nothing, data)

    if verbose
        println("Stationary computation done, starting impedance computation.")
    end
    while ω < ω1
        # solve impedance system
        solve!(UZ, isys, ω)

        # calculate approximate solution
        # obtain measurement in frequency  domain
        IL = impedance(isys, ω, steadystate, dmeas_stdy, dmeas_tran)
        if verbose && ω ≈ ω0
             @printf("Successfully calculated impedance (prints only for the first frequency)")
        end
        # record approximate solution
        push!(allomega, ω)
        push!(allIL, IL)

        if fdtest
            # compute reference using finite difference approximation
            amplitude = 1.0e-6
            data_perturbed = (R = R, D = D, C = C, ω = ω)
            #change boundary condition to reflect the perturbation
            bc_cos = function (f, u, node, data)
                p = parameters(u)
                boundary_dirichlet!(f, u, node; species = 2, region = 1, value = 1.0)
                boundary_dirichlet!(f, u, node; species = 1, region = excited_bc, value = 1 + amplitude * cos(data.ω * node.time))
                return boundary_dirichlet!(f, u, node; species = 1, region = meas_bc, value = 0.0)
            end

            sys_cos = VoronoiFVM.System(
                grid; unknown_storage = unknown_storage,
                data = data_perturbed,
                flux = flux,
                storage = storage,
                reaction = reaction,
                bcondition = bc_cos,
                nparams = 0, #we no longer track derivative with respect to parameters
                assembly = assembly
            )
            enable_species!(sys_cos, 1, [1])
            enable_species!(sys_cos, 2, [1])
            inival = unknowns(sys_cos)
            inival[1, :] .= steadystate[1, :]
            inival[2, :] .= steadystate[2, :]


            #same for the sine perturbation
            bc_sin = function (f, u, node, data)
                p = parameters(u)
                boundary_dirichlet!(f, u, node; species = 2, region = 1, value = 1.0)
                boundary_dirichlet!(f, u, node; species = 1, region = excited_bc, value = 1 + amplitude * sin(data.ω * node.time))
                return boundary_dirichlet!(f, u, node; species = 1, region = meas_bc, value = 0.0)
            end

            sys_sin = VoronoiFVM.System(
                grid; unknown_storage = unknown_storage,
                data = data_perturbed,
                flux = flux,
                storage = storage,
                reaction = reaction,
                bcondition = bc_sin,
                nparams = 0,
                assembly = assembly
            )
            enable_species!(sys_sin, 1, [1])
            enable_species!(sys_sin, 2, [1])

            dt = (2 * π / ω) / Ndt
            period = 2 * π / ω
            control = VoronoiFVM.SolverControl(Δt_max = dt, Δt_min = dt, Δt = dt, Δu_opt = 1.0e10)

            # Either run a fixed number of periods or use a fixed-point iteration
            # on the period map (automatic mode).
            tsol_cos = nothing
            tsol_sin = nothing
            if N_periods == :auto
                max_periods = 200
                # Use several periods before checking convergence.
                min_periods = 3
                cycle_rtol = 1.0e-7
                cycle_atol = 1.0e-9

                u0_cos = copy(steadystate)
                u0_sin = copy(steadystate)
                converged = false

                for iperiod in 1:max_periods
                    tsol_cos = solve(
                        sys_cos; inival = u0_cos, times = (0.0, period), force_first_step = true, control = control
                    )
                    tsol_sin = solve(
                        sys_sin; inival = u0_sin, times = (0.0, period), force_first_step = true, control = control
                    )

                    u1_cos = tsol_cos[end]
                    u1_sin = tsol_sin[end]

                    if iperiod >= min_periods
                        conv_cos = isapprox(dofs(u1_cos), dofs(u0_cos), rtol = cycle_rtol, atol = cycle_atol)
                        conv_sin = isapprox(dofs(u1_sin), dofs(u0_sin), rtol = cycle_rtol, atol = cycle_atol)
                        if conv_cos && conv_sin
                            converged = true
                            break
                        end
                    end

                    u0_cos .= u1_cos
                    u0_sin .= u1_sin
                end
                @assert converged "Automatic period selection did not converge for ω=$ω within $max_periods periods."
            else
                tend = (N_periods + 1) * period
                # Compute a sufficiently long transient and evaluate the impedance on the last period.
                tsol_cos = solve(
                    sys_cos; steadystate, times = (0.0, tend), force_first_step = true, control = control
                )
                tsol_sin = solve(
                    sys_sin; steadystate, times = (0.0, tend), force_first_step = true, control = control
                )
            end


            #and use the results to compute the impedance using finite difference approximation

            time_impedance = zeros(ComplexF64, Ndt)
            if ω ≈ ω0 && verbose
                @show length(tsol_cos.t) Ndt + 1 period / dt
                @show (tsol_cos.t[2:end] - tsol_cos.t[1:(end - 1)])[(end - 10):end]
                @show dt
            end

            j_last_period = length(tsol_cos.t) - Ndt
            for i in 1:Ndt
                j = j_last_period + i
                time = tsol_cos.t[j]

                @assert isapprox(time, tsol_sin.t[j], rtol = 1.0e-5)

                u_cos = tsol_cos.u[j]
                u_sin = tsol_sin.u[j]

                #compute flux at the boundary for both solutions and subtract the reference flux to get the flux perturbation

                endcos = view(u_cos, :, lastedge)
                endsin = view(u_sin, :, lastedge)

                flux(outflux_cos, endcos, nothing, data)
                flux(outflux_sin, endsin, nothing, data)

                outflux_cos .-= outflux_ref
                outflux_sin .-= outflux_ref
                tau = 1 / (X[end] - X[end - 1])


                time_impedance[i] = (outflux_cos[1] * tau + 1im * outflux_sin[1] * tau) / (amplitude * exp(1im * ω * time))

            end
            IxL = length(time_impedance) / sum(time_impedance)
            if verbose
                 ratio = IL / IxL

                @printf(
                    "Finite difference approximation of impedance at ω = %10.5g: %10.5g%+10.5gi, calculated impedance: %10.5g%+10.5gi, ratio distance to one: %10.5g\n",
                    ω,
                    real(IxL), imag(IxL),
                    real(IL), imag(IL),
                    abs(ratio - 1.0)
                )
            end
            push!(allIxL, IxL)
        end

        # increase omega
        ω = ω * ω_incfactor
    end

    vis = GridVisualizer(; Plotter = Plotter, legend = :rt)
    if fdtest
        scalarplot!(
            vis, real(allIxL), imag(allIxL); label = "finite difference", color = :red, linestyle = :dot
        )
    end
    scalarplot!(
        vis, real(allIL), imag(allIL); label = "calc", show = true, clear = false, color = :blue, linestyle = :solid
    )
    if fdtest
        println("Ratio of calculated impedance to finite difference impedance, should be close to one: ")
        avg_ratio = sum(allIL ./ allIxL) / length(allomega)
        @printf("Minimum distance to one: %.3e, Average value: %.3g%+.3gi, Maximum distance to one: %.3e\n", minimum(abs.(allIL ./ allIxL .- 1)), real(avg_ratio), imag(avg_ratio), maximum(abs.(allIL ./ allIxL .- 1)))
        #@show minimum(abs.(allIL ./ allIxL)), sum(allIL ./ allIxL) / length(allomega), maximum(abs.(allIL ./ allIxL))
    end
    return sum(allIL)
end


using Test

function runtests()
    testval = 50.960361928838 + 4.510584656768053im
    for unknown_storage in (:sparse, :dense)
        for assembly in (:edgewise, :cellwise)
            @test main(; unknown_storage, assembly) ≈ testval
        end
    end
    return
end

end #end of module
