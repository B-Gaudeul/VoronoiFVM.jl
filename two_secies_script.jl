using VoronoiFVM
using ExtendableGrids: geomspace, simplexgrid
using GridVisualize
using OrdinaryDiffEqSDIRK

using Plots

nref = 0
Plotter = Plots
verbose = false
unknown_storage = :sparse
assembly = :edgewise
time_embedding = :none

L = 1.0
R = 1.0
D = 1.0
C = 1.0
ω0 = 1.0e-3
ω1 = 5.0e1

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
    f[1] = data.D * (u[1, 1]*u[2, 1] - u[1, 2]*u[2, 2])
    f[2] = data.D * (u[2, 1] - u[2, 2])
end

storage = function (f, u, node, data)
    f[1] = data.C * u[1]
    f[2] = data.C * u[2]
end

reaction = function (f, u, node, data)
    f[1] = data.R * u[1]*u[2] - u[2]
    f[2] = data.R * u[2]* u[1]
end

excited_bc= 1
excited_bcval = 1.0
excited_spec = 1
meas_bc = 2

bc = function (f, u, node, data)
    p = parameters(u)
    boundary_dirichlet!(f, u, node; species = 2, region = 1, value = 1.0)
    boundary_dirichlet!(f, u, node; species = 1, region = excited_bc, value = p[1])
    boundary_dirichlet!(f, u, node; species = 1, region = meas_bc, value = 0.0)
end

sys = VoronoiFVM.System(grid; unknown_storage = unknown_storage,
                            data = data,
                            flux = flux,
                            storage = storage,
                            reaction = reaction,
                            bcondition = bc,
                            nparams = 1,
                            assembly = assembly)

enable_species!(sys, 1, [1])
enable_species!(sys, 2, [1])
factory = TestFunctionFactory(sys)
measurement_testfunction = testfunction(factory, [excited_bc], [meas_bc])


tend=1.0
# Transient solvers 
#tsol=solve(sys; inival = 0.0, params = [1.0], times=(0.0,tend),force_first_step=true)

inival=unknowns(sys,inival=0)
problem = ODEProblem(sys,inival,(0,tend); params = [1.0])
odesol = solve(problem,ImplicitEuler())
tsol=reshape(odesol,sys)
steadystate=tsol.u[end]