"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapGmsh
using GridapPETSc


# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
include("./AnalyticalSolution.jl")

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Define the finite element space: Raviart-Thomas of order 1
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["Boundaries"], vector_type=Vector{ComplexF64})

# Define the trial function with null Dirichlet boundary conditions
uD_null = VectorValue(0.0, 0.0, 0.0)

U = TrialFESpace(V, [uD_null])

degree = 2

Ω = Triangulation(model) # Computational domain
Ω_physical = Triangulation(model, tags=["Fluid"])

dΩ = Measure(Ω, degree)
dΩ_physical = Measure(Ω_physical, degree)

xp = get_physical_coordinate(Ω)

Γ = BoundaryTriangulation(model, tags="Cube")
dΓ = Measure(Γ, degree)
nb = get_normal_vector(Γ) # Normal vector to the source boundary

# Define the analytical solution
u_analytical(x) = spherical_field(x, x₀, y₀, z₀, r, k_F(ω), ρ_F(ω), P_0, ω)
p_analytical(x) = spherical_pressure_field(x, x₀, y₀, z₀, r, k_F(ω), P_0)

# Define the tensors H, H^-1 and the Jacobian for the porous PML in the horizontal and vertical direction (quadratic profile)
γ_1(x) = abs(x[1]) < L/2 ? 1. : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2
γ_2(x) = x[2] < 0 ? 1. + 1im/k_F(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-H) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-H)^2/d_PML^2 : 1.)
γ_3(x) = abs(x[3]) < w/2 ? 1. : 1. + 1im/k_F(ω) * σ_0 * (abs(x[3]) - w/2)^2/d_PML^2
Hm(x) = TensorValue{3,3,ComplexF64}(γ_1(x), 0+0im, 0+0im, 0+0im, γ_2(x), 0+0im, 0+0im, 0+0im, γ_3(x))
Hinv(x) = inv(Hm(x))
J(x) = det(Hm(x))
Jinv(x) = 1/J(x)

# Bilinear term
a(u, v) = ∫( ρ_F(ω) * c_F(ω)^2 * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ-
          ∫( ρ_F(ω) * (Jinv ∘ xp) * (((Hm ∘ xp) ⋅ u) ⋅ ((Hm ∘ xp) ⋅ v)) * ω^2 )*dΩ

# Source term
b(v) = ∫((v⋅nb) * (p_analytical))dΓ # Look that the sign has been changed due to the normal vector points inside the boundary instead of outside it.

# Assembly the system
op = AffineFEOperator(a, b, U, V)

# Define the PETSc options for check the order of magnitude of the residual, only perform 10 itrations and save the solution to then use them as initial guess

options = "-pc_type gamg -ksp_type gmres -ksp_atol 0.0 -ksp_rtol 0.0 -ksp_max_it 10 -ksp_monitor"

println("Performing 10 iterations to get the order of magnitude of the initial residual")

GridapPETSc.with(args=split(options)) do
    ls = PETScLinearSolver()
    fesolver = FESolver(ls)
    xh = zero(U)
    xh = solve!(xh, fesolver, op)

    Uh, cache = xh
    global initial_guess = FEFunction(U, get_free_dof_values(Uh))
    # Initialize a reference container for last residual iteration
    global r_tol = Ref{PetscReal}()
    # Retrieve the total number of KSP iterations
    @check_error_code GridapPETSc.PETSC.KSPGetResidualNorm(cache.ksp[], r_tol)
end

# Maximum number of iteration for the final solution
max_it = 50000

# Define the PETSc options to take the previous order of magnitude of residual as stop criterium as well as the initial guess
options = "-pc_type gamg -ksp_type gmres -ksp_gmres_restart $(Int(0.1*max_it)) -ksp_atol 0.0 -ksp_rtol $(r_tol[]*1e-3) -ksp_max_it $(max_it) -ksp_initial_guess_nonzero true  -ksp_monitor"

println("Solving the system using the previous residual as stop criterium as well as the initial guess")

GridapPETSc.with(args=split(options)) do
        # println(length(Uh.cell_dof_values))
        ls = PETScLinearSolver()
        fesolver = FESolver(ls)
        xh = initial_guess
        xh = solve!(xh, fesolver, op)

        Uh, _ = xh

        # Post-processing
        uh = (Jinv ∘ xp) * ((Hm ∘ xp) ⋅ Uh)
        u = u_analytical ∘ xp
        # error = u - uh

        # Isolate component by component to compute the error
        uh_x = (u->u[1]) ∘ uh
        uh_y = (u->u[2]) ∘ uh
        uh_z = (u->u[3]) ∘ uh
        u_x = (u->u[1]) ∘ u
        u_y = (u->u[2]) ∘ u
        u_z = (u->u[3]) ∘ u

        writevtk(Ω_physical,"./results/results_physical_" * "$name" * ".vtu", cellfields=[ "Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
                                                                                           "Re(u_analytical)"=>real(u), "Im(u_analytical)"=>imag(u),
                                                                                           "Re(error)"=>real(error),  "Im(error)"=>imag(error)])
        
        writevtk(Ω,"./results/results_" * "$name" * ".vtu", cellfields=[ "Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
                                                                         "Re(u_analytical)"=>real(u), "Im(u_analytical)"=>imag(u),
                                                                         "Re(error)"=>real(error),  "Im(error)"=>imag(error)])

        return 100 * sqrt(∑(∫(abs2((u_x-uh_x)))*dΩ_physical)/∑(∫(abs2((u_x)))*dΩ_physical)+
                            ∑(∫(abs2((u_y-uh_y)))*dΩ_physical)/∑(∫(abs2((u_y)))*dΩ_physical)+
                            ∑(∫(abs2((u_z-uh_z)))*dΩ_physical)/∑(∫(abs2((u_z)))*dΩ_physical))
end
