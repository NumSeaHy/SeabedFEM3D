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
using Revise

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
includet("./AnalyticalIncident.jl")
using .AnalyticalIncident

# Load the mesh
model = GmshDiscreteModel("data/mesh.msh")

# Define the incident wave
u_incident(x) = analytical_incident(x, ω, P0)

# Define the finite element space: Raviart-Thomas of order 1
order = 1
reffe = ReferenceFE(raviart_thomas, Float64, order)
V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["Boundaries", "Animals"], vector_type=Vector{ComplexF64})

# Define the trial function with null Dirichlet boundary conditions
uD_null = VectorValue(0.0, 0.0, 0.0)

U = TrialFESpace(V, [uD_null, u_incident])

degree = 2

Ω = Triangulation(model) # Computational domain
Ω_physical = Triangulation(model, tags=["Fluid", "Porous"])

dΩ = Measure(Ω, degree)
dΩ_physical = Measure(Ω_physical, degree)

xp = get_physical_coordinate(Ω)

# Define the properties of each of the mediums
ρ(x) = x[2] - t_P < 0 ? ρ_P(ω) : ρ_F(ω)
K(x) = x[2] - t_P < 0 ? K_P(ω) : K_F(ω)

# Define the tensors H, H^-1 and the Jacobian for the porous PML in the horizontal and vertical direction (quadratic profile)
γ_1(x) = abs(x[1]) < L/2  ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2)
γ_2(x) = x[2] < 0 ? 1. + 1im/k_P(ω) * σ_0 * x[2]^2/d_PML^2 : ((x[2]-(t_P+t_F)) > 0 ? 1. + 1im/k_F(ω) * σ_0 * (x[2]-(t_P+t_F))^2/d_PML^2 : 1.)
γ_3(x) = abs(x[3]) < w/2 ? 1. : (x[2]-t_P <= 0  ? 1. + 1im/k_P(ω) * σ_0 * (abs(x[3]) - w/2)^2/d_PML^2 : 1. + 1im/k_F(ω) * σ_0 * (abs(x[3]) - w/2)^2/d_PML^2)
Hm(x) = TensorValue{3,3,ComplexF64}(γ_1(x), 0+0im, 0+0im, 0+0im, γ_2(x), 0+0im, 0+0im, 0+0im, γ_3(x))
Hinv(x) = inv(Hm(x))
J(x) = det(Hm(x))
Jinv(x) = 1/J(x)

# Bilinear term
a(u, v) = ∫( (K ∘ xp) * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ-
          ∫( (ρ ∘ xp) * (Jinv ∘ xp) * (((Hm ∘ xp) ⋅ u) ⋅ ((Hm ∘ xp) ⋅ v)) * ω^2 )*dΩ

# Source term
b(v) = 0.0 + 0.0im

# Assembly the system
op = AffineFEOperator(a, b, U, V)

println("The number of degrees of freedom is $(get_matrix(op).m)")

# Define the PETSc options for check the order of magnitude of the residual, only perform 10 itrations and save the solution to then use them as initial guess

# options = "-ksp_converged_reason -pc_factor_mat_solver_type mumps"
# options = "-ksp_converged_reason -ksp_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
options = "-ksp_converged_reason -ksp_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"

println("Using MUMPS solver to solve the linear system of equation")

GridapPETSc.with(args=split(options)) do
    ls = PETScLinearSolver()
    Uh = solve(ls, op)
    # Post-processing
    uh = (Jinv ∘ xp) * ((Hm ∘ xp) ⋅ Uh)
    writevtk(Ω,"./results/result.vtu", cellfields=[ "Re(u)"=>real(uh), "Im(u)"=>imag(uh)])
      
end

