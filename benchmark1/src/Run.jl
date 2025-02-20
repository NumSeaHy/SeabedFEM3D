"""
Compute the solution of the Helmholtz equation in a fluid and porous domain with PML
using the Finite Element Method in Gridap.jl
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry
using  GridapGmsh

# Include the files with the functions to define the mesh and the physical parameters
include("./Configuration.jl")
include("./AnalyticalSolution.jl")

function RunFEM(name)
    
    # Load the mesh
    model = GmshDiscreteModel("data/" * "$name" * ".msh")

    # Define the finite element space: Raviart-Thomas of order 1
    order = 1
    reffe = ReferenceFE(raviart_thomas, Float64, order)
    V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["Boundaries"], vector_type=Vector{ComplexF64})

    # Define the trial function with null Dirichlet boundary conditions
    uD_null = VectorValue(0.0, 0.0, 0.0)

    U = TrialFESpace(V, [uD_null])

    degree = 2

    Ω = Triangulation(model) # Computational domain
    dΩ = Measure(Ω, degree)

    xp = get_physical_coordinate(Ω)


    Γ = BoundaryTriangulation(model, tags="Sphere")
    dΓ = Measure(Γ, degree)
    nb = get_normal_vector(Γ) # Normal vector to the source boundary

    # Define the analytical solution
    u_analytical(x) = spherical_field(x, 0, H/2, 0, r, k, ρ, P_0, ω)
    AnalyticalBox(x) =  (abs(x[1]) < L/2 && abs(x[3]) < w/2 && x[2] > 0 && (x[2]-H) < 0) ? 1 : 0

    # Define the tensors H, H^-1 and the Jacobian for the fluid PML only in the horizontal direction (quadratic profile)
    k = ω/c

    # Define the tensors H, H^-1 and the Jacobian for the porous PML in the horizontal and vertical direction (quadratic profile)
    γ_1(x) = abs(x[1]) < L/2 ? 1. : 1. + 1im/k * σ_0 * (abs(x[1]) - L/2)^2/d_PML^2
    γ_2(x) = x[2] < 0 ? 1. + 1im/k * σ_0 * x[2]^2/d_PML^2 : ((x[2]-H) > 0 ? 1. + 1im/k * σ_0 * (x[2]-H)^2/d_PML^2 : 1.)
    γ_3(x) = abs(x[3]) < w/2 ? 1. : 1. + 1im/k * σ_0 * (abs(x[3]) - w/2)^2/d_PML^2
    Hm(x) = TensorValue{3,3,ComplexF64}(γ_1(x), 0+0im, 0+0im, 0+0im, γ_2(x), 0+0im, 0+0im, 0+0im, γ_3(x))
    Hinv(x) = inv(Hm(x))
    J(x) = det(Hm(x))
    Jinv(x) = 1/J(x)

    # Bilinear term
    a(u, v) = ∫( ρ * c^2 * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ-
            ∫( ρ * (Jinv ∘ xp) * (((Hm ∘ xp) ⋅ u) ⋅ ((Hm ∘ xp) ⋅ v)) * ω^2 )*dΩ

    # Source term
    b(v) = ∫((v⋅nb) * (-P_0))dΓ 

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    Uh = solve(op)
    uh = (Jinv ∘ xp) * ((Hm ∘ xp) ⋅ Uh)
    u = u_analytical ∘ xp
    λ = AnalyticalBox ∘ xp
    analytical_solution = u_analytical ∘ xp
    error = λ * (uh - analytical_solution)

    # Isolate component by component to comput the error
    uh_x = (u->u[1]) ∘ uh
    uh_y = (u->u[2]) ∘ uh
    uh_z = (u->u[3]) ∘ uh
    u_x = (u->u[1]) ∘ u
    u_y = (u->u[2]) ∘ u
    u_z = (u->u[3]) ∘ u

    writevtk(Ω,"./results/results2" * "$name" * ".vtu", cellfields=[ "Re(uh)"=>real(uh), "Im(uh)"=>imag(uh),
                                                                     "Re(u_analytical)"=>real(analytical_solution), "Im(u_analytical)"=>imag(analytical_solution),
                                                                     "Re(error)"=>real(error),  "Im(error)"=>imag(error)])

    return 100 * sqrt(sum(∫(abs2((λ*u_x-λ*uh_x)))*dΩ)/sum(∫(abs2((λ*u_x)))*dΩ)+
                      sum(∫(abs2((λ*u_y-λ*uh_y)))*dΩ)/sum(∫(abs2((λ*u_y)))*dΩ)+
                      sum(∫(abs2((λ*u_z-λ*uh_z)))*dΩ)/sum(∫(abs2((λ*u_z)))*dΩ))

end 


