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

# Load the mesh
model = GmshDiscreteModel("data/rectangle_sphere.msh")

# Define the tags of the mesh boundary
dimension = 3
labels = get_face_labeling(model)
tags = get_face_tag(labels, dimension)

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


Γ = BoundaryTriangulation(model, tags="Sphere Surface")
dΓ = Measure(Γ, degree)
nb = get_normal_vector(Γ) # Normal vector to the source boundary


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
# a(u, v) = ∫( ρ * c^2 * (Jinv ∘ xp) * divergence(u) * divergence(v) )*dΩ-
#           ∫( ρ * (Jinv ∘ xp) * (((Hm ∘ xp) ⋅ u) ⋅ ((Hm ∘ xp) ⋅ v)) * ω^2 )*dΩ

a(u, v) = ∫(  ρ * c^2  * (J ∘ xp) * ((Hinv ∘ xp) ⊙ ∇(u)) ⋅ ((Hinv ∘ xp) ⊙ ∇(v)) )*dΩ -
          ∫(  ρ * (J ∘ xp) * (u⋅v) * ω^2 )*dΩ


# Source term
b(v) = ∫((J ∘ xp) * P_0 *(v ⋅ nb))*dΓ
# b(v) = ∫((v⋅nb) * P_0)dΓ # Look that the sign has been changed due to the normal vector points inside the boundary instead of outside it.

# Assembly the system
op = AffineFEOperator(a, b, U, V)
# Solve the system
uh = solve(op)
# uh = (Jinv ∘ xp) * ((Hm ∘ xp) ⋅ Uh)


writevtk(Ω,"./results/results_oldformulation.vtu", cellfields=[ "Re(uh)"=>real(uh), "Im(uh)"=>imag(uh)])