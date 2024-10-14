using SpecialFunctions
using Gridap
using GridapGmsh

# Include the files with the functions to define the mesh and the physical parameters
# include("./Configuration.jl")

function distance(x, x₀, y₀, z₀)
    return sqrt((x[1]-x₀)^2 + (x[2]-y₀)^2 + (x[3]-z₀)^2)
end

function spherical_field(x, x₀, y₀, z₀, r, kF, ρF, P0, ω)
    C = P0 * kF *  1/(2 * ρF * ω^2 * hankelh1(1/2, kF*r) * sqrt(pi/(2*kF*r)))
    
    ux = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) - hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) - hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * (x[1]-x₀)/distance(x, x₀, y₀, z₀)
    
    uy = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) - hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) - hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * (x[2]-y₀)/distance(x, x₀, y₀, z₀)
    
    uz = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) - hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) - hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * (x[3]-z₀)/distance(x, x₀, y₀, z₀)

    return VectorValue(ux, uy, uz)
    
end

# # Load the mesh
# model = GmshDiscreteModel("data/rectangle_sphere.msh")

# Ω = Triangulation(model) # Computational domain

# xp = get_physical_coordinate(Ω)

# k = ω/c

# # Define the fields
# u_analytical(x) = spherical_field(x, x₀, y₀, z₀, r, k, ρ, P_0, ω)

# u = u_analytical ∘ xp

# writevtk(Ω,"./results/resultsanalytical.vtu", cellfields=[ "Re(uh)"=>real(u), "Im(uh)"=>imag(u)])