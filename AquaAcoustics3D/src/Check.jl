using Revise
includet("MarineForms.jl")
using .MarineForm
include("Configuration.jl")
using Distributions

# Definition of the objects of the bottom
definition = [
    (RigidSphere, Dict(
    :N => N_rigid_spheres,
    :r_distribution => Normal(r, σ_r),
    :x_range => (-L/2 + (r + 4*σ_r + tol_sphere), L/2 - (r + 4*σ_r + tol_sphere)),
    :y_range => (0+(r + 4*σ_r + + tol_sphere), t_P-(r + 4*σ_r + + tol_sphere)),
    :z_range => (-w/2 + (r + 4*σ_r + tol_sphere), w/2 - (r + 4*σ_r + tol_sphere)))),
    
    (PorousSphere, Dict(
    :N => N_porous_spheres,
    :r_distribution => Normal(r, σ_r),
    :x_range => (-L/2 + (r + 4*σ_r + tol_sphere), L/2 - (r + 4*σ_r + tol_sphere)),
    :y_range => (0+(r + 4*σ_r + + tol_sphere), t_P-(r + 4*σ_r + + tol_sphere)),
    :z_range => (-w/2 + (r + 4*σ_r + tol_sphere), w/2 - (r + 4*σ_r + tol_sphere)))),

    (RigidCockle, Dict(
        :cockle_type => :open, 
        :N => N_open_cockles, 
        :r_distribution   => Normal(by_default_radius, σ_r_cockle),
        :by_default_radius=> by_default_radius,
        :x_range          => (-L/2 + (by_default_radius + 4*σ_r_cockle), L/2 - (by_default_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_radius + 4*σ_r_cockle), t_P - (by_default_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_radius + 4*σ_r_cockle), w/2 - (by_default_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    )),
    
    (RigidCockle, Dict(
        :cockle_type => :closed, 
        :N => N_closed_cockles,
        :r_distribution   => Normal(by_default_radius, σ_r_cockle),
        :by_default_radius=> by_default_radius,
        :x_range          => (-L/2 + (by_default_radius + 4*σ_r_cockle), L/2 - (by_default_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_radius + 4*σ_r_cockle), t_P - (by_default_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_radius + 4*σ_r_cockle), w/2 - (by_default_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    ))
]

# Generate the animals
animals = generate_animals(definition, max_iterations)

center, radius = get_evolving_sphere("./clam_geometries/EntireClam.brep")
center