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
        :name => "Cockle",
        :type => :open, 
        :N => N_open_cockles, 
        :r_distribution   => Normal(by_default_cockle_radius, σ_r_cockle),
        :by_default_radius=> by_default_cockle_radius,
        :x_range          => (-L/2 + (by_default_cockle_radius + 4*σ_r_cockle), L/2 - (by_default_cockle_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_cockle_radius + 4*σ_r_cockle), t_P - (by_default_cockle_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_cockle_radius + 4*σ_r_cockle), w/2 - (by_default_cockle_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    )),
    
    (RigidCockle, Dict(
        :name => "Cockle",
        :type => :closed, 
        :N => N_closed_cockles,
        :r_distribution   => Normal(by_default_cockle_radius, σ_r_cockle),
        :by_default_radius=> by_default_cockle_radius,
        :x_range          => (-L/2 + (by_default_cockle_radius + 4*σ_r_cockle), L/2 - (by_default_cockle_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_cockle_radius + 4*σ_r_cockle), t_P - (by_default_cockle_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_cockle_radius + 4*σ_r_cockle), w/2 - (by_default_cockle_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    )),


    (RigidQueenScallop, Dict(
        :name => "QueenScallop",
        :type => :closed, 
        :N => N_closed_queenscallops,
        :r_distribution   => Normal(by_default_queenscallop_radius, σ_r_cockle),
        :by_default_radius=> by_default_queenscallop_radius,
        :x_range          => (-L/2 + (by_default_queenscallop_radius + 4*σ_r_cockle), L/2 - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_queenscallop_radius + 4*σ_r_cockle), t_P - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_queenscallop_radius + 4*σ_r_cockle), w/2 - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    )),

    (RigidQueenScallop, Dict(
        :name => "QueenScallop",
        :type => :closed, 
        :N => N_closed_queenscallops,
        :r_distribution   => Normal(by_default_queenscallop_radius, σ_r_cockle),
        :by_default_radius=> by_default_queenscallop_radius,
        :x_range          => (-L/2 + (by_default_queenscallop_radius + 4*σ_r_cockle), L/2 - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :y_range          => (0 + (by_default_queenscallop_radius + 4*σ_r_cockle), t_P - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :z_range          => (-w/2 + (by_default_queenscallop_radius + 4*σ_r_cockle), w/2 - (by_default_queenscallop_radius + 4*σ_r_cockle)),
        :α_range          => (0, 2π),
        :β_range          => (0, π/2),
        :γ_range          => (0, 2π)
    ))


]

# Generate the animals
animals = generate_animals(definition, max_iterations)
