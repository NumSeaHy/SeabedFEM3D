module MarineForm

# Loading the packages
using Distributions
using Gmsh

# Export the functions of this module that are going to be used in other modules
export generate_animals, PorousSphere, RigidSphere, generate_geometry

"Definition of the abstract type MarineForms."
abstract type MarineForms end

"Definition of the abstract type Sphere."
abstract type Sphere <: MarineForms end

"Definition of the abstract type Cockle"
abstract type Cockle <: MarineForms end

"Definition of the Porous Sphere object."
struct PorousSphere <: Sphere
    x::Float64 # x-coordinate of the center of the Sphere
    y::Float64 # y-coordinate of the center of the Sphere
    z::Float64 # Z-coordinate of the center of the sphere
    r::Float64 # Radius of the Sphere
    tol::Float64 # Tolerance of the Sphere to avoid collisions with the boundaries of the physical domain
end

"Definition of the Rigid Sphere object."
struct RigidSphere <: Sphere
    x::Float64 # x-coordinate of the center of the Sphere
    y::Float64 # y-coordinate of the center of the Sphere
    z::Float64 # z-coordinate of the center of the Sphere
    r::Float64 # Radius of the Sphere
    tol::Float64 # Tolerance of the Sphere to avoid collisions with the boundaries of the physical domain
end

"Definition of the rigid Cockle object."
struct RigidCockle <: Cockle
    brep_path::String # Path to the Cockle Brep file
    x::Float64 # x-coordinate of the center of the Cockle
    y::Float64 # y-coordinate of
    z::Float64 # z-coordinate of the center of the Cockle
    θ::Float64 # Spherical coordinate θ to manage the orientation of the Cockle
    φ::Float64 # Spherical coordinate φ to manage the orientation of the Cockle
    r::Float64 # Radius of spherical coordinate to manage the scale of the Cockle
end

"This function defines a template for the generation of marine forms and display an error message if the type is not implemented."
# function generate(::Type{T}, params::Dict{Symbol,Any}) where T <: MarineForms
#     error("Generate not implemented for type ", T)
# end

"This function generates a Sphere object based on the given parameters.
The first argument is the type of the object to generate, and the second argument
is a dictionary with the parameters of the object."
function generate(::Type{RigidSphere}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    z_range = get(params, :z_range, (0.0, 1.0))  
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
    
    # Needed to compute the tolerance
    r = rand(r_distribution)
    
    return RigidSphere(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        z_range[1] + (z_range[2] - z_range[1]) * rand(),  
        r,
        0.5 * r  # Adjust tolerance, assuming it's based on radius
    )
end

function generate(::Type{PorousSphere}, params::Dict{Symbol,Any})
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    z_range = get(params, :z_range, (0.0, 1.0))  # Add z_range for the 3D case
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
    
    # Needed to compute the tolerance
    r = rand(r_distribution)
    
    return PorousSphere(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        z_range[1] + (z_range[2] - z_range[1]) * rand(),  # Generate the z coordinate
        r,
        0.5 * r  # Adjust tolerance, assuming it's based on radius
    )
end

function generate_geometry(animal::PorousSphere)
    
    gmsh.model.occ.addSphere(animal.x, animal.y, animal.z, animal.r)     # Add a sphere
end

function generate_geometry(animal::RigidSphere, h::Float64)
    
    gmsh.model.occ.addSphere(animal.x, animal.y, animal.z, animal.r)     # Add a sphere
end

# Function to generate animals based on a definition
function generate_animals(definition)
    animals = Vector{MarineForms}()
    for (form_type, params) in definition
        for _ in 1:params[:N]
            flag = false
            while !flag
                new_animal = generate(form_type, params)
                if is_collision_free(new_animal, animals)
                    push!(animals, new_animal)
                    flag = true
                end
            end
        end
    end
    return animals
end

"This function checks if exists collision between two marine forms.
Since at this time we only have Clamps and Spheres, we can use the exterior radius
to check the collision and the entry of the function can be an abstract animal. At the 
moment we introduce new animals, we can't use this directly, but we can use the concrete types."
# function check_collisions(form1::MarineForms, form2::MarineForms)
#     dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2)
    
#     return dist < (form1.re + form2.re + max(form1.tol, form2.tol) )  # Condition to have a collision
# end

function check_collisions(form1::Sphere, form2::Sphere)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2 + (form1.z - form2.z)^2)
    
    return dist < (form1.r + form2.r + max(form1.tol, form2.tol))  # Condition to have a collision
end



# + max(form1.tol, form2.tol)
"""
This function checks if exists collision between the animal created and the previously created animals.
Due to that the array initially is filled with nothing, the function checks if the animal is not nothing,
and then checks if the animal collides with the existing animals.
"""
function is_collision_free(new_animal::MarineForms, existing_animals::Vector{MarineForms})
    for animal in existing_animals
        if check_collisions(new_animal, animal)
            return false  # collision
        end
    end
    return true  
end

end 


