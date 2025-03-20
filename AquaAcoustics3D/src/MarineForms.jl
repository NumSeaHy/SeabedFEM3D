module MarineForm

# Loading the packages
using Distributions
using PyCall
using BoundingSphere
using Gmsh

# Export the functions of this module that are going to be used in other modules
export generate_animals, PorousSphere, RigidSphere, RigidCockle, RigidQueenScallop, RigidScallop, generate_geometry

mutable struct CommonGeometry
    brep_path::String
    x::Float64
    y::Float64
    z::Float64
    scaled_radius::Float64
    by_default_radius::Float64
    α::Float64
    β::Float64
    γ::Float64
    bounding_box::NTuple{6, Float64}  
    evolving_sphere_center::Vector{Float64}
    evolving_sphere_radius::Float64
end

struct CommonSphere
    x::Float64
    y::Float64
    z::Float64
    r::Float64
    tol::Float64
end

# Base abstract type for all marine animal geometries.
abstract type MarineAnimalGeometry end

# Define an abstract subtype for scanned geometries.
abstract type ScannedMarineAnimalGeometry <: MarineAnimalGeometry end

# Define all the scanned geometries as subtypes of ScannedMarineAnimalGeometry.
abstract type Cockle <: ScannedMarineAnimalGeometry end

abstract type QueenScallop <: ScannedMarineAnimalGeometry end

abstract type Scallop <: ScannedMarineAnimalGeometry end

# Define an abstract subtype for synthetic (gmsh-based) geometries.
abstract type SyntheticMarineAnimalGeometry <: MarineAnimalGeometry end

# Define all the synthetic geometries as subtypes of SyntheticMarineAnimalGeometry.
abstract type Sphere <: SyntheticMarineAnimalGeometry end

# Global cache: maps a brep_path to a tuple of (center, radius)
const EVOLVING_SPHERE_CACHE = Dict{String, Tuple{Vector{Float64}, Float64}}()

# This function checks if the value is cached; if not, computes and caches it.
function get_cached_center_and_radius(brep_path::String)
    # Check is the previous bounding sphere computation is cached
    if haskey(EVOLVING_SPHERE_CACHE, brep_path)
        return EVOLVING_SPHERE_CACHE[brep_path]
    else
        # Call the expensive computation
        center, radius = get_evolving_sphere(brep_path)
        EVOLVING_SPHERE_CACHE[brep_path] = (center, radius)
        return center, radius
    end
end

struct RigidCockle <: Cockle
    geometry::CommonGeometry
    function RigidCockle(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box)
        center, radius = get_cached_center_and_radius(brep_path)
        geo = CommonGeometry(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box, center, radius)
        new(geo)
    end
end

struct RigidQueenScallop <: QueenScallop
    geometry::CommonGeometry
    function RigidQueenScallop(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box)
        center, radius = get_cached_center_and_radius(brep_path)
        geo = CommonGeometry(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box, center, radius)
        new(geo)
    end
end

struct RigidScallop <: Scallop
    geometry::CommonGeometry
    function RigidScallop(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box)
        center, radius = get_cached_center_and_radius(brep_path)
        geo = CommonGeometry(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box, center, radius)
        new(geo)
    end
end


"Definition of the Porous Sphere object."
struct PorousSphere <: Sphere
    geometry::CommonSphere
    function PorousSphere(x, y, z, r, tol)
        geo = CommonSphere(x, y, z, r, tol)
        new(geo)
    end
end

"Definition of the Rigid Sphere object."
struct RigidSphere <: Sphere
    geometry::CommonSphere
    function RigidSphere(x, y, z, r, tol)
        geo = CommonSphere(x, y, z, r, tol)
        new(geo)
    end
end

"Function to apply yaw, pitch and roll rotations to the object"
function apply_rotation(α::Float64, β::Float64, γ::Float64)
    return [
        cos(α) -sin(α) 0
        sin(α) cos(α) 0
        0 0 1
    ] * [
        cos(β) 0 sin(β)
        0 1 0
        -sin(β) 0 cos(β)
    ] * [
        1 0 0
        0 cos(γ) -sin(γ)
        0 sin(γ) cos(γ)
    ]
end

"This function defines a template for the generation of marine forms and display an error message if the type is not implemented."
# function generate(::Type{T}, params::Dict{Symbol,Any}) where T <: MarineForms
#     error("Generate not implemented for type ", T)
# end

"This function generates a Sphere object based on the given parameters.
The first argument is the type of the object to generate, and the second argument
is a dictionary with the parameters of the object."
function generate(::Type{T}, params::Dict{Symbol,Any}) where {T <: Sphere}
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    z_range = get(params, :z_range, (0.0, 1.0)) 
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))

    # Needed to compute the tolerance
    r = rand(r_distribution)

    return T(
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        z_range[1] + (z_range[2] - z_range[1]) * rand(), 
        r,
        0.5 * r  # Adjust tolerance, assuming it's based on radius
    )
    
end

function  generate(::Type{T}, params::Dict{Symbol,Any}) where {T <: ScannedMarineAnimalGeometry}
     # Read the desired cockle type from params. For example, it could be :open or :closed.
     type = get(params, :type, :open)
     # Get the name of the object, needed to load the correct brep file.
     name = get(params, :name, "NoName")  

     # Choose the brep_path based on the type.
     brep_path = type == :open ? "./geometry/"*lowercase(name)*"_geometries/Open"*name*".brep" : "./geometry/"*lowercase(name)*"_geometries/Closed"*(name)*".brep"

    if isfile(brep_path)
        nothing
    else
        error("Animal type not found: $name in the folder $brep_path")
    end
    
     # Retrieve spatial and distribution parameters.
     x_range        = get(params, :x_range, (0.0, 1.0))
     y_range        = get(params, :y_range, (0.0, 1.0))
     z_range        = get(params, :z_range, (0.0, 1.0))
     r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
     by_default_radius = get(params, :by_default_radius, 1.0)
     α_range        = get(params, :α_range, (0.0, 0.0))
     β_range        = get(params, :β_range, (0.0, 0.0))
     γ_range        = get(params, :γ_range, (0.0, 0.0))
 
     # Generate the random values (deterministic if RNG is seeded)
     x = x_range[1] + (x_range[2] - x_range[1]) * rand()
     y = y_range[1] + (y_range[2] - y_range[1]) * rand()
     z = z_range[1] + (z_range[2] - z_range[1]) * rand()
     scaled_radius = rand(r_distribution)
     α = α_range[1] + (α_range[2] - α_range[1]) * rand()
     β = β_range[1] + (β_range[2] - β_range[1]) * rand()
     γ = γ_range[1] + (γ_range[2] - γ_range[1]) * rand()
 
    # By default bounding box, one the gmsh object is generated this value is updated (reason why is a mutable object)
     bounding_box = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
 
     # Use the appropriate constructor based on cockle_type.
     return T(brep_path, x, y, z, scaled_radius, by_default_radius, α, β, γ, bounding_box)

end

"Geometry generation functions"
function generate_geometry(animal::ScannedMarineAnimalGeometry)
    
    # Import the solid from the brep file
    object = gmsh.model.occ.importShapes(animal.geometry.brep_path)
    
    # Since the CAD objet is imported in milimeters from FreeCAD is a must scale it to meters
    scale_factor = 0.001
    
    # Use the by_default_radius to scale the object to the desired size
    scale_factor *= animal.geometry.scaled_radius / animal.geometry.by_default_radius

    # Apply the affine transformation for scaling rotation[] and translation
    scaling_matrix = [
        scale_factor, 0, 0, 0,  # X scale
        0, scale_factor, 0, 0,  # Y scale
        0, 0, scale_factor, 0,  # Z scale            # Homogeneous component for affine transformations
    ]

    # Apply the scaling transformation
    gmsh.model.occ.affineTransform([(object[1][1], object[1][2])], scaling_matrix)

    # Apply the rotation. Constructing the rotation matrix based on yaw, pth, roll angles of the object
    rotation_matrix = apply_rotation(animal.geometry.α, animal.geometry.β, animal.geometry.γ)
    
    # Based on the previous rotation matrix, we need to create a 1D array to pass it to the gmsh API
    rotation_matrix_gmsh = [
        rotation_matrix[1, 1], rotation_matrix[1, 2], rotation_matrix[1, 3], 0,  # X axis
        rotation_matrix[2, 1], rotation_matrix[2, 2], rotation_matrix[2, 3], 0,  # Y axis
        rotation_matrix[3, 1], rotation_matrix[3, 2], rotation_matrix[3, 3], 0,  # Z axis
    ]

    # Apply the rotation by using the affine transformation
    gmsh.model.occ.affineTransform([(object[1][1], object[1][2])], rotation_matrix_gmsh)

    # Get the bounding box of the object to make the translation
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

    # Translate the current object based on the desired animal.geometry position
    dx = animal.geometry.x - (xmax + xmin) / 2
    dy = animal.geometry.y - (ymax + ymin) / 2
    dz = animal.geometry.z - (zmax + zmin) / 2

    # Apply the translation
    gmsh.model.occ.translate([(object[1][1], object[1][2])], dx, dy, dz)

    # Get the bounding box of the object after all the operations to update the bounding box of the object and then detect the object
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

    # Update the bounding box of the object
    animal.geometry.bounding_box = (xmin, ymin, zmin, xmax, ymax, zmax)

    return object[1][2]

end

function generate_geometry(animal::Sphere)
    gmsh.model.occ.addSphere(animal.x, animal.y, animal.z, animal.r)     # Add a sphere
end

"Function to generate animals based on a definition"
function generate_animals(definition, max_it)
    animals = Vector{MarineAnimalGeometry}()
    iterations = 0
    for (form_type, params) in definition
        for _ in 1:params[:N]
            flag = false
            while !flag && iterations < max_it
                new_animal = generate(form_type, params)
                if is_collision_free(new_animal, animals)
                    push!(animals, new_animal)
                    flag = true
                end
                iterations += 1
            end
            if iterations >= max_it
                error("Max iterations reached, maybe the objects are too big or the domain is too small")
            end
        end
    end
    return animals
end

"Overlapping check functions"
function check_collisions(form1::Sphere, form2::Sphere)
    dist = sqrt((form1.x - form2.x)^2 + (form1.y - form2.y)^2 + (form1.z - form2.z)^2)
    
    return dist < (form1.r + form2.r + max(form1.tol, form2.tol))  # Condition to have a collision
end

function check_collisions(form1::ScannedMarineAnimalGeometry, form2::ScannedMarineAnimalGeometry)
    # Translate the center of the form1 to the Cockle desired position
    center1 = [form1.geometry.evolving_sphere_center[1] + (form1.geometry.x-form1.geometry.evolving_sphere_center[1]), form1.geometry.evolving_sphere_center[2] +  (form1.geometry.y-form1.geometry.evolving_sphere_center[2]), form1.geometry.evolving_sphere_center[3] + (form1.geometry.z-form1.geometry.evolving_sphere_center[3])]
    # Scale the radius of the form1 to the desired size taking into account the desire radius of the Cockle
    radius1 = form1.geometry.evolving_sphere_radius * form1.geometry.scaled_radius/form1.geometry.by_default_radius 

    # Translate the center of the form2 to the Cockle desired position
    center2 = [form2.geometry.evolving_sphere_center[1] + (form2.geometry.x-form2.geometry.evolving_sphere_center[1]), form2.geometry.evolving_sphere_center[2] +  (form2.geometry.y-form2.geometry.evolving_sphere_center[2]), form2.geometry.evolving_sphere_center[3] + (form2.geometry.z-form2.geometry.evolving_sphere_center[3])]
    # Scale the radius of the form2 to the desired size taking into account the desire radius of the Cockle
    radius2 = form2.geometry.evolving_sphere_radius * form2.geometry.scaled_radius/form2.geometry.by_default_radius

    return sqrt((center1[1] - center2[1])^2 + (center1[2] - center2[2])^2 + (center1[3] - center2[3])^2) < (radius1 + radius2 + max(radius1/2, radius2/2))  # 
    
end

function check_collisions(form1::ScannedMarineAnimalGeometry, form2::Sphere)
    evolving_sphere_center = [form1.geometry.evolving_sphere_radius[1] + (form1.geometry.x-form1.geometry.evolving_sphere_radius[1]), form1.geometry.evolving_sphere_radius[2] +  (form1.geometry.y-form1.geometry.evolving_sphere_radius[2]), form1.geometry.evolving_sphere_radius[3] + (form1.geometry.z-form1.geometry.evolving_sphere_radius[3])]
    
    animal_sphere_radius = form1.geometry.evolving_sphere_radius * form1.geometry.scaled_radius/form1.geometry.by_default_radius 

    return sqrt((form2.x - evolving_sphere_center[1])^2 + (form2.y - evolving_sphere_center[2])^2 + (form2.z - evolving_sphere_center[3])^2) < (form2.r + animal_sphere_radius)  # Condition to have a collision 
end

"""
This function checks if exists collision between the animal created and the previously created animals.
Due to that the array initially is filled with nothing, the function checks if the animal is not nothing,
and then checks if the animal collides with the existing animals.
"""
function is_collision_free(new_animal::MarineAnimalGeometry, existing_animals::Vector{MarineAnimalGeometry})
    for animal in existing_animals
        if check_collisions(new_animal, animal)
            return false  # collision
        end
    end
    return true  
end

function get_evolving_sphere(brep_path::String)
    
    # Import the solid from the brep file using OpenCascade Python API
    breptools = pyimport("OCC.Core.BRepTools")
    occ = pyimport("OCC.Core")

    # Create objects to read the BRep file
    builder = occ.BRep.BRep_Builder()
    shape = occ.TopoDS.TopoDS_Shape()
    breptools.breptools_Read(shape, brep_path, builder)

    # Import required OCC modules
    TopExp = pyimport("OCC.Core.TopExp")
    TopAbs = pyimport("OCC.Core.TopAbs")
    TopoDS = pyimport("OCC.Core.TopoDS")
    BRep = pyimport("OCC.Core.BRep")

    # Create an explorer for vertices in the shape
    explorer = TopExp.TopExp_Explorer(shape, TopAbs.TopAbs_VERTEX)

    # Initialize an empty Julia array to hold the points
    points = []

    # Iterate over all vertices in the shape
    while explorer.More()
        # Convert the current item to a vertex
        vertex = TopoDS.topods_Vertex(explorer.Current())
        # Get the point coordinates of the vertex
        pnt = BRep.BRep_Tool.Pnt(vertex)
        # Append the (x, y, z) tuple to the points array
        push!(points, [pnt.X(), pnt.Y(), pnt.Z()])
        # Move to the next vertex
        explorer.Next()
    end

    points_correct_type_mm = [v for v in points] .* 1e-3 # Convert the points to meters and to the data format expected by the  BoundingSphere package
    
    # Compute the bounding sphere of the points
    center, radius = boundingsphere(points_correct_type_mm)

    return center, radius
    
end

end 


