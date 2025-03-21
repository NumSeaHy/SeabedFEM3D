module MarineForm

# Loading the packages
using Distributions
using PyCall
using BoundingSphere
using Gmsh

# Export the functions of this module that are going to be used in other modules
export generate_animals, PorousSphere, RigidSphere, RigidCockle, generate_geometry

"Definition of the abstract type MarineForms."
abstract type MarineForms end

"Definition of the abstract type Sphere."
abstract type Sphere <: MarineForms end

"Definition of the abstract type Cockle"
abstract type Cockle <: MarineForms end

struct EnvolvingSphere
    path::String
    center::Tuple{Float64, Float64, Float64}
    radius::Float64
end

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
mutable struct RigidCockle <: Cockle
    brep_path::String # Path to the Cockle Brep file
    x::Float64 # x-coordinate of the center of the Cockle
    y::Float64 # y-coordinate of the center of the Cockle
    z::Float64 # z-coordinate of the center of the Cockle
    scaled_radius::Float64 # Desired radius of the cockle
    by_default_radius::Float64 # Default radius of the Cockle measure in Real world (Let us see if is or not necessary.........)
    α::Float64 # Yaw angle to manage the orientation of the Cockle
    β::Float64 # Pitch angle to manage the orientation of the Cockle
    γ::Float64 # Roll angle to manage the orientation of the Cockle
    bounding_box::Tuple{Float64, Float64, Float64, Float64, Float64, Float64} # Bounding box of the Cockle
end

"Auxiliary functions to define the yaw, pitch and roll rotation angles"
function Rz(α)
    return [
        cos(α) -sin(α) 0
        sin(α) cos(α) 0
        0 0 1
    ]
end

function Ry(β)
    return [
        cos(β) 0 sin(β)
        0 1 0
        -sin(β) 0 cos(β)
    ]
end

function Rx(γ)
    return [
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

function generate(::Type{RigidCockle}, params::Dict{Symbol, Any})
    brep_path = get(params, :brep_path, "")
    x_range = get(params, :x_range, (0.0, 1.0))
    y_range = get(params, :y_range, (0.0, 1.0))
    z_range = get(params, :z_range, (0.0, 1.0))  
    r_distribution = get(params, :r_distribution, Normal(1.0, 0.2))
    by_default_radius = get(params, :by_default_radius, 1.0)
    α_range = get(params, :α_range, (0.0, 0.0))
    β_range = get(params, :β_range, (0.0, 0.0))
    γ_range = get(params, :γ_range, (0.0, 0.0))
    

    return RigidCockle(
        brep_path,
        x_range[1] + (x_range[2] - x_range[1]) * rand(), 
        y_range[1] + (y_range[2] - y_range[1]) * rand(),
        z_range[1] + (z_range[2] - z_range[1]) * rand(), 
        rand(r_distribution),
        by_default_radius, 
        α_range[1] + (α_range[2] - α_range[1]) * rand(),  # TODO Yaw angle
        β_range[1] + (β_range[2] - β_range[1]) * rand(), # TODO Pitch angle
        γ_range[1] + (γ_range[2] - γ_range[1]) * rand(),  # TODO Roll angle
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)  # Bounding box
    )
    
end


"Geometry generation functions"
function generate_geometry(animal::RigidCockle)
    
    # Import the solid from the brep file
    object = gmsh.model.occ.importShapes(animal.brep_path)
    
    # Since the CAD objet is imported in milimeters from FreeCAD is a must scale it to meters
    scale_factor = 0.001
    
    # Use the by_default_radius to scale the object to the desired size
    scale_factor *= animal.scaled_radius / animal.by_default_radius

    # Apply the affine transformation for scaling rotation[] and translation
    scaling_matrix = [
        scale_factor, 0, 0, 0,  # X scale
        0, scale_factor, 0, 0,  # Y scale
        0, 0, scale_factor, 0,  # Z scale            # Homogeneous component for affine transformations
    ]

    # Apply the scaling transformation
    gmsh.model.occ.affineTransform([(object[1][1], object[1][2])], scaling_matrix)

    # Get the bounding box of the object to make the translation
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])
    
    # Translate the current object based on the desired animal position
    dx = animal.x - (xmax + xmin) / 2
    dy = animal.y - (ymax + ymin) / 2
    dz = animal.z - (zmax + zmin) / 2

    # Apply the translation
    gmsh.model.occ.translate([(object[1][1], object[1][2])], dx, dy, dz)

    # Apply the rotation. Constructing the rotation matrix based on yaw, pth, roll angles of the object
    rotation_matrix = Rz(animal.α) * Ry(animal.β) * Rx(animal.γ)
    
    # Based on the previous rotation matrix, we need to create a 1D array to pass it to the gmsh API
    rotation_matrix_gmsh = [
        rotation_matrix[1, 1], rotation_matrix[1, 2], rotation_matrix[1, 3], 0,  # X axis
        rotation_matrix[2, 1], rotation_matrix[2, 2], rotation_matrix[2, 3], 0,  # Y axis
        rotation_matrix[3, 1], rotation_matrix[3, 2], rotation_matrix[3, 3], 0,  # Z axis
    ]

    # Apply the rotation by using the affine transformation
    gmsh.model.occ.affineTransform([(object[1][1], object[1][2])], rotation_matrix_gmsh)

    # Get the bounding box of the object after all the operations to update the bounding box of the object and then detect the object
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

    # Update the bounding box of the object
    animal.bounding_box = (xmin, ymin, zmin, xmax, ymax, zmax)

    return object[1][2]

end

function generate_geometry(animal::PorousSphere)
    
    gmsh.model.occ.addSphere(animal.x, animal.y, animal.z, animal.r)     # Add a sphere
end

function generate_geometry(animal::RigidSphere)
    
    gmsh.model.occ.addSphere(animal.x, animal.y, animal.z, animal.r)     # Add a sphere
end

# Function to generate animals based on a definition
function generate_animals(definition, max_it)
    animals = Vector{MarineForms}()
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

function check_collisions(form1::Cockle, form2::Cockle)
    center1, radius1 = get_evolving_sphere(form1)
    center2, radius2 = get_evolving_sphere(form2)

    return sqrt((center1[1] - center2[1])^2 + (center1[2] - center2[2])^2 + (center1[3] - center2[3])^2) < (radius1 + radius2)  # Condition to have a collision 
end


function check_collisions(form1::Cockle, form2::Sphere)
    cockle_sphere_center, cockle_sphere_radius = get_evolving_sphere(form1)

    return sqrt((form2.x - cockle_sphere_center[1])^2 + (form2.y - cockle_sphere_center[2])^2 + (form2.z - cockle_sphere_center[3])^2) < (form2.r + cockle_sphere_radius)  # Condition to have a collision 
end


function check_collisions(form1::Sphere, form2::Cockle)
    cockle_sphere_center, cockle_sphere_radius = get_evolving_sphere(form2)
    
    return sqrt((form1.x - cockle_sphere_center[1])^2 + (form1.y - cockle_sphere_center[2])^2 + (form1.z - cockle_sphere_center[3])^2) < (form1.r + cockle_sphere_radius)  # Condition to have a collision 
    
end

"Auxiliary functions composing the module"

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


function get_evolving_sphere(form::RigidCockle)
    
    # Import the solid from the brep file using OpenCascade Python API
    breptools = pyimport("OCC.Core.BRepTools")
    occ = pyimport("OCC.Core")

    # Create objects to read the BRep file
    builder = occ.BRep.BRep_Builder()
    shape = occ.TopoDS.TopoDS_Shape()
    breptools.breptools_Read(shape, form.brep_path, builder)

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

    # Translate the center to the Cockle desired position
    center = [center[1] + (form.x-center[1]), center[2] +  (form.y-center[2]), center[3] + (form.z-center[3])]
    
    # Scale the sphere to the desired size taking into account the desire radius of the Cockle
    scale_factor = form.r / radius
    radius *= scale_factor
    
    return center, radius 

end

function set_brep_evolving_sphere(brep_path::String)
     
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
 
     # Translate the center to the Cockle desired position
     center = [center[1] + (form.x-center[1]), center[2] +  (form.y-center[2]), center[3] + (form.z-center[3])]
     
     # Scale the sphere to the desired size taking into account the desire radius of the Cockle
     scale_factor = form.r / radius
     radius *= scale_factor
     
     return center, radius 
 
end


end 


