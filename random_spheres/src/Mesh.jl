# Load the packages
using Random
Random.seed!(1234)
using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh
using Distributions
using Revise

includet("./MarineFormsDict.jl")
includet("Boxes.jl")
using .Boxes
using .MarineForm

include("Configuration.jl") 

# Initialize Gmsh and add a model
gmsh.initialize()
gmsh.option.setNumber("Mesh.Algorithm3D", 10) # Number 10 is the Delaunay algorithm for 3D meshing that is more robust than the frontal one at least from the point of view of Gridap. Since Gmsh was able to construct the mesh using the frontal algorithm, but Gridap was not able to handle it correctly.
gmsh.model.add("RectangleSpheres")

# Mesh size parameter (15 elements per wavelength)
h_f = c_F(ω) / (15 * f)
h_p = c_P(ω) / (15 * f)
h_spheres = h_p

# Fluid domain
fluid = add_box_to_gmsh(-L/2, t_P, -w/2, L, t_F, w, "Fluid")

# xx PMLs
fluid_right_pml_xx = add_box_to_gmsh(L/2, t_P, -w/2, d_PML, t_F, w, "fluid_right_pml_xx")
fluid_left_pml_xx = add_box_to_gmsh(-L/2 - d_PML, t_P, -w/2, d_PML, t_F, w, "fluid_left_pml_xx")

# yy PMLs
fluid_top_pml_yy = add_box_to_gmsh(-L/2, t_P + t_F, -w/2, L, d_PML, w, "fluid_top_pml_yy")

# zz PMLs
fluid_front_pml_zz = add_box_to_gmsh(-L/2, t_P, w/2, L, t_F, d_PML, "fluid_front_pml_zz")
fluid_back_pml_zz = add_box_to_gmsh(-L/2, t_P, -w/2 - d_PML, L, t_F, d_PML, "fluid_back_pml_zz")

## xz PMLs
fluid_right_front_pml_xz = add_box_to_gmsh(L/2, t_P, w/2, d_PML, t_F, d_PML, "fluid_right_front_pml_xz")
fluid_left_front_pml_xz = add_box_to_gmsh(-L/2 - d_PML, t_P, w/2, d_PML, t_F, d_PML, "fluid_left_front_pml_xz")
fluid_right_back_pml_xz = add_box_to_gmsh(L/2, t_P, -w/2 - d_PML, d_PML, t_F, d_PML, "fluid_right_back_pml_xz")
fluid_left_back_pml_xz = add_box_to_gmsh(-L/2 - d_PML, t_P, -w/2 - d_PML, d_PML, t_F, d_PML, "fluid_left_back_pml_xz")

## yz PMLs
fluid_front_top_pml_yz = add_box_to_gmsh(-L/2, t_P + t_F, w/2, L, d_PML, d_PML, "fluid_front_top_pml_yz")
fluid_back_top_pml_yz = add_box_to_gmsh(-L/2, t_P + t_F, -w/2 - d_PML, L, d_PML, d_PML, "fluid_back_top_pml_yz")

## xy PMLs
fluid_right_front_pml_xy = add_box_to_gmsh(L/2, t_P + t_F, w/2, d_PML, d_PML, -w, "fluid_right_front_pml_xy")
fluid_left_front_pml_xy = add_box_to_gmsh(-L/2 - d_PML, t_P + t_F, w/2, d_PML, d_PML, -w, "fluid_left_front_pml_xy")

## xyz PMLs
fluid_right_front_top_pml_xyz = add_box_to_gmsh(L/2, t_P + t_F, w/2, d_PML, d_PML, d_PML, "fluid_right_front_top_pml_xyz")
fluid_left_front_top_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, t_P + t_F, w/2, d_PML, d_PML, d_PML, "fluid_left_front_top_pml_xyz")
fluid_right_back_top_pml_xyz = add_box_to_gmsh(L/2, t_P + t_F, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_right_back_top_pml_xyz")
fluid_left_back_top_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, t_P + t_F, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_left_back_top_pml_xyz")

# Porous domain
porous = add_box_to_gmsh(-L/2, 0, -w/2, L, t_P, w, "Porous")

# xx porous PMLs
porous_right_pml_xx = add_box_to_gmsh(L/2, 0, -w/2, d_PML, t_P, w, "porous_right_pml_xx")
porous_left_pml_xx = add_box_to_gmsh(-L/2 - d_PML, 0, -w/2, d_PML, t_P, w, "porous_left_pml_xx")

# yy porous PMLs
porous_ground_pml_yy = add_box_to_gmsh(-L/2, -d_PML, -w/2, L, d_PML, w, "porous_ground_pml_yy")

# zz porous PMLs
porous_front_pml_zz = add_box_to_gmsh(-L/2, 0, w/2, L, t_P, d_PML, "porous_front_pml_zz")
porous_back_pml_zz = add_box_to_gmsh(-L/2, 0, -w/2 - d_PML, L, t_P, d_PML, "porous_back_pml_zz")

## xz porous PMLs
porous_right_front_pml_xz = add_box_to_gmsh(L/2, 0, w/2, d_PML, t_P, d_PML, "porous_right_front_pml_xz")
porous_left_front_pml_xz = add_box_to_gmsh(-L/2 - d_PML, 0, w/2, d_PML, t_P, d_PML, "porous_left_front_pml_xz")
porous_right_back_pml_xz = add_box_to_gmsh(L/2, 0, -w/2 - d_PML, d_PML, t_P, d_PML, "porous_right_back_pml_xz")
porous_left_back_pml_xz = add_box_to_gmsh(-L/2 - d_PML, 0, -w/2 - d_PML, d_PML, t_P, d_PML, "porous_left_back_pml_xz")

# yz porous PMLs
porous_front_bottom_pml_yz = add_box_to_gmsh(-L/2, -d_PML, w/2, L, d_PML, d_PML, "porous_front_bottom_pml_yz")
porous_back_bottom_pml_yz = add_box_to_gmsh(-L/2, -d_PML, -w/2 - d_PML, L, d_PML, d_PML, "porous_back_bottom_pml_yz")

# xy porous PMLs
porous_right_bottom_pml_xy = add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, -w, "porous_right_bottom_pml_xy")
porous_left_bottom_pml_xy = add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, -w, "porous_left_bottom_pml_xy")

# xyz porous PMLs
porous_right_front_bottom_pml_xyz = add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, d_PML, "porous_right_front_bottom_pml_xyz")
porous_left_front_bottom_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, d_PML, "porous_left_front_bottom_pml_xyz")
porous_right_back_bottom_pml_xyz = add_box_to_gmsh(L/2, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "porous_right_back_bottom_pml_xyz")
porous_left_back_bottom_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "porous_left_back_bottom_pml_xyz")

gmsh.model.occ.removeAllDuplicates()  # Clean up duplicates if needed

gmsh.model.occ.synchronize()  # Single synchronize after cleaning duplicates

# Definition of the sphere objects of the bottom (at the moment only rigid spheres)
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
    :z_range => (-w/2 + (r + 4*σ_r + tol_sphere), w/2 - (r + 4*σ_r + tol_sphere))))
]

# Generate the animals
animals = generate_animals(definition)

rigid_spheres = [animal for animal in animals if isa(animal, RigidSphere)]
porous_spheres = [animal for animal in animals if isa(animal, PorousSphere)]

v_rigid_spheres = [generate_geometry(i) for i in rigid_spheres]
v_porous_spheres = [generate_geometry(i) for i in porous_spheres]

gmsh.model.occ.synchronize()

# Perform a boolean subtraction (subtract the sphere from the box)
porous_operation = gmsh.model.occ.cut([(3, porous.tag)], [(3, i) for i in v_rigid_spheres])
porous.tag = porous_operation[1][1][2]

# Synchronize after the boolean operation
gmsh.model.occ.synchronize()

# First construct the (dim,tag) array of tuples for the porous spheres
porous_spheres_tuple = [(3, i) for i in v_porous_spheres]

# Now perform the fragmentization of the porous domain together with the porous spheres
if length(v_porous_spheres) > 0
    ov, ovv = gmsh.model.occ.fragment([(3, porous.tag)], porous_spheres_tuple)
    
    # The porous tag is now the last element of the ov array
    porous.tag = last(ov)[2]

    # Synchronize after the boolean operation
    gmsh.model.occ.synchronize()
end

# Identify each of the surfaces
surfaces = gmsh.model.occ.getEntities(2)
surface_rigid_spheres_markers = []
surface_porous_spheres_markers = []
surface_boundaries_marker = []

for surface in surfaces
    com = gmsh.model.occ.getCenterOfMass(surface[1], surface[2])
    # Boolean flag to track if the surface is inside any sphere
     inside_any_sphere = false
    
     # Check the surface against all srigid pheres
     for sphere in rigid_spheres
        # Compute the distance from the surface's center of mass to the center of the sphere
        dist = sqrt((com[1] - sphere.x)^2 + (com[2] - sphere.y)^2 + (com[3] - sphere.z)^2)
        # Check if the surface is inside the sphere, considering tolerance
        if dist < (sphere.r)
            push!(surface_rigid_spheres_markers, surface[2])
            inside_any_sphere = true
            break  # Exit the loop once inside any sphere
        end
    end
    
    # Check the surface against all porous spheres
    for sphere in porous_spheres
        # Compute the distance from the surface's center of mass to the center of the sphere
        dist = sqrt((com[1] - sphere.x)^2 + (com[2] - sphere.y)^2 + (com[3] - sphere.z)^2)
        # Check if the surface is inside the sphere, considering tolerance
        if dist < (sphere.r)
            push!(surface_porous_spheres_markers, surface[2])
            inside_any_sphere = true
            break  # Exit the loop once inside any sphere
        end
    end
    
    # If the surface is not inside any sphere, push to the boundary marker that 
    if !inside_any_sphere && ((abs(com[1]) - L/2 - 0.99 * d_PML) > 0 || (abs(com[3]) - w/2 - 0.99 * d_PML) > 0 || (com[2] + 0.99 * d_PML) < 0 || (com[2] - t_P - t_F - 0.99 * d_PML) > 0)
        push!(surface_boundaries_marker, surface[2])
    end
end

# All the volumes corresponding to the fluid domain
all_fluid_domains = [fluid, fluid_right_pml_xx, fluid_left_pml_xx, fluid_top_pml_yy, fluid_front_pml_zz, fluid_back_pml_zz, fluid_right_front_pml_xz, fluid_left_front_pml_xz, fluid_right_back_pml_xz, fluid_left_back_pml_xz, fluid_front_top_pml_yz, fluid_back_top_pml_yz, fluid_right_front_pml_xy, fluid_left_front_pml_xy, fluid_right_front_top_pml_xyz, fluid_left_front_top_pml_xyz, fluid_right_back_top_pml_xyz, fluid_left_back_top_pml_xyz]

# All the volumes corresponding to the porous domain
all_porous_domains = [porous, porous_right_pml_xx, porous_left_pml_xx, porous_ground_pml_yy, porous_front_pml_zz, porous_back_pml_zz, porous_right_front_pml_xz, porous_left_front_pml_xz, porous_right_back_pml_xz, porous_left_back_pml_xz, porous_front_bottom_pml_yz, porous_back_bottom_pml_yz, porous_right_bottom_pml_xy, porous_left_bottom_pml_xy, porous_right_front_bottom_pml_xyz, porous_left_front_bottom_pml_xyz, porous_right_back_bottom_pml_xyz, porous_left_back_bottom_pml_xyz]

# 3D groups for the resulting fragments
pml_groups = Dict(
    "fluid_pml_xx" => [fluid_right_pml_xx.tag, fluid_left_pml_xx.tag],
    "fluid_pml_yy" => [fluid_top_pml_yy.tag],
    "fluid_pml_zz" => [fluid_front_pml_zz.tag, fluid_back_pml_zz.tag],
    
    "fluid_pml_yz" => [
        fluid_front_top_pml_yz.tag, fluid_back_top_pml_yz.tag,
    ],
    "fluid_pml_xz" => [
        fluid_right_front_pml_xz.tag, fluid_left_front_pml_xz.tag,
        fluid_right_back_pml_xz.tag, fluid_left_back_pml_xz.tag
    ],
    "fluid_pml_xy" => [
        fluid_right_front_pml_xy.tag, fluid_left_front_pml_xy.tag
    ],
    
    "fluid_pml_xyz" => [
        fluid_right_front_top_pml_xyz.tag, fluid_left_front_top_pml_xyz.tag,
        fluid_right_back_top_pml_xyz.tag, fluid_left_back_top_pml_xyz.tag
    ],
    
    "porous_pml_xx" => [porous_right_pml_xx.tag, porous_left_pml_xx.tag],
    "porous_pml_yy" => [porous_ground_pml_yy.tag],
    "porous_pml_zz" => [porous_front_pml_zz.tag, porous_back_pml_zz.tag],

    "porous_pml_xz" => [
        porous_right_front_pml_xz.tag, porous_left_front_pml_xz.tag,
        porous_right_back_pml_xz.tag, porous_left_back_pml_xz.tag
    ],
    "porous_pml_yz" => [
        porous_front_bottom_pml_yz.tag, porous_back_bottom_pml_yz.tag,
    ],
    "porous_pml_xy" => [
        porous_right_bottom_pml_xy.tag, porous_left_bottom_pml_xy.tag,
    ],
    "porous_pml_xyz" => [
        porous_right_front_bottom_pml_xyz.tag, porous_left_front_bottom_pml_xyz.tag,
        porous_right_back_bottom_pml_xyz.tag, porous_left_back_bottom_pml_xyz.tag,
    ]
)

# Set the physical groups for the fluid domains
fluid_tag = gmsh.model.addPhysicalGroup(3, [fluid.tag])
fluid_pml_xx_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xx"])
fluid_pml_yy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yy"])
fluid_pml_zz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_zz"])
fluid_pml_xy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xy"])
fluid_pml_yz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yz"])
fluid_pml_xz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xz"])
fluid_pml_xyz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xyz"])

# Set the physical groups for the porous domains
porous_tag = gmsh.model.addPhysicalGroup(3, [porous.tag])
porous_pml_xx_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xx"])
porous_pml_yy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_yy"])
porous_pml_zz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_zz"])
porous_pml_xy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xy"])
porous_pml_yz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_yz"])
porous_pml_xz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xz"])
porous_pml_xyz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xyz"])

# Set the physical groups for the porous spheres
porous_spheres_tag = gmsh.model.addPhysicalGroup(3, v_porous_spheres)

# Set the physical groups for the surfaces
surface_rigid_spheres_markers = gmsh.model.addPhysicalGroup(2, surface_rigid_spheres_markers)
surface_porous_spheres_markers = gmsh.model.addPhysicalGroup(2, surface_porous_spheres_markers)
surface_boundaries_marker = gmsh.model.addPhysicalGroup(2, surface_boundaries_marker)

# Volumes of the fluid  and pml fluid domains
gmsh.model.setPhysicalName(3, fluid_tag, "Fluid")
gmsh.model.setPhysicalName(3, fluid_pml_xx_tag, "Fluid PML_xx")
gmsh.model.setPhysicalName(3, fluid_pml_yy_tag, "Fluid PML_yy")
gmsh.model.setPhysicalName(3, fluid_pml_zz_tag, "Fluid PML_zz")
gmsh.model.setPhysicalName(3, fluid_pml_xy_tag, "Fluid PML_xy")
gmsh.model.setPhysicalName(3, fluid_pml_yz_tag, "Fluid PML_yz")
gmsh.model.setPhysicalName(3, fluid_pml_xz_tag, "Fluid PML_xz")
gmsh.model.setPhysicalName(3, fluid_pml_xyz_tag, "Fluid PML_xyz")

# Volumes of the porous and pml porous domains
gmsh.model.setPhysicalName(3, porous_tag, "Porous")
gmsh.model.setPhysicalName(3, porous_pml_xx_tag, "Porous PML_xx")
gmsh.model.setPhysicalName(3, porous_pml_yy_tag, "Porous PML_yy")
gmsh.model.setPhysicalName(3, porous_pml_zz_tag, "Porous PML_zz")
gmsh.model.setPhysicalName(3, porous_pml_xy_tag, "Porous PML_xy")
gmsh.model.setPhysicalName(3, porous_pml_yz_tag, "Porous PML_yz")
gmsh.model.setPhysicalName(3, porous_pml_xz_tag, "Porous PML_xz")
gmsh.model.setPhysicalName(3, porous_pml_xyz_tag, "Porous PML_xyz")

# Volumes of the porous spheres
gmsh.model.setPhysicalName(3, porous_spheres_tag, "Porous Spheres")

# Set 2D-physical names for the physical groups
gmsh.model.setPhysicalName(2, surface_rigid_spheres_markers, "RigidSpheresBoundary")
gmsh.model.setPhysicalName(2, surface_porous_spheres_markers, "PorousSpheresBoundary")
gmsh.model.setPhysicalName(2, surface_boundaries_marker, "Boundaries")

# Set the mesh size
ov_f = gmsh.model.getBoundary([(3, fluid_domain.tag) for fluid_domain in all_fluid_domains], false, false, true)
ov_p = gmsh.model.getBoundary([(3, porous_domain.tag) for porous_domain in all_porous_domains], false, false, true)
gmsh.model.mesh.setSize(ov_f, h_f)
gmsh.model.mesh.setSize(ov_p, h_p)
gmsh.model.mesh.setSize(gmsh.model.getBoundary([(2, i) for i in 7:N_rigid_spheres+6], true, true, true), h_spheres)

# Synchronize to register the physical group
gmsh.model.occ.synchronize()

# Generate a 3D mesh
gmsh.model.mesh.generate(3)

# Save the resulting mesh to a file
gmsh.write("./data/multiple_spheres.msh")

# Finalize Gmsh session
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/multiple_spheres.msh")

# Write the mesh to a vtk file
writevtk(model,"./results/multiple_spheres")
