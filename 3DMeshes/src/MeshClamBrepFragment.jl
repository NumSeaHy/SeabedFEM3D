using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh
using Revise

include("../src/Configuration.jl")
includet("Boxes.jl")
using .Boxes

# Initialize Gmsh and add a model
gmsh.initialize()
gmsh.option.setNumber("Mesh.Algorithm3D", 10) # Number 10 is the Delaunay algorithm for 3D meshing that is more robust than the frontal one
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 3) 
gmsh.model.add("Rectangle_and_Sphere")


# Mesh size parameter (ensure variables c and f are defined in Configuration.jl)
h1 = 0.04

# Define the fluid domain
fluid = add_box_to_gmsh(-L/2, 0, -w/2, L, H, w, "Fluid")

# xx PMLs
fluid_right_pml_xx = add_box_to_gmsh(L/2, 0, -w/2, d_PML, H, w, "fluid_right_pml_xx")
fluid_left_pml_xx = add_box_to_gmsh(-L/2 - d_PML, 0, -w/2, d_PML, H, w, "fluid_left_pml_xx")

# yy PMLs
fluid_top_pml_yy = add_box_to_gmsh(-L/2, H, -w/2, L, d_PML, w, "fluid_top_pml_yy")
fluid_bottom_pml_yy = add_box_to_gmsh(-L/2, -d_PML, -w/2, L, d_PML, w, "fluid_bottom_pml_yy")

# zz PMLs
fluid_front_pml_zz = add_box_to_gmsh(-L/2, 0, w/2, L, H, d_PML, "fluid_front_pml_zz")
fluid_back_pml_zz = add_box_to_gmsh(-L/2, 0, -w/2 - d_PML, L, H, d_PML, "fluid_back_pml_zz")

## xz PMLs
fluid_right_front_pml_xz = add_box_to_gmsh(L/2, 0, w/2, d_PML, H, d_PML, "fluid_right_front_pml_xz")
fluid_left_front_pml_xz = add_box_to_gmsh(-L/2 - d_PML, 0, w/2, d_PML, H, d_PML, "fluid_left_front_pml_xz")
fluid_right_back_pml_xz = add_box_to_gmsh(L/2, 0, -w/2 - d_PML, d_PML, H, d_PML, "fluid_right_back_pml_xz")
fluid_left_back_pml_xz = add_box_to_gmsh(-L/2 - d_PML, 0, -w/2 - d_PML, d_PML, H, d_PML, "fluid_left_back_pml_xz")

## yz PMLs
fluid_front_top_pml_yz = add_box_to_gmsh(-L/2, H, w/2, L, d_PML, d_PML, "fluid_front_top_pml_yz")
fluid_back_top_pml_yz = add_box_to_gmsh(-L/2, H, -w/2 - d_PML, L, d_PML, d_PML, "fluid_back_top_pml_yz")
fluid_front_ground = add_box_to_gmsh(-L/2, -d_PML, w/2, L, d_PML, d_PML, "fluid_front_ground")
fluid_back_ground = add_box_to_gmsh(-L/2, -d_PML, -w/2 - d_PML, L, d_PML, d_PML, "fluid_back_ground")

## xy PMLs
fluid_right_front_pml_xy = add_box_to_gmsh(L/2, H, w/2, d_PML, d_PML, -w, "fluid_right_front_pml_xy")
fluid_left_front_pml_xy = add_box_to_gmsh(-L/2 - d_PML, H, w/2, d_PML, d_PML, -w, "fluid_left_front_pml_xy")
fluid_left_ground_pml_xy = add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, -w, "fluid_left_ground_pml_xy")
fluid_right_ground_pml_xy = add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, -w, "fluid_right_ground_pml_xy")


## xyz PMLs
fluid_right_front_top_pml_xyz = add_box_to_gmsh(L/2, H, w/2, d_PML, d_PML, d_PML, "fluid_right_front_top_pml_xyz")
fluid_left_front_top_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, H, w/2, d_PML, d_PML, d_PML, "fluid_left_front_top_pml_xyz")
fluid_right_back_top_pml_xyz = add_box_to_gmsh(L/2, H, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_right_back_top_pml_xyz")
fluid_left_back_top_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, H, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_left_back_top_pml_xyz")
fluid_right_front_ground_pml_xyz = add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, d_PML, "fluid_right_front_ground_pml_xyz")
fluid_left_front_ground_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, d_PML, "fluid_left_front_ground_pml_xyz")
fluid_right_back_ground_pml_xyz = add_box_to_gmsh(L/2, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_right_back_ground_pml_xyz")
fluid_left_back_ground_pml_xyz = add_box_to_gmsh(-L/2 - d_PML, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "fluid_left_back_ground_pml_xyz")

gmsh.model.occ.removeAllDuplicates()  # Clean up duplicates if needed

gmsh.model.occ.synchronize()  # Single synchronize after loop

# Open the .brep geometry file
object = gmsh.model.occ.importShapes("./ClamGeometries/EntireClam.brep")
# # Since the geometry is imported in milimeters from FreeCAD is a must scale it to meters
gmsh.model.occ.synchronize()

# Translate the geometry to the center of the domain
# # Define the scaling factor (0.001 to convert from millimeters to meters)
scale_factor = 0.001

# Apply the affine transformation for scaling rotation[] and translation
scaling_matrix = [
    scale_factor, 0, 0, 0,  # X scale
    0, scale_factor, 0, 0,  # Y scale
    0, 0, scale_factor, 0,  # Z scale            # Homogeneous component for affine transformations
]

# # Apply scaling to all entities
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

# Apply the scaling transformation
gmsh.model.occ.affineTransform([(object[1][1], object[1][2])], scaling_matrix)

# Compute the new bounding box
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

# Translate the geometry to the center of the domain
dx = (-xmax-xmin)/2
dy = 1/2*(H - ymax - ymin)
dz = (-zmax-zmin)/2
gmsh.model.occ.translate([(3, object[1][2])], dx, 0, 0)
gmsh.model.occ.translate([(3, object[1][2])], 0, dy, 0)
gmsh.model.occ.translate([(3, object[1][2])], 0, 0, -dz)


# Compute the new bounding box
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

# Perform a boolean subtraction (subtract the scattering object from the porous domain)

# porous.tag = last(ov)[2]
# porous.tag = ov[1end-1][2]

# Perform a boolean subtraction (subtract the sphere from the box)
ov, ovv = gmsh.model.occ.fragment([(3, fluid.tag)], [(object[1][1], object[1][2])])

# fluid = gmsh.model.occ.cut([(3, fluid.tag)], [(object[1][1], object[1][2])])
gmsh.model.occ.synchronize()
volume = ov[end][2]

# gmsh.model.occ.removeAllDuplicates()  # Clean up duplicates if needed
gmsh.model.occ.synchronize()  # Single synchronize after loop

# Identify each of the surfaces
surfaces = gmsh.model.occ.getEntities(2)
surface_clam_marker = []
surface_boundaries_marker = []

for surface in surfaces
    com = gmsh.model.occ.getCenterOfMass(surface[1], surface[2])
    if xmin < com[1] < xmax && ymin < com[2] < ymax && zmin < com[3] < zmax
        push!(surface_clam_marker, surface[2])
    elseif (abs(com[1]) - L/2 - 0.99 * d_PML) > 0 || (abs(com[3]) - w/2 - 0.99 * d_PML) > 0 || (com[2] + 0.99 * d_PML) < 0 || (com[2] - H - 0.99*d_PML)> 0
        push!(surface_boundaries_marker, surface[2])
    end  
end

# 3D groups for the resulting fragments
pml_groups = Dict(
    "fluid_pml_xx" => [fluid_right_pml_xx.tag, fluid_left_pml_xx.tag],
    "fluid_pml_yy" => [fluid_top_pml_yy.tag, fluid_bottom_pml_yy.tag],
    "fluid_pml_zz" => [fluid_front_pml_zz.tag, fluid_back_pml_zz.tag],
    
    "fluid_pml_yz" => [fluid_front_top_pml_yz.tag, fluid_back_top_pml_yz.tag, fluid_front_ground.tag, fluid_back_ground.tag],
    "fluid_pml_xz" => [fluid_right_front_pml_xz.tag, fluid_left_front_pml_xz.tag,fluid_right_back_pml_xz.tag, fluid_left_back_pml_xz.tag],
    "fluid_pml_xy" => [fluid_right_front_pml_xy.tag, fluid_left_front_pml_xy.tag, fluid_left_ground_pml_xy.tag, fluid_right_ground_pml_xy.tag],
    
    "fluid_pml_xyz" => [fluid_right_front_top_pml_xyz.tag, fluid_left_front_top_pml_xyz.tag, fluid_right_back_top_pml_xyz.tag,
                        fluid_left_back_top_pml_xyz.tag, fluid_right_front_ground_pml_xyz.tag, fluid_left_front_ground_pml_xyz.tag,fluid_right_back_ground_pml_xyz.tag, fluid_left_back_ground_pml_xyz.tag]
)


# # Set 3D-physical groups for the resulting fragments
object_tag = gmsh.model.addPhysicalGroup(3, [object[1][2]])
fluid_tag = gmsh.model.addPhysicalGroup(3, [volume])
pml_xx_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xx"])
pml_yy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yy"])
pml_zz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_zz"])
pml_xy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xy"])
pml_yz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yz"])
pml_xz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xz"])
pml_xyz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xyz"])

# Set 2D-Physical groups for the resulting fragments
fluid_boundaries_tag = gmsh.model.addPhysicalGroup(2, surface_boundaries_marker)
cube_little_tag = gmsh.model.addPhysicalGroup(2, surface_clam_marker)

# Set 3D-physical names for the physical groups
gmsh.model.setPhysicalName(3, fluid_tag, "Fluid")
gmsh.model.setPhysicalName(3, object_tag, "Object")
gmsh.model.setPhysicalName(3, pml_xx_tag, "PML_xx")
gmsh.model.setPhysicalName(3, pml_yy_tag, "PML_yy")
gmsh.model.setPhysicalName(3, pml_zz_tag, "PML_zz")
gmsh.model.setPhysicalName(3, pml_xy_tag, "PML_xy")
gmsh.model.setPhysicalName(3, pml_yz_tag, "PML_yz")
gmsh.model.setPhysicalName(3, pml_xz_tag, "PML_xz")
gmsh.model.setPhysicalName(3, pml_xyz_tag, "PML_xyz")

# Set 2D-physical names for the physical groups
gmsh.model.setPhysicalName(2, fluid_boundaries_tag, "Boundaries")
gmsh.model.setPhysicalName(2, cube_little_tag, "Clam")

# Synchronize the model
gmsh.model.occ.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h1)

# Generate a 3D mesh
gmsh.model.mesh.generate(3)


# Save the resulting mesh to a file
gmsh.write("./data/fluid_clam_brep_fragment.msh")

# Finalize Gmsh session
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/fluid_clam_brep_fragment.msh")
# Write the mesh to a vtk file

Γ = BoundaryTriangulation(model, tags="Object")
normals = get_normal_vector(Γ)

writevtk(model, "./results/fluid_clam_brep_fragment")
writevtk(Γ, "./results/fluid_clam_brep_fragment_boundary", cellfields = ["Normals"=>normals])

