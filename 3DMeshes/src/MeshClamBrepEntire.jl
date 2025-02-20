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
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 2) 
gmsh.model.add("BuriedClams")


# Mesh size parameter (ensure variables c and f are defined in Configuration.jl)
h_f = c_F(ω) / (15 * f)
h_p = c_P(ω) / (15 * f)

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

gmsh.model.occ.synchronize()  # Single synchronize after loop

# Open the .brep geometry file
object = gmsh.model.occ.importShapes("./clam_geometries/EntireClam.brep")
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
dy = 1/2*(t_P - ymax - ymin)
dz = (-zmax-zmin)/2
gmsh.model.occ.translate([(3, object[1][2])], dx, 0, 0)
gmsh.model.occ.translate([(3, object[1][2])], 0, dy, 0)
gmsh.model.occ.translate([(3, object[1][2])], 0, 0, dz)

# Compute the new bounding box
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(object[1][1], object[1][2])

# Perform a boolean subtraction (subtract the sphere from the box)
cut_operation = gmsh.model.occ.cut([(3, porous.tag)], [(object[1][1], object[1][2])])

gmsh.model.occ.synchronize()

porous.tag = cut_operation[1][1][2]

gmsh.model.occ.removeAllDuplicates()  # Clean up duplicates if needed
gmsh.model.occ.synchronize()  # Single synchronize after loop

# Identify each of the surfaces
surfaces = gmsh.model.occ.getEntities(2)
surface_clam_marker = []
surface_boundaries_marker = []

for surface in surfaces
    com = gmsh.model.occ.getCenterOfMass(surface[1], surface[2])
    if xmin < com[1] < xmax && ymin < com[2] < ymax && zmin < com[3] < zmax
        push!(surface_clam_marker, surface[2])
    elseif (abs(com[1]) - L/2 - 0.99 * d_PML) > 0 || (abs(com[3]) - w/2 - 0.99 * d_PML) > 0 || (com[2] + 0.99 * d_PML) < 0 || (com[2] - t_P - t_F - 0.99 * d_PML) > 0
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


# Set 3D-physical groups for the fluid domain together with the fluid PMLs
fluid_domain = gmsh.model.addPhysicalGroup(3, [fluid.tag])
pml_xx_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xx"])
pml_yy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yy"])
pml_zz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_zz"])
pml_xy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xy"])
pml_yz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_yz"])
pml_xz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xz"])
pml_xyz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["fluid_pml_xyz"])

# Set 3D-physical groups for the porous domain together with the porous PMLs
porous_tag = gmsh.model.addPhysicalGroup(3, [porous.tag])
porous_pml_xx_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xx"])
porous_pml_yy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_yy"])
porous_pml_zz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_zz"])
porous_pml_xy_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xy"])
porous_pml_yz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_yz"])
porous_pml_xz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xz"])
porous_pml_xyz_tag = gmsh.model.addPhysicalGroup(3, pml_groups["porous_pml_xyz"])

# Set 2D-Physical groups for the resulting fragments
fluid_boundaries_tag = gmsh.model.addPhysicalGroup(2, surface_boundaries_marker)
clam_objects_surface = gmsh.model.addPhysicalGroup(2, surface_clam_marker)

# Set 3D-physical names for the fluid physical groups
gmsh.model.setPhysicalName(3, fluid_domain, "Fluid")
gmsh.model.setPhysicalName(3, pml_xx_tag, "fluid_PML_xx")
gmsh.model.setPhysicalName(3, pml_yy_tag, "fluid_PML_yy")
gmsh.model.setPhysicalName(3, pml_zz_tag, "fluid_PML_zz")
gmsh.model.setPhysicalName(3, pml_xy_tag, "fluid_PML_xy")
gmsh.model.setPhysicalName(3, pml_yz_tag, "fluid_PML_yz")
gmsh.model.setPhysicalName(3, pml_xz_tag, "fluid_PML_xz")
gmsh.model.setPhysicalName(3, pml_xyz_tag, "fluid_PML_xyz")

# Set 3D-physical names for the porous physical groups
gmsh.model.setPhysicalName(3, porous_tag, "Porous")
gmsh.model.setPhysicalName(3, porous_pml_xx_tag, "porous_PML_xx")
gmsh.model.setPhysicalName(3, porous_pml_yy_tag, "porous_PML_yy")
gmsh.model.setPhysicalName(3, porous_pml_zz_tag, "porous_PML_zz")
gmsh.model.setPhysicalName(3, porous_pml_xy_tag, "porous_PML_xy")
gmsh.model.setPhysicalName(3, porous_pml_yz_tag, "porous_PML_yz")
gmsh.model.setPhysicalName(3, porous_pml_xz_tag, "porous_PML_xz")
gmsh.model.setPhysicalName(3, porous_pml_xyz_tag, "porous_PML_xyz")

# Set 2D-physical names for the physical groups
gmsh.model.setPhysicalName(2, fluid_boundaries_tag, "Boundaries")
gmsh.model.setPhysicalName(2, clam_objects_surface, "Clam")

# Synchronize the model
gmsh.model.occ.synchronize()

# Get the points that constitutes each of the domains
ov_f = gmsh.model.getBoundary([(3, fluid_domain.tag) for fluid_domain in all_fluid_domains], false, false, true)
ov_p = gmsh.model.getBoundary([(3, porous_domain.tag) for porous_domain in all_porous_domains], false, false, true)

# Set a specific mesh size to each of that points
gmsh.model.mesh.setSize(ov_f, h_f)
gmsh.model.mesh.setSize(ov_p, h_p)

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h1)

# Generate a 3D mesh
gmsh.model.mesh.generate(3)

# Save the resulting mesh to a file
gmsh.write("./data/fluid_clam_brep.msh")

# Finalize Gmsh session
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/fluid_clam_brep.msh")

# Write the mesh to a vtk file
writevtk(model,"./results/fluid_clam_brep")

# Extract the boundary triangulation
Γ = BoundaryTriangulation(model, tags="Clam")

# Extract the normal vectors
normals = get_normal_vector(Γ)

# Write the boundary triangulation to a vtk file
writevtk(Γ, "./results/fluid_clam_brep_boundary", cellfields = ["Normals"=>normals])


