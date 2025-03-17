using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh
using Distributions
using Revise

include("Configuration.jl")
includet("Boxes.jl")
includet("MarineForms.jl")
using .Boxes
using .MarineForm

# Define the meshing alforithms available in Gmsh to then choose
meshing_algorithms = Dict(
    "Delaunay"           => 1, # Best one, but should be used together with the option "Mesh.OptimizeNetgen = 1", if not Gridap will throw an erro
    "Frontal"            => 4, # Good local mesh around the cockles, but the global mesh is not good (some very very small elements det ≈ 0 ?)
    "Frontal Delaunay"   => 5, # Deprecated seems by looking at the Gmsh.pdf documentation
    "Frontal Hex"        => 6, # I think that not make sense use a hexahedral mesh for this problem... and also seems to be deprecated
    "MMG3D"              => 7,
    "R-tree"             => 9, # Almost impossible to use with so complicated geometry... it takes a lot of time to generate the tree
    "HXT"                => 10, # Similar to Delaunay but with a different algorithm that use redundant elements from my point of view
)

algorithm = "Delaunay"

# Initialize Gmsh and add a model
gmsh.initialize()
# gmsh.option.setString("General.LogFile", "gmsh.log")
gmsh.option.setNumber("Mesh.Algorithm3D", meshing_algorithms[algorithm])
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 6)
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
gmsh.option.setNumber("Geometry.Tolerance", 1e-9)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
# gmsh.option.setNumber("Geometry.OCCThruSectionsDegree", 1)
gmsh.model.add("BuriedClams")


# Mesh size parameter (ensure variables c and f are defined in Configuration.jl)
h_f = real(c_F(ω)) / (15 * f)
h_p = real(c_P(ω)) / (15 * f)
h_cockle = h_p/2

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

println("...Generating the scenarios...")
animals = generate_animals(definition, max_iterations)

rigid_cockle = [animal for animal in animals if isa(animal, RigidCockle)]

println("...Generating the geometries...")
v_cockle = [generate_geometry(animal) for animal in rigid_cockle]

# Perform a boolean subtraction (subtract the sphere from the box)
cut_operation = gmsh.model.occ.cut([(3, porous.tag)], [(3, i) for i in v_cockle])

gmsh.model.occ.synchronize()  # Single synchronize after loop

# After the cut operation, the porous domain is the first element of the resulting list since it is the only volume that remains (cockle vanishes!)
porous.tag = cut_operation[1][1][2]

# Identify each of the surfaces
surfaces = gmsh.model.occ.getEntities(2)
surface_cockles_marker = []
surface_boundaries_marker = []


for surface in surfaces
    com = gmsh.model.occ.getCenterOfMass(surface[1], surface[2])

    # Identify the surfaces where the cockles are located to apply the boundary condition of the incident field arriving
    for cockle in rigid_cockle
        xmin, ymin, zmin, xmax, ymax, zmax = cockle.bounding_box
        if xmin < com[1] < xmax && ymin < com[2] < ymax && zmin < com[3] < zmax
            push!(surface_cockles_marker, surface[2])
        end
    end
    
    # Identify the boundaries of the PML domain where the homogeneuous dirichlet boundary condition will be applied
    if (abs(com[1]) - L/2 - 0.99 * d_PML) > 0 || (abs(com[3]) - w/2 - 0.99 * d_PML) > 0 || (com[2] + 0.99 * d_PML) < 0 || (com[2] - t_P - t_F - 0.99 * d_PML) > 0
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



# Set the 2D-physical groups for the surfaces
surface_cockles_marker = gmsh.model.addPhysicalGroup(2, surface_cockles_marker)
surface_boundaries_marker = gmsh.model.addPhysicalGroup(2, surface_boundaries_marker)

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
gmsh.model.setPhysicalName(2, surface_boundaries_marker, "Boundaries")
gmsh.model.setPhysicalName(2, surface_cockles_marker, "Cockles")

# Synchronize the model
gmsh.model.occ.synchronize()

# Get the points that constitutes each of the domains
ov_f = gmsh.model.getBoundary([(3, fluid_domain.tag) for fluid_domain in all_fluid_domains], false, false, true)
ov_p = gmsh.model.getBoundary([(3, porous_domain.tag) for porous_domain in all_porous_domains], false, false, true)
ov_cockles = gmsh.model.getBoundary([(2, surface_cockles_marker[i]) for i in eachindex(surface_cockles_marker)], false, false, true)

# Set a specific mesh size to each of that points
gmsh.model.mesh.setSize(ov_f, h_f)
gmsh.model.mesh.setSize(ov_p, h_p)
gmsh.model.mesh.setSize(ov_cockles, h_cockle)

# Generate a 3D mesh
gmsh.model.mesh.generate(3)

# Save the resulting mesh to a file
gmsh.write("./data/fluid_clam_brep_"*algorithm*".msh")

# Finalize Gmsh session
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./data/fluid_clam_brep_"*algorithm*".msh")

# Write the mesh to a vtk file
writevtk(model,"./results/fluid_clam_brep_"*algorithm)

# # Extract the boundary triangulation
# Γ = BoundaryTriangulation(model, tags="Clam")

# # Extract the normal vectors
# normals = get_normal_vector(Γ)

# # Write the boundary triangulation to a vtk file
# writevtk(Γ, "./results/fluid_clam_brep_boundary", cellfields = ["Normals"=>normals])


