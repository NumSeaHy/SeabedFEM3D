"""
This file generates the mesh of the PML domain using Gmsh 
The mesh includes the fluid domain, the porous domain and the PML and it is centered at the origin.
The mesh is saved in the folder data in the format .msh and .json
"""
# Load the packages
using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh
using LinearAlgebra
using SparseArrays
using LinearAlgebra

function compute_tangent_normal(points::Vector{Vector{Float64}})
    n = length(points)
    tangents = Vector{Vector{Float64}}(undef, n)
    normals = Vector{Vector{Float64}}(undef, n)

    # Define a reference direction for normals (e.g., pointing upwards)
    reference_normal = [0.0, 1.0]

    for i in 2:n-1
        # Compute the forward and backward differences
        forward_diff = [points[i+1][1] - points[i][1], points[i+1][2] - points[i][2]]
        backward_diff = [points[i][1] - points[i-1][1], points[i][2] - points[i-1][2]]

        # Average the differences to approximate the tangent vector
        tangent = normalize([(forward_diff[1] + backward_diff[1]) / 2,
                             (forward_diff[2] + backward_diff[2]) / 2])

        # Rotate the tangent vector by 90 degrees counterclockwise to get the normal vector
        normal = [-tangent[2], tangent[1]]

        # Ensure the normal vector points in the same direction as the reference normal
        if dot(normal, reference_normal) < 0
            normal = [-normal[1], -normal[2]]
        end

        # Store the unit tangent and normal vectors
        tangents[i] = tangent
        normals[i] = normal
    end

    # Handle the first point using the forward difference
    tangents[1] = normalize([points[2][1] - points[1][1], points[2][2] - points[1][2]])
    normals[1] = [-tangents[1][2], tangents[1][1]]
    if dot(normals[1], reference_normal) < 0
        normals[1] = [-normals[1][1], -normals[1][2]]
    end

    # Handle the last point using the backward difference
    tangents[n] = normalize([points[n][1] - points[n-1][1], points[n][2] - points[n-1][2]])
    normals[n] = [-tangents[n][2], tangents[n][1]]
    if dot(normals[n], reference_normal) < 0
        normals[n] = [-normals[n][1], -normals[n][2]]
    end

    return normals
end
# Initialize gmsh
gmsh.initialize()
gmsh.model.add("fluid_and_porous_PML")

# Mesh size depending on the angular frequency
h = 0.01
h_clam = h/6

# Input of the physical parameters neccesary to calculate the mesh size
L = 1 # Horizontal length of the domain [m]
t_P = 0.5 # Thickness of the porous domain [m]
# Clam parametrization
a = 0.03 # Major semi-axis of the ellipse [m]
b = a/2 # Minor semi-axis of the ellipse [m]
e = 0.002 # Maximum espesor of the ellipse [-]
θ₀ = pi/15 # Minimum angle of the ellipse [rad]
x = 0
y = t_P/2
α = 0 # Angle of rotation of the ellipse [rad]
R(α) = [cos(α) -sin(α); sin(α) cos(α)] # Rigid-body rotation matrix

# x_c = 0
# y_c = 0.25
# r_in = 0.05  # Interior radius
# r_out = 0.1 # Exterior radius
# θ_b = 180 * pi/180  # Bisector angle of the opening in radians
# θ_o = 10 * pi/180 # Opening angle 

# Define points based on animal attributes
x1 = [a*cos(θ₀/2), b*sin(θ₀/2)]
x2 = [0, b]
x3 = [-a, 0]
x4 = [0, -b]
x5 = [a*cos(θ₀/2), -b*sin(θ₀/2)]
x6 = [(a+e)*cos(θ₀/2), (b+e)*sin(θ₀/2)]
x7 = [0, b+e]
x8 = [-(a+e), 0]
x9 = [0, -(b+e)]
x10 = [(a+e)*cos(θ₀/2), -(b+e)*sin(θ₀/2)]    

# Transform points with rotation and translation
aux1 = R(α) * x1 .+ [x, y]
aux2 = R(α) * x2 .+ [x, y]
aux3 = R(α) * x3 .+ [x, y]
aux4 = R(α) * x4 .+ [x, y]
aux5 = R(α) * x5 .+ [x, y]
aux6 = R(α) * x6 .+ [x, y]
aux7 = R(α) * x7 .+ [x, y]
aux8 = R(α) * x8 .+ [x, y]
aux9 = R(α) * x9 .+ [x, y]
aux10 = R(α) * x10 .+ [x, y]
    


# Corners of Porous Domain
p_P1 = gmsh.model.geo.addPoint(-L/2, 0, 0, h)
p_P2 = gmsh.model.geo.addPoint(L/2, 0, 0, h)
p_P3 = gmsh.model.geo.addPoint(L/2, t_P, 0, h) 
p_P4 = gmsh.model.geo.addPoint(-L/2, t_P, 0, h)


# Elliptical Clam
cp = gmsh.model.geo.addPoint(x, y, 0, h)
e_p1 = gmsh.model.geo.addPoint(aux1[1], aux1[2], 0, h_clam)
e_p2 = gmsh.model.geo.addPoint(aux2[1], aux2[2], 0, h_clam)
e_p3 = gmsh.model.geo.addPoint(aux3[1], aux3[2], 0, h_clam)
e_p4 = gmsh.model.geo.addPoint(aux4[1], aux4[2], 0, h_clam)
e_p5 = gmsh.model.geo.addPoint(aux5[1], aux5[2], 0, h_clam)
e_p6 = gmsh.model.geo.addPoint(aux6[1], aux6[2], 0, h_clam)
e_p7 = gmsh.model.geo.addPoint(aux7[1], aux7[2], 0, h_clam)
e_p8 = gmsh.model.geo.addPoint(aux8[1], aux8[2], 0, h_clam)
e_p9 = gmsh.model.geo.addPoint(aux9[1], aux9[2], 0, h_clam)
e_p10 = gmsh.model.geo.addPoint(aux10[1], aux10[2], 0, h_clam)

# Edges of the Porous
l_P1 = gmsh.model.geo.add_line(p_P1, p_P2)
l_P2 = gmsh.model.geo.add_line(p_P2, p_P3)
l_P3 = gmsh.model.geo.add_line(p_P3, p_P4)
l_P4 = gmsh.model.geo.add_line(p_P4, p_P1)


# Edges of the clam
el1 = gmsh.model.geo.addEllipseArc(e_p1, cp, e_p3, e_p2)
el2 = gmsh.model.geo.addEllipseArc(e_p2, cp, e_p3, e_p3)
el3 = gmsh.model.geo.addEllipseArc(e_p3, cp, e_p3, e_p4)
el4 = gmsh.model.geo.addEllipseArc(e_p4, cp, e_p3, e_p5)
el5 = gmsh.model.geo.addEllipseArc(e_p6, cp, e_p8, e_p7)
el6 = gmsh.model.geo.addEllipseArc(e_p7, cp, e_p8, e_p8)
el7 = gmsh.model.geo.addEllipseArc(e_p8, cp, e_p8, e_p9)
el8 = gmsh.model.geo.addEllipseArc(e_p9, cp, e_p8, e_p10)
lel1 = gmsh.model.geo.addLine(e_p6, e_p1)
lel2 = gmsh.model.geo.addLine(e_p5, e_p10)


# Curve Loops
cl_P = gmsh.model.geo.add_curve_loop([l_P1, l_P2, l_P3, l_P4])
cl_rock = gmsh.model.geo.addCurveLoop([el1, el2, el3, el4, lel2, -el8, -el7, -el6, -el5,  lel1])

# Surfaces
s_P = gmsh.model.geo.addPlaneSurface([cl_P, cl_rock])
s_clam = gmsh.model.geo.addPlaneSurface([cl_rock])

# Synchronize model before meshing
gmsh.model.geo.synchronize()

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Set physical groups for 1D entities
l_rock = gmsh.model.addPhysicalGroup(1, [el1, el2, el3, el4, lel2, el8, el7, el6, el5,  lel1])

# Set physical groups for 2D entities
f_Ph = gmsh.model.addPhysicalGroup(2, [s_P])
f_clam = gmsh.model.addPhysicalGroup(2, [s_clam])
# Set physical names for 1D entities

gmsh.model.setPhysicalName(1, l_rock, "rock")

# Set physical names for 2D entities
gmsh.model.setPhysicalName(2, f_Ph, "physical_domain")
gmsh.model.setPhysicalName(2, f_clam, "clam")

# Get the domain points that are inside the bounding box of the object
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, f_Ph)
all_nodes_tags, all_nodes, all_params = gmsh.model.mesh.getNodes(2, -1, false, true)
nodes_in_bbox_tags = []

for (j,i) in enumerate(1:3:length(all_nodes))
    x, y, z = all_nodes[i], all_nodes[i+1], all_nodes[i+2]
    if xmin <= x <= xmax && ymin <= y <= ymax && zmin <= z <= zmax
        push!(nodes_in_bbox_tags, all_nodes_tags[j])
    end
end

nodes_associated_tags = [gmsh.model.mesh.getNode(i)[1] for i in nodes_in_bbox_tags]

elements_by_coordinates_inside_bbox = gmsh.model.mesh.getElementsByCoordinates(nodes_associated_tags[1][1], nodes_associated_tags[1][2], nodes_associated_tags[1][3])

decimal_nodes_inbox_tag = [Int64(hex) for hex in nodes_in_bbox_tags]

f_element = gmsh.model.addPhysicalGroup(2, elements_by_coordinates_inside_bbox)
gmsh.model.geo.addPlaneSurface(elements_by_coordinates_inside_bbox)
gmsh.model.setPhysicalName(2, f_element, "toitriangle")

# lines_rock = gmshnodes_in_bbox_tags.model.getEntitiesForPhysicalGroup(1, l_rock)
# lines_rock_1 = [lines_rock[1], lines_rock[3], lines_rock[5], lines_rock[7]]
# # Get the mesh nodes on the surface, including those on the bounda


# total_coord = []
# for i in eachindex(lines_rock_1)
#     tags, coord, param = gmsh.model.mesh.getNodes(1, lines_rock[i], true)
#     for j in 1:3:length(coord)
#         push!(total_coord, coord[j])
#         push!(total_coord, coord[j+1])
#         push!(total_coord, coord[j+2])
#     end
# end

# # Extracting the x and y compoenents of the nodes in a vector
# points = [[total_coord[i], total_coord[i+1]] for i in 1:3:length(total_coord)]
# # Get the surface normals on all the points on the surface corresponding to
# # the parametric coordinates of the nodes
# # norm = gmsh.model.getNormal(lines_rock[1], param)
# normals = compute_tangent_normal(points)

# coords_normals = []
# # For each surface in the model:
# for i in eachindex(points)
#         push!(coords_normals, points[i][1])
#         push!(coords_normals, points[i][2])
#         push!(coords_normals, 0.0)
#         push!(coords_normals, normals[i][1])
#         push!(coords_normals, normals[i][2])
#         push!(coords_normals, 0.0)
# end
# # Create a list-based vector view on points to display the normals, and a scalar
# # view on points to display the curvatures
# vn = gmsh.view.add("normals")
# gmsh.view.addListData(vn, "VP", length(coords_normals) // 6, coords_normals)
# gmsh.view.option.setNumber(vn, "ShowScale", 0)
# gmsh.view.option.setNumber(vn, "ArrowSizeMax", 30)
# gmsh.view.option.setNumber(vn, "ColormapNumber", 19)


# gmsh.fltk.run()



# Write mesh to file
gmsh.write("./erase/mesh.msh")
gmsh.finalize()

# Convert the mesh to the Gridap format
model = GmshDiscreteModel("./erase/mesh.msh")

# Write the mesh to a vtk file
writevtk(model,"./erase/mesh")

Γ = BoundaryTriangulation(model, tags="rock")
n_nodes = get_physical_coordinate(Γ)
n_vectors = get_normal_vector(Γ)
cell_p = get_cell_points(Γ)
n_vectors(cell_p)
# Function to compute unit tangent and normal vectors for a 1D curve in 2D space

