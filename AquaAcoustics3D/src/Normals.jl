using Gmsh

gmsh.initialize()

# The API provides access to geometrical data in a CAD kernel agnostic manner.

# Let's create a simple CAD model by fusing a sphere and a cube, then mesh the
# surfaces:
gmsh.model.add("x5")
s = gmsh.model.occ.addSphere(0, 0, 0, 1)
b = gmsh.model.occ.addBox(0.5, 0, 0, 1.3, 2, 3)
gmsh.model.occ.fuse([(3, s)], [(3, b)])
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)

# We can for example retrieve the exact normals and the curvature at all the
# mesh nodes (i.e. not normals and curvatures computed from the mesh, but
# directly evaluated on the geometry), by querying the CAD kernels at the
# corresponding parametric coordinates.
normals = []

# For each surface in the model:
for e in gmsh.model.getEntities(2)
    # Retrieve the surface tag
    s = e[2]

    # Get the mesh nodes on the surface, including those on the boundary
    # (contrary to internal nodes, which store their parametric coordinates,
    # boundary nodes will be reparametrized on the surface in order to compute
    # their parametric coordinates, the result being different when
    # reparametrized on another adjacent surface)
    tags, coord, param = gmsh.model.mesh.getNodes(2, s, true)

    # Get the surface normals on all the points on the surface corresponding to
    # the parametric coordinates of the nodes
    norm = gmsh.model.getNormal(s, param)

    # Store the normals and the curvatures so that we can display them as
    # list-based post-processing views
    for i in 1:3:length(coord)
        push!(normals, coord[i])
        push!(normals, coord[i + 1])
        push!(normals, coord[i + 2])
        push!(normals, norm[i])
        push!(normals, norm[i + 1])
        push!(normals, norm[i + 2])
    end
end
# Create a list-based vector view on points to display the normals, and a scalar
# view on points to display the curvatures
vn = gmsh.view.add("normals")
gmsh.view.addListData(vn, "VP", length(normals) ÷ 6, normals)
gmsh.view.option.setNumber(vn, "ShowScale", 0)
gmsh.view.option.setNumber(vn, "ArrowSizeMax", 30)
gmsh.view.option.setNumber(vn, "ColormapNumber", 19)

gmsh.fltk.run()

gmsh.finalize()
