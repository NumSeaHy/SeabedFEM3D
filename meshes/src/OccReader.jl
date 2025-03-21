using PyCall

occ = pyimport("OCC.Core")
# Import required OCC modules
TopExp = pyimport("OCC.Core.TopExp")
TopAbs = pyimport("OCC.Core.TopAbs")
topods = pyimport("OCC.Core.TopoDS")
BRep = pyimport("OCC.Core.BRep")
breptools = pyimport("OCC.Core.BRepTools")


builder = occ.BRep.BRep_Builder()
shape = occ.TopoDS.TopoDS_Shape()
filename = "./clam_geometries/EntireClam.brep"
breptools.breptools_Read(shape, filename, builder)

bnd = pyimport("OCC.Core.Bnd")
brepbndlib = pyimport("OCC.Core.BRepBndLib")
box = bnd.Bnd_Box()

# The third argument (True) ensures that the box is computed with the maximum bounding range.
brepbndlib.brepbndlib_Add(shape, box, true)

# Retrieve the bounding box limits.
xmin, ymin, zmin, xmax, ymax, zmax = box.Get()



# Assume 'shape' is already defined (e.g., loaded from a BREP file)
# Create an explorer for vertices in the shape
explorer = TopExp.TopExp_Explorer(shape, TopAbs.TopAbs_VERTEX)

# Initialize an empty Julia array to hold the points
points = []

# Iterate over all vertices in the shape
while explorer.More()
    # Convert the current item to a vertex
    vertex = topods.topods_Vertex(explorer.Current())
    # Get the point coordinates of the vertex
    pnt = BRep.BRep_Tool.Pnt(vertex)
    # Append the (x, y, z) tuple to the points array
    push!(points, [pnt.X(), pnt.Y(), pnt.Z()])
    # Move to the next vertex
    explorer.Next()
end

points_correct_type = [v for v in points]

using BoundingSphere

# Compute the bounding sphere of the points
center, radius = boundingsphere(points_correct_type)

# Import the necessary OCC modules
gp = pyimport("OCC.Core.gp")
BRepPrimAPI = pyimport("OCC.Core.BRepPrimAPI")

# Define the center of the sphere (here at the origin) and its radius.
center = gp.gp_Pnt(center[1], center[2], center[3])
radius = radius

# Create the sphere using BRepPrimAPI_MakeSphere
sphere_shape = BRepPrimAPI.BRepPrimAPI_MakeSphere(center, radius).Shape()

BRepTools = pyimport("OCC.Core.BRepTools")

# Define a filename and path. Here, we save it in the current directory.
filename = "./clam_geometries/evolving_sphere_entire_clam.brep"

# Write the BRep file.
BRepTools.breptools_Write(sphere_shape, filename)

BRepBuilderAPI = pyimport("OCC.Core.BRepBuilderAPI")

# Define the translation vector, e.g., move by dx, dy, dz
dx = 100.0
dy = 50.0
dz = -30.0

# Create a transformation object and set it to a translation
trsf = gp.gp_Trsf()
trsf.SetTranslation(gp.gp_Vec(dx, dy, dz))

# Apply the transformation to the shape
transformer = BRepBuilderAPI.BRepBuilderAPI_Transform(shape, trsf)
translated_shape = transformer.Shape()

# Define the output file name/path
filename = "clam_geometries/translated_shape.brep"


# Write the translated shape to file
breptools.breptools_Write(translated_shape, filename)

println("Translated BREP saved as: ", filename)

builder = occ.BRep.BRep_Builder()
shape = occ.TopoDS.TopoDS_Shape()
filename = "./clam_geometries/translated_shape.brep"
breptools.breptools_Read(shape, filename, builder)


# Assume 'shape' is already defined (e.g., loaded from a BREP file)
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

points_correct_type = [v for v in points]

# Compute the bounding sphere of the points
center, radius = boundingsphere(points_correct_type)

# Import the necessary OCC modules
gp = pyimport("OCC.Core.gp")
BRepPrimAPI = pyimport("OCC.Core.BRepPrimAPI")

# Define the center of the sphere (here at the origin) and its radius.
center = gp.gp_Pnt(center[1], center[2], center[3])
radius = radius

# Create the sphere using BRepPrimAPI_MakeSphere
sphere_shape = BRepPrimAPI.BRepPrimAPI_MakeSphere(center, radius).Shape()


# Define a filename and path. Here, we save it in the current directory.
filename = "./clam_geometries/evolving_sphere_entire_clam2_t.brep"

# Write the BRep file.
breptools.breptools_Write(sphere_shape, filename)


