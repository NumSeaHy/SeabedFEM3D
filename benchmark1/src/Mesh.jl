using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

include("Configuration.jl") 

function FunMesh(L, H, w, r, d_PML, x₀, y₀, z₀, c, f, N, name)

    # Initialize Gmsh and add a model
    gmsh.initialize()
    gmsh.model.add("Rectangle_and_Sphere")

    # Number of elements in the mesh
    h1 = c/(N*f)
    
    # Use the OCC kernel (default for boolean operations)
    rect_tag = gmsh.model.occ.addBox(-L/2-d_PML, -d_PML, -w/2-d_PML, L+2*d_PML, H+2*d_PML, w+2*d_PML)  # Add a rectangular domain (box)
    sphere_tag = gmsh.model.occ.addSphere(0, H/2, 0, r)      # Add a sphere

    # Perform a boolean subtraction (subtract the sphere from the box)
    fluid = gmsh.model.occ.cut([(3, rect_tag)], [(3, sphere_tag)])
    volume = fluid[1]

    # Synchronize after the boolean operation
    gmsh.model.occ.synchronize()

    # Identify each of the surfaces
    surfaces = gmsh.model.occ.getEntities(2)
    surface_sphere_marker = []
    surface_boundaries_marker = []

    for surface in surfaces
        com = gmsh.model.occ.getCenterOfMass(surface[1], surface[2])
        if sqrt((com[1]-x₀)^2 + (com[2]-y₀)^2 + (com[3]-z₀)^2) < r
            push!(surface_sphere_marker, surface[2])
        else
            push!(surface_boundaries_marker, surface[2])
        end  
    end

    # Create a physical group for the fluid domain and assign the remaining volume
    f_volume = gmsh.model.addPhysicalGroup(3, [volume[1][2]])
    surface_sphere_marker = gmsh.model.addPhysicalGroup(2, surface_sphere_marker)
    surface_boundaries_marker = gmsh.model.addPhysicalGroup(2, surface_boundaries_marker)


    # Set the physical name no the physical group
    gmsh.model.setPhysicalName(3, f_volume, "Fluid Domain")
    gmsh.model.setPhysicalName(2, surface_sphere_marker, "Sphere Surface")
    gmsh.model.setPhysicalName(2, surface_boundaries_marker, "Boundaries")

    # Set the mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h1)

    # Synchronize to register the physical group
    gmsh.model.occ.synchronize()

    # Generate a 3D mesh
    gmsh.model.mesh.generate(3)

    # Save the resulting mesh to a file
    gmsh.write("./data/" * "$name" * ".msh")

    # Finalize Gmsh session
    gmsh.finalize()

    # Convert the mesh to the Gridap format
    model = GmshDiscreteModel("./data/" * "$name" * ".msh")
    # Write the mesh to a vtk file
    writevtk(model,"./results/"* "$name")
end

NumberElements = [6, 9, 12, 15]
Meshes = ["coarse", "medium", "fine", "extrafine"]

for (i, data) in enumerate(NumberElements)
    FunMesh(L, H, w, r, d_PML, x₀, y₀, z₀, c_F(ω), f, data, Meshes[i])
end