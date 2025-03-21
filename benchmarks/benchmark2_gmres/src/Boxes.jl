module Boxes

export Box, get_box_vertices, get_box_surfaces, boxes_share_surface, create_pml_boxes, add_box_to_gmsh

using Gmsh  # Import Gmsh to use `gmsh.model.occ.addBox`

mutable struct Box
    x₀::Float64  # x-coordinate of the origin
    y₀::Float64  # y-coordinate of the origin
    z₀::Float64  # z-coordinate of the origin
    DX::Float64  # Increment in the x-coordinate with respect to x₀
    DY::Float64  # Increment in the y-coordinate with respect to y₀
    DZ::Float64  # Increment in the z-coordinate with respect to z₀
    tag::Int     # Tag of the box in Gmsh
    name::String  # Name of the box
end

# Function to initialize and tag a Box in Gmsh
function add_box_to_gmsh(x₀, y₀, z₀, DX, DY, DZ, name)
    box = Box(x₀, y₀, z₀, DX, DY, DZ, 0, name)
    box.tag = gmsh.model.occ.addBox(x₀, y₀, z₀, DX, DY, DZ)
    return box
end

# Function to create PML boxes with specified parameters
function create_pml_boxes(L, H, w, d_PML)
    # Define PML boxes by type and return them as a dictionary for easier reference
    pml_boxes = Dict(
        "right_pml_xx" => add_box_to_gmsh(L/2, 0, -w/2, d_PML, H, w, "right_pml_xx"),
        "left_pml_xx" => add_box_to_gmsh(-L/2 - d_PML, 0, -w/2, d_PML, H, w, "left_pml_xx"),
        
        "ground_pml_yy" => add_box_to_gmsh(-L/2, -d_PML, -w/2, L, d_PML, w, "ground_pml_yy"),
        "top_pml_yy" => add_box_to_gmsh(-L/2, H, -w/2, L, d_PML, w, "top_pml_yy"),
        
        "front_pml_zz" => add_box_to_gmsh(-L/2, 0, w/2, L, H, d_PML, "front_pml_zz"),
        "back_pml_zz" => add_box_to_gmsh(-L/2, 0, -w/2 - d_PML, L, H, d_PML, "back_pml_zz"),
        
        # Modified entries
        "right_front_pml_xz" => add_box_to_gmsh(L/2, 0, w/2, d_PML, H, d_PML, "right_front_pml_xz"),
        "left_front_pml_xz" => add_box_to_gmsh(-L/2 - d_PML, 0, w/2, d_PML, H, d_PML, "left_front_pml_xz"),
        "right_back_pml_xz" => add_box_to_gmsh(L/2, 0, -w/2 - d_PML, d_PML, H, d_PML, "right_back_pml_xz"),
        "left_back_pml_xz" => add_box_to_gmsh(-L/2 - d_PML, 0, -w/2 - d_PML, d_PML, H, d_PML, "left_back_pml_xz"),
        
        "front_top_pml_yz" => add_box_to_gmsh(-L/2, H, w/2, L, d_PML, d_PML, "front_top_pml_yz"),
        "front_bottom_pml_yz" => add_box_to_gmsh(-L/2, -d_PML, w/2, L, d_PML, d_PML, "front_bottom_pml_yz"),
        "back_top_pml_yz" => add_box_to_gmsh(-L/2, H, -w/2 - d_PML, L, d_PML, d_PML, "back_top_pml_yz"),
        "back_bottom_pml_yz" => add_box_to_gmsh(-L/2, -d_PML, -w/2 - d_PML, L, d_PML, d_PML, "back_bottom_pml_yz"),
        
        "right_front_pml_xy" => add_box_to_gmsh(L/2, H, w/2, d_PML, d_PML, -w, "right_front_pml_xy"),
        "left_front_pml_xy" => add_box_to_gmsh(-L/2 - d_PML, H, w/2, d_PML, d_PML, -w, "left_front_pml_xy"),
        "right_back_pml_xy" => add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, -w, "right_back_pml_xy"),
        "left_back_pml_xy" => add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, -w, "left_back_pml_xy"),
        
        "right_front_top_pml_xyz" => add_box_to_gmsh(L/2, H, w/2, d_PML, d_PML, d_PML, "right_front_top_PML_xyz"),
        "left_front_top_pml_xyz" => add_box_to_gmsh(-L/2 - d_PML, H, w/2, d_PML, d_PML, d_PML, "left_front_top_PML_xyz"),
        "right_front_bottom_pml_xyz" => add_box_to_gmsh(L/2, -d_PML, w/2, d_PML, d_PML, d_PML, "right_front_bottom_PML_xyz"),
        "left_front_bottom_pml_xyz" => add_box_to_gmsh(-L/2 - d_PML, -d_PML, w/2, d_PML, d_PML, d_PML, "left_front_bottom_PML_xyz"),
        "right_back_top_pml_xyz" => add_box_to_gmsh(L/2, H, -w/2 - d_PML, d_PML, d_PML, d_PML, "right_back_top_PML_xyz"),
        "left_back_top_pml_xyz" => add_box_to_gmsh(-L/2 - d_PML, H, -w/2 - d_PML, d_PML, d_PML, d_PML, "left_back_top_PML_xyz"),
        "right_back_bottom_pml_xyz" => add_box_to_gmsh(L/2, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "right_back_bottom_PML_xyz"),
        "left_back_bottom_pml_xyz" => add_box_to_gmsh(-L/2 - d_PML, -d_PML, -w/2 - d_PML, d_PML, d_PML, d_PML, "left_back_bottom_PML_xyz")
    )
    return pml_boxes
end

# Functions to calculate vertices and surfaces
function get_box_vertices(box::Box)
    x0, y0, z0 = box.x₀, box.y₀, box.z₀
    DX, DY, DZ = box.DX, box.DY, box.DZ
    return [
        (x0, y0, z0), (x0 + DX, y0, z0), (x0 + DX, y0 + DY, z0), (x0, y0 + DY, z0),
        (x0, y0, z0 + DZ), (x0 + DX, y0, z0 + DZ), (x0 + DX, y0 + DY, z0 + DZ), (x0, y0 + DY, z0 + DZ)
    ]
end

function get_box_surfaces(vertices::Vector{Tuple{Float64, Float64, Float64}})
    return [
        Set([vertices[1], vertices[2], vertices[3], vertices[4]]),  # Bottom face
        Set([vertices[5], vertices[6], vertices[7], vertices[8]]),  # Top face
        Set([vertices[1], vertices[2], vertices[6], vertices[5]]),  # front face
        Set([vertices[3], vertices[4], vertices[8], vertices[7]]),  # Back face
        Set([vertices[1], vertices[4], vertices[8], vertices[5]]),  # Left face
        Set([vertices[2], vertices[3], vertices[7], vertices[6]])   # Right face
    ]
end

# Function to check if two boxes share a surface
function boxes_share_surface(box1::Box, box2::Box)
    vertices1 = get_box_vertices(box1)
    vertices2 = get_box_vertices(box2)
    surfaces1 = get_box_surfaces(vertices1)
    surfaces2 = get_box_surfaces(vertices2)
    
    for surface1 in surfaces1
        for surface2 in surfaces2
            if surface1 == surface2
                return true
            end
        end
    end
    return false
end

end  # End of module Boxes
