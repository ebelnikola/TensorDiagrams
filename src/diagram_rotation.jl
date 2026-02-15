#################################################
# TOC
# topic 0:           DICTIONARY ROTATING UTILITES
# topic 1:           ROTATION (NODE)
# topic 2:           ROTATION (DIAGRAM)
#################################################

#################################################
# topic 0: DICTIONARY ROTATING UTILITES
#################################################
"""
    rotate_boundary_vectors_ccw(d::Dict{String, Vector{Int}})

Rotates boundary vectors (legs or position indices) 90 degrees CCW.
Applies necessary reversals to maintain standard ordering conventions:
- Left (Bottom->Top) -> Bottom (Left->Right): Reversed
- Right (Bottom->Top) -> Top (Left->Right): Reversed
- Top (Left->Right) -> Left (Bottom->Top): Preserved
- Bottom (Left->Right) -> Right (Bottom->Top): Preserved
"""
function rotate_boundary_vectors_ccw(d::Dict{String,Vector{Int}})
    new_d = Dict{String,Vector{Int}}()
    if haskey(d, "top")
        new_d["left"] = d["top"]
    end
    if haskey(d, "left")
        new_d["bottom"] = reverse(d["left"])
    end
    if haskey(d, "bottom")
        new_d["right"] = d["bottom"]
    end
    if haskey(d, "right")
        new_d["top"] = reverse(d["right"])
    end
    return new_d
end

"""
    rotate_boundary_scalars_ccw(d::Dict{String, Int})

Rotates scalar boundary properties (slot numbers) 90 degrees CCW.
"""
function rotate_boundary_scalars_ccw(d::Dict{String,Int})
    new_d = Dict{String,Int}()
    if haskey(d, "top")
        new_d["left"] = d["top"]
    end
    if haskey(d, "left")
        new_d["bottom"] = d["left"]
    end
    if haskey(d, "bottom")
        new_d["right"] = d["bottom"]
    end
    if haskey(d, "right")
        new_d["top"] = d["right"]
    end
    return new_d
end

#########################################################
# topic 1: ROTATION (NODE)
#########################################################

"""
    rotate(node::TensorNode, k::Int)

Rotate the `TensorNode` `k` times by 90 degrees anti-clockwise.
Updates the `rotation` field and rotates the legs accordingly.
If `k` is negative, it corresponds to clockwise rotation.

# Arguments
- `node::TensorNode`: The node to rotate.
- `k::Int`: The number of 90-degree rotations.

# Returns
- A new rotated `TensorNode`.
"""
function rotate(node::TensorNode, k::Int=1)
    new_node = copy(node)

    # 1. Update rotation counter
    # Normalize k to be between 0 and 3
    k = mod(k, 4)
    new_node.rotation = mod(new_node.rotation + k, 4)

    # 2. Rotate legs
    # Mapping for 1 CCW rotation: Top->Left, Left->Bottom, Bottom->Right, Right->Top
    # legs on Left and Right sides needs to be reversed because of the drawing order
    # we also swap width and height
    # We apply this k times
    for _ in 1:k
        new_node.legs = rotate_boundary_vectors_ccw(new_node.legs)
        new_node.width, new_node.height = new_node.height, new_node.width
    end

    return new_node
end

function unrotate(node::TensorNode)
    rotate(node, -node.rotation)
end

#########################################################
# topic 2: ROTATION (DIAGRAM)
#########################################################

"""
    rotate(diag::TensorDiagram, k::Int=1)

Rotate the entire `TensorDiagram` `k` times by 90 degrees anti-clockwise.
If `k` is negative, performs clockwise rotation.

# Arguments
- `diag::TensorDiagram`: The diagram to rotate.
- `k::Int=1`: Number of rotations.

# Returns
- A new rotated `TensorDiagram`.
"""
function rotate(diag::TensorDiagram, k::Int=1)
    diag_r = copy(diag)

    k = mod(k, 4)
    if k == 0
        return diag_r
    end

    # 1. Rotate nodes
    diag_r.nodes = rotate.(diag.nodes, k)

    # 2. Permute Boundary Properties and Coordinates

    for _ in 1:k
        diag_r.boundary_slots_num = rotate_boundary_scalars_ccw(diag_r.boundary_slots_num)
        diag_r.boundary_legs = rotate_boundary_vectors_ccw(diag_r.boundary_legs)
        diag_r.boundary_legs_posidx = rotate_boundary_vectors_ccw(diag_r.boundary_legs_posidx)

        # 3. Swap width/height
        diag_r.width, diag_r.height = diag_r.height, diag_r.width

        # 4. Rotate Coordinates (1 step CCW: (x, y) -> (-y, x))
        if !isnothing(diag_r.node_coordinates)
            diag_r.node_coordinates = map(diag_r.node_coordinates) do coord
                if ismissing(coord)
                    missing
                else
                    Vector{Float64}([-coord[2], coord[1]])
                end
            end
        end
    end

    return diag_r
end
