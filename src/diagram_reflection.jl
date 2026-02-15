#################################################
# TOC
# topic 1:           REFLECTION (NODE)
# topic 2:           REFLECTION (DIAGRAM)
#################################################


#########################################################
# topic 1: REFLECTION (NODE)
#########################################################
"""
    hor_reflection(node::TensorNode)

Create a horizontally reflected version of the TensorNode.
- Swaps left and right legs
- Reverses order of top and bottom legs
- Toggles the `hor_ref` flag

# Arguments:
- node: The TensorNode to reflect

# Returns:
- A new TensorNode with horizontal reflection applied
"""
function hor_reflection(node::TensorNode)
    new_node = copy(node)

    new_legs = Dict{String,Array{Int,1}}()
    # Swap left and right
    if haskey(node.legs, "left")
        new_legs["right"] = node.legs["left"]
    end
    if haskey(node.legs, "right")
        new_legs["left"] = node.legs["right"]
    end
    # Reverse top and bottom
    if haskey(node.legs, "top")
        new_legs["top"] = reverse(node.legs["top"])
    end
    if haskey(node.legs, "bottom")
        new_legs["bottom"] = reverse(node.legs["bottom"])
    end

    new_node.legs = new_legs
    new_node.hor_ref = !node.hor_ref

    return new_node
end
"""
    ver_reflection(node::TensorNode)

Create a vertically reflected version of the TensorNode.
- Swaps top and bottom legs
- Reverses order of left and right legs
- Toggles the `ver_ref` flag

# Arguments:
- node: The TensorNode to reflect

# Returns:
- A new TensorNode with vertical reflection applied
"""
function ver_reflection(node::TensorNode)
    new_node = copy(node)

    new_legs = Dict{String,Array{Int,1}}()
    # Swap top and bottom
    if haskey(node.legs, "top")
        new_legs["bottom"] = node.legs["top"]
    end
    if haskey(node.legs, "bottom")
        new_legs["top"] = node.legs["bottom"]
    end

    # Reverse left and right
    if haskey(node.legs, "left")
        new_legs["left"] = reverse(node.legs["left"])
    end
    if haskey(node.legs, "right")
        new_legs["right"] = reverse(node.legs["right"])
    end

    new_node.legs = new_legs
    new_node.ver_ref = !node.ver_ref

    return new_node
end
"""
    reflect(node::TensorNode; dir="vertical")

Create a reflected version of the TensorNode in the specified direction.

# Arguments
- `node::TensorNode`: The node to reflect.
- `dir::String="vertical"`: The direction of reflection. Can be "vertical" or "horizontal".

# Returns
- A new `TensorNode` with the reflection applied.

# Throws
- `ArgumentError`: If `dir` is not "vertical" or "horizontal".
"""
function reflect(node::TensorNode; dir="vertical")
    if dir == "vertical"
        return ver_reflection(node)
    elseif dir == "horizontal"
        return hor_reflection(node)
    else
        throw(ArgumentError("Invalid reflection direction: $dir. Use 'vertical' or 'horizontal'."))
    end
end


"""
    unreflect(node)

Restores the original orientation of a `TensorNode` by undoing any horizontal or vertical reflections stored in its metadata.

# Arguments
- `node`: The `TensorNode` (or similar object) to unreflect. Must support `hor_reflection` and `ver_reflection`, and have `hor_ref` and `ver_ref` boolean properties.

# Returns
- A new `TensorNode` returned to its base orientation (i.e., with `hor_ref` and `ver_ref` effectively applied to cancel out the stored state).
"""
function unreflect(node)
    if node.hor_ref
        node = hor_reflection(node)
    end
    if node.ver_ref
        node = ver_reflection(node)
    end
    return node
end


#########################################################
# topic 2: REFLECTION (DIAGRAM) 
#########################################################


"""
    reflect(diag::TensorDiagram; dir="vertical")

Reflects a `TensorDiagram` in the specified direction.
This operation involves:
1. Creating a copy of the diagram.
2. Calling `reflect` on all constituent nodes.
3. Swapping the boundary legs between opposing sides (top/bottom if vertical, left/right if horizontal).
4. Reversing the position indices of the boundary legs on the other pair of sides.
5. Updating the coordinates of the nodes if they are present.

# Arguments
- `diag::TensorDiagram`: The diagram to reflect.
- `dir::String="vertical"`: The direction of reflection. Can be "vertical" or "horizontal".

# Returns
- A new `TensorDiagram` representing the reflected diagram.
"""
function reflect(diag::TensorDiagram; dir="vertical")
    diag_r = copy(diag)

    diag_r.nodes .= reflect.(diag.nodes; dir=dir)

    flipping_sides = dir == "vertical" ? ["top", "bottom"] : ["left", "right"]
    reversing_sides = dir == "vertical" ? ["left", "right"] : ["top", "bottom"]

    diag_r.boundary_slots_num[flipping_sides[1]] = get(diag.boundary_slots_num, flipping_sides[2], 0)
    diag_r.boundary_slots_num[flipping_sides[2]] = get(diag.boundary_slots_num, flipping_sides[1], 0)
    diag_r.boundary_legs[flipping_sides[1]] = get(diag.boundary_legs, flipping_sides[2], Int[])
    diag_r.boundary_legs[flipping_sides[2]] = get(diag.boundary_legs, flipping_sides[1], Int[])
    diag_r.boundary_legs_posidx[flipping_sides[1]] = get(diag.boundary_legs_posidx, flipping_sides[2], Int[])
    diag_r.boundary_legs_posidx[flipping_sides[2]] = get(diag.boundary_legs_posidx, flipping_sides[1], Int[])


    rev_sides_num_1 = get(diag.boundary_slots_num, reversing_sides[1], 0)
    rev_sides_num_2 = get(diag.boundary_slots_num, reversing_sides[2], 0)

    diag_r.boundary_legs_posidx[reversing_sides[1]] = (rev_sides_num_1 + 1) .- get(diag.boundary_legs_posidx, reversing_sides[1], Int[])
    diag_r.boundary_legs_posidx[reversing_sides[2]] = (rev_sides_num_2 + 1) .- get(diag.boundary_legs_posidx, reversing_sides[2], Int[])

    if !isnothing(diag.node_coordinates)
        diag_r.node_coordinates .= (diag.node_coordinates) .|> x -> ismissing(x) ? missing : dir == "vertical" ? [x[1], -x[2]] : [-x[1], x[2]]
    end

    return diag_r
end