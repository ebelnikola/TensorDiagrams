#################################################
# TOC
# topic 0:           EXPORTS
# topic 1:           DIAGRAM TO CODE
#################################################


#########################################################
# topic 0: EXPORTS
#########################################################
export diagram_to_code

#########################################################
# topic 1: DIAGRAM TO CODE
#########################################################

"""
    form_name_with_flips(node::TensorNode)

Constructs the name of the tensor variable corresponding to a given `TensorNode`.
If the node has flipped indices (indicated by `node.flipped`), a suffix is appended to the base name.
The suffix is of the form "_flp" followed by the indices in `node.flipped` joined by underscores.

# Arguments
- `node::TensorNode`: The tensor node.

# Returns
- A `String` representing the tensor name.
"""
function form_name_with_flips(node::TensorNode)
    tensor_name = node.name

    flipped_legs_ind = findall(node.flipped)

    if !isempty(flipped_legs_ind)
        if endswith(tensor_name, "'")
            tensor_name = tensor_name[1:end-1]
            tensor_name *= "_flp_" * join(flipped_legs_ind, "_")
            tensor_name *= "'"
        else
            tensor_name *= "_flp_" * join(flipped_legs_ind, "_")
        end
    end

    return tensor_name
end

"""
    diagram_to_code(diag::TensorDiagram; result_name="res", order::Array{Int}=Int[], sides_order=["left", "top", "right", "bottom"])

Generates a string representation of the contraction from the diagram using @tensoropt.

# Arguments:
- `diag::TensorDiagram`: The tensor diagram to convert to code.
- `result_name::String="res"`: The name of the result variable in the generated code.
- `order::Array{Int}=Int[]`: The order of the output legs. If empty, defaults to the order specified by `sides_order`.
- `sides_order::Array{String}=["left", "top", "right", "bottom"]`: The order of sides to use when constructing the default output leg order.

# Returns:
- `code::String`: A string representing the tensor contraction code.
"""
function diagram_to_code(diag::TensorDiagram; result_name="res", order::Array{Int}=Int[], sides_order=["left", "top", "right", "bottom"], semicolon_pos=0, perm_side_based_order::Vector{Int}=Int[])
    code = "@tensoropt $result_name"

    # 1. Completing the order array if needed. 
    #    The order of legs in the final contraction, 
    #    if not provided, is:
    #
    #    Left boundary legs from bottom to top, 
    #    Top boundary legs from left to right,
    #    Right boundary legs from bottom to top,
    #    Bottom boundary legs from left to right, 
    #    
    #    If custom sides_order is provided, the order of boundaries
    #    will be taken from there, but the order of legs on each side 
    #    remains the same as described above.

    if length(order) == 0
        groups = [
            get(diag.boundary_legs, sides_order[i], Int[])[
                sortperm(get(diag.boundary_legs_posidx, sides_order[i], Int[]))
            ]
            for i in eachindex(sides_order)
        ]
        order = vcat(groups...)
        if !isempty(perm_side_based_order)
            order = order[perm_side_based_order]
        end
    end

    # 2. Adding the ordered indices to the left hand side of the code string
    if length(order) > 0
        if semicolon_pos == 0
            code *= "[" * join(order, ", ") * "]"
        else
            code *= "[" * join(order[1:semicolon_pos], ", ") * "; " * join(order[semicolon_pos+1:end], ", ") * "]"
        end
    end

    # 3. Adding the right hand side
    code *= ":="

    for i in eachindex(diag.contraction_pattern)
        node = diag.nodes[i]
        tensor_name = form_name_with_flips(node)
        pattern = diag.contraction_pattern[i]
        code *= i > 1 ? " * $tensor_name" : "$tensor_name"
        if length(pattern) > 0
            code *= "[" * join(pattern, ", ") * "]"
        end
    end

    # 4. Factor multiplication
    code *= "*(" * diag.factor * ")"

    return code
end
