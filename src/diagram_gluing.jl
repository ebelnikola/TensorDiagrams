#################################################
# TOC
# topic 1:           GLUING
#################################################

function _opposed_side(side::String)
    if side == "left"
        return "right"
    elseif side == "right"
        return "left"
    elseif side == "top"
        return "bottom"
    elseif side == "bottom"
        return "top"
    else
        throw(ArgumentError("Invalid side: $side. Supported sides are 'left', 'right', 'top', 'bottom'."))
    end
end

"""
    glue_diagrams(d1::TensorDiagram, d2::TensorDiagram, sides::Vector{String})

Glue two `TensorDiagram`s along the specified sides.

# Arguments
- `d1::TensorDiagram`: The first diagram.
- `d2::TensorDiagram`: The second diagram.
- `sides::Vector{String}`: A two-element vector specifying which sides to glue. 
  For example, `["right", "left"]` glues the right side of `d1` to the left side of `d2`.
  Supported values: "left", "right", "top", "bottom".

# Returns
- `TensorDiagram`: The resulting glued diagram.
"""
function glue_diagrams(d1::TensorDiagram, d2::TensorDiagram, sides::Vector{String}; reflect_second=false)

    # 0. Making copies so that we can mutate diagrams in place
    d1 = copy(d1)
    d2 = copy(d2)

    if reflect_second
        d2 = reflect(d2; dir=sides[2] in ["left", "right"] ? "vertical" : "horizontal")
    end

    # 1. Geometry and compatibility:
    #    diagrams orientation, 
    #    spaces compatibility, 
    #    dangling indices dropping,
    #    nodes positioning 

    # 1.1 Orientation

    sides_order = ["left", "top", "right", "bottom"]
    orientation_1 = findfirst(==(sides[1]), sides_order)
    orientation_2 = findfirst(==(sides[2]), sides_order)

    if mod(orientation_1 - orientation_2, 4) != 2
        d2 = rotate(d2, 2 + (orientation_2 - orientation_1))
        sides[2] = _opposed_side(sides[1])
    end

    ishorizontal_1 = sides[1] in ["left", "right"]

    # 1.2 Spaces compatibility check

    space_id_1 = get_space_id(d1; side=sides[1])
    space_id_2 = get_space_id(d2; side=sides[2])

    if space_id_1 != space_id_2
        throw(ArgumentError("Cannot glue $(sides[1]) and $(sides[2]) sides of diagrams: space IDs do not match (d1 has $space_id_1, d2 has $space_id_2)."))
    end

    # 1.3 Dropping dangling indices on the gluing sides
    d1 = drop_dangling_indices!(d1, sides[1])
    d2 = drop_dangling_indices!(d2, sides[2])

    # 1.4 Nodes positioning

    shift_x_1, shift_y_1 = 0.0, 0.0
    shift_x_2, shift_y_2 = 0.0, 0.0

    if ishorizontal_1
        origin_shift_1 = sides[1] == "right" ? d2.width / 2 : -d2.width / 2
        origin_shift_2 = sides[2] == "left" ? -d1.width / 2 : d1.width / 2

        shift_x_1 = -origin_shift_1
        shift_x_2 = -origin_shift_2
    else
        origin_shift_1 = sides[1] == "top" ? d2.height / 2 : -d2.height / 2
        origin_shift_2 = sides[2] == "bottom" ? -d1.height / 2 : d1.height / 2

        shift_y_1 = -origin_shift_1
        shift_y_2 = -origin_shift_2
    end

    apply_shift(coords, sx, sy) = isnothing(coords) ? nothing : [(ismissing(c) ? missing : c .+ [sx, sy]) for c in coords]

    coords1 = apply_shift(d1.node_coordinates, shift_x_1, shift_y_1)
    coords2 = apply_shift(d2.node_coordinates, shift_x_2, shift_y_2)

    new_coordinates = if isnothing(coords1) && isnothing(coords2)
        nothing
    elseif isnothing(coords1)
        vcat(fill(missing, length(d1.nodes)), coords2)
    elseif isnothing(coords2)
        vcat(coords1, fill(missing, length(d2.nodes)))
    else
        vcat(coords1, coords2)
    end

    # 2. Relabeling diagrams 

    # 2.1 Creating generators of new indices     
    pattern_ind_of_d1 = unique(vcat(d1.contraction_pattern...))
    boundary_ind_of_d1 = unique(vcat(values(d1.boundary_legs)...))
    all_ind_of_d1 = unique(vcat(pattern_ind_of_d1, boundary_ind_of_d1))


    new_pos_index_cnt = [max(maximum(all_ind_of_d1), 0) + 1]
    new_neg_index_cnt = [min(minimum(all_ind_of_d1), 0) - 1]

    npi() = generate_nfi!(new_pos_index_cnt; dir="positive")
    nni() = generate_nfi!(new_neg_index_cnt; dir="negative")

    # 2.2. ordering legs on the gluing sides to match each other
    perm = sortperm(d1.boundary_legs_posidx[sides[1]])
    d1.boundary_legs[sides[1]] = d1.boundary_legs[sides[1]][perm]
    d1.boundary_legs_posidx[sides[1]] = d1.boundary_legs_posidx[sides[1]][perm]

    perm = sortperm(d2.boundary_legs_posidx[sides[2]])
    d2.boundary_legs[sides[2]] = d2.boundary_legs[sides[2]][perm]
    d2.boundary_legs_posidx[sides[2]] = d2.boundary_legs_posidx[sides[2]][perm]

    # 2.3. creating relabeling dicts
    #      we need two dictionaries, one for each diagram 
    #      for diagram d1 we relabel indices on the gluing side only, so 
    #      we initialize the dictionary with idx=>idx for the rest of indices
    #      for diagram d2 we relabel all indices to avoid conflicts, so
    #      we initialize it empty
    gluing_side_inds_d1 = d1.boundary_legs[sides[1]]
    gluing_side_inds_d2 = d2.boundary_legs[sides[2]]

    all_ind_of_d1_except_gluing = setdiff(all_ind_of_d1, gluing_side_inds_d1)

    relabeling_dict_d1 = Dict{Int,Int}(idx => idx for idx in all_ind_of_d1_except_gluing)
    relabeling_dict_d2 = Dict{Int,Int}()
    for (old_idx_1, old_idx_2) in zip(gluing_side_inds_d1, gluing_side_inds_d2)
        new_idx = npi()
        relabeling_dict_d1[old_idx_1] = new_idx
        relabeling_dict_d2[old_idx_2] = new_idx
    end

    pattern_ind_of_d2 = unique(vcat(d2.contraction_pattern...))
    boundary_ind_of_d2 = unique(vcat(values(d2.boundary_legs)...))
    all_ind_of_d2 = unique(vcat(pattern_ind_of_d2, boundary_ind_of_d2))

    for idx in all_ind_of_d2
        if !haskey(relabeling_dict_d2, idx)
            new_idx = idx > 0 ? npi() : nni()
            relabeling_dict_d2[idx] = new_idx
        end
    end

    # 2.4 Applying the relabeling dicts 
    relabel!(d1, relabeling_dict_d1)
    relabel!(d2, relabeling_dict_d2)

    # 3. Creating the new diagram 

    # 3.1 New nodes and contraction pattern 
    new_nodes = vcat(d1.nodes, d2.nodes)
    new_contraction_pattern = vcat(d1.contraction_pattern, d2.contraction_pattern)

    # 3.2 New boundary legs dictionaries 

    # 3.2.1. removing positive indices from old boundary dictionaries

    delete!(d1.boundary_legs, sides[1])
    delete!(d2.boundary_legs, sides[2])

    delete!(d1.boundary_legs_posidx, sides[1])
    delete!(d2.boundary_legs_posidx, sides[2])

    # 3.2.2 constructing new boundary dicts



    d1_pos = sides[2] # At this stage sides are necessarily different 
    d2_pos = sides[1] # we took care of this above by rotations 
    #                   Now, imagine sides = left, right, which means that 
    #                   we glue left boundary of d1 with right boundary of d2
    #                   So, d1_pos will be in the right part of the diagram (sides[2]), 
    #                   and d2_pos will be in the left part of the diagram (sides[1]) 


    if ishorizontal_1
        # Horizontal gluing: d1 and d2 are side-by-side
        # New "vertical" slots (top/bottom) are sum of previous
        new_slots_top = get(d1.boundary_slots_num, "top", 0) + get(d2.boundary_slots_num, "top", 0)
        new_slots_bottom = get(d1.boundary_slots_num, "bottom", 0) + get(d2.boundary_slots_num, "bottom", 0)

        # New "horizontal" slots (left/right) come from the outer boundaries

        is_d1_left = (d1_pos == "left")
        # If d1 is on the left, new left is d1's left boundary, new right is d2's right boundary

        new_slots_left = is_d1_left ? get(d1.boundary_slots_num, "left", 0) : get(d2.boundary_slots_num, "left", 0)
        new_slots_right = is_d1_left ? get(d2.boundary_slots_num, "right", 0) : get(d1.boundary_slots_num, "right", 0)

    else
        # Vertical gluing: d1 and d2 are stacked
        # New "horizontal" slots (left/right) are sum of previous
        new_slots_left = get(d1.boundary_slots_num, "left", 0) + get(d2.boundary_slots_num, "left", 0)
        new_slots_right = get(d1.boundary_slots_num, "right", 0) + get(d2.boundary_slots_num, "right", 0)

        # New "vertical" slots (top/bottom) come from outer boundaries
        is_d1_bottom = (d1_pos == "bottom")
        # If d1 is on the bottom, new bottom is d1's bottom boundary, new top is d2's top boundary

        new_slots_bottom = is_d1_bottom ? get(d1.boundary_slots_num, "bottom", 0) : get(d2.boundary_slots_num, "bottom", 0)
        new_slots_top = is_d1_bottom ? get(d2.boundary_slots_num, "top", 0) : get(d1.boundary_slots_num, "top", 0)
    end

    new_boundary_slots_num = Dict(
        "top" => new_slots_top,
        "bottom" => new_slots_bottom,
        "left" => new_slots_left,
        "right" => new_slots_right
    )

    new_boundary_legs = Dict{String,Vector{Int}}()
    new_boundary_legs_posidx = Dict{String,Vector{Int}}()


    # Now, if we do horizontal gluing (dL dR) we should shift the posidx on 
    # top and bottom sides for the diagram that will be on the right side.
    # The shift amount is the number of top/bottom slots of the left diagram.

    for side in ["top", "bottom", "left", "right"]
        legs_1 = get(d1.boundary_legs, side, Int[])
        legs_2 = get(d2.boundary_legs, side, Int[])
        posidx_1 = get(d1.boundary_legs_posidx, side, Int[])
        posidx_2 = get(d2.boundary_legs_posidx, side, Int[])

        current_shift_1 = 0
        current_shift_2 = 0

        if ishorizontal_1 && (side == "top" || side == "bottom")
            # Horizontal Gluing. 
            # If d1 is on the LEFT (d1_pos == "left"), d2 is on the RIGHT. 
            # d2 needs to be shifted by d1's width (slots)
            if d1_pos == "left"
                current_shift_2 = get(d1.boundary_slots_num, side, 0)
            end

            # If d2 is on the LEFT (d2_pos == "left"), d1 is on the RIGHT.
            # d1 needs to be shifted by d2's width (slots)
            if d2_pos == "left"
                current_shift_1 = get(d2.boundary_slots_num, side, 0)
            end
        elseif !ishorizontal_1 && (side == "left" || side == "right")
            # Vertical Gluing.
            # If d1 is on the BOTTOM (d1_pos == "bottom"), d2 is on the TOP.
            # d2 needs to be shifted by d1's height (slots)
            if d1_pos == "bottom"
                current_shift_2 = get(d1.boundary_slots_num, side, 0)
            end
            # If d2 is on the BOTTOM (d2_pos == "bottom"), d1 is on the TOP.
            # d1 needs to be shifted by d2's height (slots)
            if d2_pos == "bottom"
                current_shift_1 = get(d2.boundary_slots_num, side, 0)
            end
        end

        # Apply shifts
        # Note: If side is not one of the shifted ones, shift is 0.
        posidx_1_shifted = [pos + current_shift_1 for pos in posidx_1]
        posidx_2_shifted = [pos + current_shift_2 for pos in posidx_2]

        new_boundary_legs[side] = vcat(legs_1, legs_2)
        new_boundary_legs_posidx[side] = vcat(posidx_1_shifted, posidx_2_shifted)
    end

    # 3.3 New labels dict
    new_labels = merge(d1.labels, d2.labels)

    # 3.4 New factor

    new_factor = d1.factor == "1.0" ? d2.factor : (d2.factor == "1.0" ? d1.factor : "(" * d1.factor * "*" * d2.factor * ")")

    # 3.5 New node coordinates 

    # already computed and stored in new_coordinates

    # 3.6 New width and height

    new_width = ishorizontal_1 ? d1.width + d2.width : max(d1.width, d2.width)
    new_height = ishorizontal_1 ? max(d1.height, d2.height) : d1.height + d2.height

    # 4. return
    return TensorDiagram(
        nodes=new_nodes,
        contraction_pattern=new_contraction_pattern,
        boundary_legs=new_boundary_legs,
        boundary_legs_posidx=new_boundary_legs_posidx,
        boundary_slots_num=new_boundary_slots_num,
        labels=new_labels,
        factor=new_factor,
        node_coordinates=new_coordinates,
        width=new_width,
        height=new_height
    )
end

