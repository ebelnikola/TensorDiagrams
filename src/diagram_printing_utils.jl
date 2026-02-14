"""
    pretty_format(diag::TensorDiagram; pad_to::Int=0)

Returns a formatted string representation of a `TensorDiagram`.
Aligns the output keys for better readability and improves dictionary visualization.

# Arguments
- `diag`: The `TensorDiagram` to format.
- `pad_to`: Minimum height (number of newlines) of the output string. Padding is distributed roughly equally.
"""
function pretty_format(diag::TensorDiagram; pad_to::Int=0)
    io = IOBuffer()

    println(io, "TensorDiagram:")

    # 0. Consistency Check
    consistency_results = check_diagram_consistency(diag)
    is_consistent = all(values(consistency_results))

    if is_consistent
        println(io, "  Status: Consistent")
    else
        println(io, "  Status: Inconsistent! Run check_diagram_consistency(diag) for details.")
    end

    # 1. Nodes & Contraction Pattern
    println(io, "  Nodes & Contraction Pattern:")
    node_names = [n.name for n in diag.nodes]
    pattern_strs = [string(p) for p in diag.contraction_pattern]
    _print_aligned_columns(io, [node_names, pattern_strs])

    # 2. Boundary Legs Num
    h_num = get(diag.boundary_slots_num, "horizontal", 0)
    v_num = get(diag.boundary_slots_num, "vertical", 0)
    println(io, "  Hor. boundary legs num: ", h_num)
    println(io, "  Ver. boundary legs num: ", v_num)
    println(io, "")

    # 3. Boundary Structure (Table only)
    println(io, "  Boundary Structure:")
    sides = ["left", "top", "right", "bottom"]
    valid_sides = [side for side in sides if haskey(diag.boundary_legs, side) && !isempty(diag.boundary_legs[side])]

    if isempty(valid_sides)
        println(io, "    (No boundary legs defined)\n")
    else
        side_col = valid_sides
        legs_col = [string(diag.boundary_legs[s]) for s in valid_sides]
        pos_col = [haskey(diag.boundary_legs_posidx, s) ? string(diag.boundary_legs_posidx[s]) : "-" for s in valid_sides]

        _print_aligned_columns(io, [side_col, legs_col, pos_col], ["Side:", "Legs:", "Positions:"])
    end

    # 4. Labels Table
    println(io, "  Labels:")
    sorted_keys = sort(collect(keys(diag.labels)))
    if isempty(sorted_keys)
        println(io, "    (No labels)\n")
    else
        lbl_col = [string(k) for k in sorted_keys]
        val_col = [string(diag.labels[k]) for k in sorted_keys]
        _print_aligned_columns(io, [lbl_col, val_col], ["Pattern Index:", "Label:"])
    end

    # 5. Factor
    println(io, "  Factor: ", diag.factor)
    println(io, "")

    # 6. Node Coordinates
    if !isnothing(diag.node_coordinates) && !isempty(diag.node_coordinates) && length(diag.node_coordinates) == length(diag.nodes)
        println(io, "  Node Coordinates:")
        node_data = []
        missing_nodes = String[]
        for (i, node) in enumerate(diag.nodes)
            c = diag.node_coordinates[i]

            if !ismissing(c)
                coord_str = string(c) # Representation exactly as is
                # Extract coordinates for sorting
                x = length(c) >= 1 ? c[1] : 0.0
                y = length(c) >= 2 ? c[2] : 0.0
                push!(node_data, (y=y, x=x, str="$(node.name)$coord_str"))
            else
                push!(missing_nodes, "$(node.name)[-,-]")
            end
        end

        # Sort by Y descending, then X ascending
        sort!(node_data, by=d -> (-d.y, d.x))

        if !isempty(node_data)
            current_y = node_data[1].y
            line_nodes = String[]

            for data in node_data
                if !isapprox(data.y, current_y; atol=1e-5)
                    println(io, "    " * join(line_nodes, " "))
                    current_y = data.y
                    empty!(line_nodes)
                end
                push!(line_nodes, data.str)
            end
            println(io, "    " * join(line_nodes, " "))
        end

        if !isempty(missing_nodes)
            println(io, "    " * join(missing_nodes, " "))
        end
    else
        if !isnothing(diag.node_coordinates) && !isempty(diag.node_coordinates)
            # This means coordinates exist but length does not match nodes
            println(io, "  Node Coordinates: [Inconsistent data (length mismatch between nodes and node_coordinates)]")
        else
            println(io, "  Node Coordinates: None")
        end
    end
    println(io, "")

    # 7. Dimensions
    println(io, "  Dimensions: width=$(diag.width), height=$(diag.height)")

    content = String(take!(io))
    num_newlines = count(==('\n'), content)

    if num_newlines < pad_to
        diff = pad_to - num_newlines
        top = div(diff, 2)
        bottom = diff - top
        return repeat("\n", top) * content * repeat("\n", bottom)
    else
        return content
    end
end

"""
    pretty_print(diag::TensorDiagram; kwargs...)

Prints a formatted string representation of a `TensorDiagram` to standard output.
Passes all keyword arguments to `pretty_format`.
"""
function pretty_print(diag::TensorDiagram; kwargs...)
    print(pretty_format(diag; kwargs...))
end

"""
    pretty_format(node::TensorNode; pad_to::Int=0)

Returns a formatted string representation of a `TensorNode`.
Aligns the output keys for better readability and improves dictionary visualization.

# Arguments
- `node`: The `TensorNode` to format.
- `pad_to`: Minimum height (number of newlines) of the output string. Padding is distributed roughly equally.
"""
function pretty_format(node::TensorNode; pad_to::Int=0)
    io = IOBuffer()

    println(io, "TensorNode: \"$(node.name)\"")

    # 0. Consistency Check
    consistency_results = check_node_consistency(node)
    is_consistent = all(values(consistency_results))

    if is_consistent
        println(io, "  Status: Consistent")
    else
        println(io, "  Status: Inconsistent! Run check_node_consistency(node) for details.")
    end

    # 1. Legs Structure
    println(io, "  Legs:")
    sides = ["left", "top", "right", "bottom"]
    valid_sides = [side for side in sides if haskey(node.legs, side) && !isempty(node.legs[side])]

    if isempty(valid_sides)
        println(io, "    (No legs defined)\n")
    else
        side_col = valid_sides
        legs_col = [string(node.legs[s]) for s in valid_sides]
        _print_aligned_columns(io, [side_col, legs_col], ["Side:", "Indices:"])
    end

    # 2. Allowed Labels
    println(io, "  Allowed Labels:")
    if isempty(node.allowed_labels)
        println(io, "    (No allowed labels info)\n")
    else
        idx_col = [string(i) for i in 1:length(node.allowed_labels)]
        lbls_col = [string(lbls) for lbls in node.allowed_labels]
        _print_aligned_columns(io, [idx_col, lbls_col], ["Leg Index:", "Allowed Labels:"])
    end

    # 3. Label Combinations
    if !isempty(node.allowed_label_combinations)
        println(io, "  Allowed Label Combinations:")
        comb_col = [join(c, " ") for c in node.allowed_label_combinations]
        _print_aligned_columns(io, [comb_col], ["Combination:"])
    end

    if !isempty(node.forbidden_label_combinations)
        println(io, "  Forbidden Label Combinations:")
        comb_col = [join(c, " ") for c in node.forbidden_label_combinations]
        _print_aligned_columns(io, [comb_col], ["Combination:"])
    end

    # 4. Properties
    println(io, "  Properties:")
    props_row1 = ["reflections", "adjoint", "out_legs_num", "dual", "flipped"]
    props_row2 = [
        "hor: $(node.hor_ref)",
        string(node.adjoint),
        string(node.out_legs_num),
        string(Int.(node.dual)),
        string(Int.(node.flipped))
    ]
    props_row3 = [
        "ver: $(node.ver_ref)",
        "",
        "",
        "",
        ""
    ]
    _print_aligned_columns(io, [props_row1, props_row2, props_row3])

    # 5. Visuals
    println(io, "  Visuals:")
    println(io, "    width: $(isnothing(node.width) ? "auto" : node.width)")
    println(io, "    height: $(isnothing(node.height) ? "auto" : node.height)")
    c = node.color
    println(io, "    color: RGBA($(round(c.r, digits=3)), $(round(c.g, digits=3)), $(round(c.b, digits=3)), $(round(c.alpha, digits=3)))")

    content = String(take!(io))
    num_newlines = count(==('\n'), content)

    if num_newlines < pad_to
        diff = pad_to - num_newlines
        top = div(diff, 2)
        bottom = diff - top
        return repeat("\n", top) * content * repeat("\n", bottom)
    else
        return content
    end
end

"""
    pretty_print(node::TensorNode; kwargs...)

Prints a formatted string representation of a `TensorNode` to standard output.
Passes all keyword arguments to `pretty_format`.
"""
function pretty_print(node::TensorNode; kwargs...)
    print(pretty_format(node; kwargs...))
end

function _print_aligned_columns(io::IO, cols::Vector{Vector{String}}, headers::Vector{String}=String[])
    isempty(cols) && return
    num_rows = length(cols[1])
    num_cols = length(cols)
    if num_rows == 0
        return
    end

    # Calculate width per item in the list
    col_widths = Int[] # Width per LIST ITEM index
    for i in 1:num_rows
        w = 0
        for col_idx in 1:num_cols
            w = max(w, length(cols[col_idx][i]))
        end
        push!(col_widths, w + 2)
    end

    # Calculate wrapped lines
    max_line_width = 100
    current_width = 0
    chunk_indices = Int[]

    for i in 1:num_rows
        current_width += col_widths[i]
        if current_width > max_line_width
            push!(chunk_indices, i - 1)
            current_width = col_widths[i]
        end
    end
    push!(chunk_indices, num_rows)

    current_start = 1
    for end_idx in chunk_indices
        if end_idx < current_start
            continue
        end
        range = current_start:end_idx

        # Print each type of column (header row)
        for col_idx in 1:num_cols
            prefix = "    "
            if col_idx <= length(headers)
                prefix *= rpad(headers[col_idx], 16)
            else
                prefix *= "    "
            end

            row_str = prefix
            for i in range
                row_str *= rpad(cols[col_idx][i], col_widths[i])
            end
            println(io, row_str)
        end
        println(io, "")
        current_start = end_idx + 1
    end
end
