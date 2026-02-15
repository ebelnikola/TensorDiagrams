#################################################
# TOC
# topic 1:           PLOTTING LEG MARKERS ALONG SEGMENTS
# topic 2:           PLOTTING TENSOR NODES
# topic 3:           PLOTTING TENSOR DIAGRAMS
# topic 4:           DISPLAY UTILITIES
#################################################

using CairoMakie
using NetworkLayout

#################################################
# topic 1: PLOTTING LEG MARKERS ALONG SEGMENTS  
#################################################

"""
    draw_legs_on_segment!(ax, start_coordinates, end_coordinates, n_slots_total, leg_labels, activity; marker_size, fontsize, strokewidth=2)

Draw legs along a line segment with active legs showing labels and inactive legs as gray squares.

Arguments

Positional:
- ax: The Makie axis to draw on
- start_coordinates: Tuple (x, y) for the start position of the segment
- end_coordinates: Tuple (x, y) for the end position of the segment
- n_slots_total: Total number of leg positions along the segment
- leg_labels: Array of labels for active legs (e.g., [1, 2, 3])
- activity: Boolean array of length n_slots_total indicating whether each position is active (e.g., [true, false, true, false])

Keyword:
- marker_size: Size of the leg markers. If `nothing` (default), it is calculated automatically based on `fontsize` and the length of leg labels.
- fontsize: Font size for the leg labels (default: 14)
- strokewidth: Width of the marker stroke (default: 2)
- leg_background_colors: Optional vector of colors (RGBA) to use for the active leg markers. If empty (default), colors are determined by activity (active=white, inactive=gray).

Returns:
- nothing
"""
function draw_legs_on_segment!(ax, start_coordinates, end_coordinates, n_slots_total, leg_labels, activity; perm=collect(1:length(leg_labels)), marker_size=nothing, fontsize=14, strokewidth=2, leg_background_colors=RGBA[])
    x_start, y_start = start_coordinates
    x_end, y_end = end_coordinates

    t = collect(1:n_slots_total) ./ (n_slots_total + 1)
    x = x_start .+ t * (x_end - x_start)
    y = y_start .+ t * (y_end - y_start)
    if isempty(leg_background_colors)
        colors = (x -> x ? :white : :gray).(activity)
    else
        colors = leg_background_colors
    end

    # Calculate marker_size if not provided
    # We assume fontsize is in the same units as marker_size (screen units by default).
    # A character typically has width ~0.6 * fontsize.
    # A standard square marker with side 2.0 * fontsize provides comfortable padding for 1-2 digits.
    # For longer strings, we scale up.
    ms = if isnothing(marker_size)
        max_chars = isempty(leg_labels) ? 0 : maximum(length.(string.(leg_labels)))
        max(2.0 * fontsize, max_chars * fontsize)
    else
        marker_size
    end

    scatter!(ax, x, y, color=colors, markersize=ms,
        marker=:rect, strokecolor=:black, strokewidth=strokewidth)
    text!(ax, x[activity][perm], y[activity][perm]; text=string.(leg_labels),
        align=(:center, :center), fontsize=fontsize)
    return nothing
end

#################################################
# topic 2: PLOTTING TENSOR NODES 
#################################################

"""
    node_aspect_ratio(node::TensorNode)

Calculate the aspect ratio (height/width) for a TensorNode based on its leg configuration.

Arguments

Positional:
- node: The TensorNode to calculate aspect ratio for

Returns:
- Float64: The aspect ratio (height/width) = 1.0 if no legs, golden ratio φ if only horizontal legs, 1/φ if only vertical legs, or horizontal_legs/vertical_legs otherwise.
"""
function node_aspect_ratio(node::TensorNode)
    if !isnothing(node.width) && !isnothing(node.height)
        return node.height / node.width
    end
    # Count legs on each side (handle missing keys)
    n_left = haskey(node.legs, "left") ? length(node.legs["left"]) : 0
    n_right = haskey(node.legs, "right") ? length(node.legs["right"]) : 0
    n_top = haskey(node.legs, "top") ? length(node.legs["top"]) : 0
    n_bottom = haskey(node.legs, "bottom") ? length(node.legs["bottom"]) : 0

    # Calculate dimensions
    vertical_legs = max(n_top, n_bottom)
    horizontal_legs = max(n_left, n_right)
    φ = (1 + sqrt(5)) / 2

    # Aspect ratio: height / width = horizontal_legs / vertical_legs
    if vertical_legs == 0 && horizontal_legs == 0
        aspect_ratio = 1.0
    elseif vertical_legs == 0
        aspect_ratio = φ  # Use golden ratio if no vertical legs
    elseif horizontal_legs == 0
        aspect_ratio = 1 / φ  # Use inverse golden ratio if no horizontal legs
    else
        aspect_ratio = horizontal_legs / vertical_legs
    end

    return aspect_ratio
end

"""
    default_node_width(node::TensorNode)

Calculate the default width for a TensorNode based on its properties.

Arguments:
- node: The TensorNode

Returns:
- Float64: The calculated width
"""
function default_node_width(node::TensorNode)
    aspect_ratio = node_aspect_ratio(node)
    vert_legs_num = max(length(get(node.legs, "top", Int[])), length(get(node.legs, "bottom", Int[])))

    if !isnothing(node.width)
        return node.width
    elseif !isnothing(node.height)
        return node.height / aspect_ratio
    else
        return float(vert_legs_num + 2) * 0.55
    end
end


"""
    draw_tensor_node!(ax, node::TensorNode, center_x=0.0, center_y=0.0, width=nothing; leg_fontsize=14, strokewidth=2, name_fontsize=28)

Draw a TensorNode on the given axis at a specified position and scale. 
Leg colors are determined by the `dual` and `flipped` fields of the node (green for "in", red for "out" / dual).

Arguments

Positional:
- ax: The Makie axis to draw on
- node: The TensorNode to draw
- center_x: X coordinate of the node center (default: 0.0)
- center_y: Y coordinate of the node center (default: 0.0)
- width: Width of the node (default: 1.0)

Keyword:
- leg_fontsize: Font size for the leg labels (default: 14)
- strokewidth: Width of the marker stroke (default: 2)
- name_fontsize: Font size for tensor name (default: 28)

Returns:
- nothing
"""
function draw_tensor_node!(ax, node::TensorNode, center_x=0.0, center_y=0.0, width=nothing; leg_fontsize=14, strokewidth=2, name_fontsize=28)

    if isnothing(width)
        width = default_node_width(node)
    end

    aspect_ratio = node_aspect_ratio(node)
    height = width * aspect_ratio

    # Rectangle corners (centered at given position)
    node_x_min = center_x - width / 2
    node_x_max = center_x + width / 2
    node_y_min = center_y - height / 2
    node_y_max = center_y + height / 2

    # Draw filled rectangle with black outline
    poly!(ax,
        Point2f[(node_x_min, node_y_min), (node_x_max, node_y_min),
            (node_x_max, node_y_max), (node_x_min, node_y_max)],
        color=node.color,
        strokecolor=:black,
        strokewidth=strokewidth
    )

    # Draw legs on the node boundaries
    sides_to_endpoints = Dict(
        "left" => [(node_x_min, node_y_min), (node_x_min, node_y_max)],
        "right" => [(node_x_max, node_y_min), (node_x_max, node_y_max)],
        "top" => [(node_x_min, node_y_max), (node_x_max, node_y_max)],
        "bottom" => [(node_x_min, node_y_min), (node_x_max, node_y_min)]
    )

    dual = copy(node.dual)
    dual[node.out_legs_num+1:end] .= .!dual[node.out_legs_num+1:end]
    flipped = node.flipped

    for side in ["left", "right", "top", "bottom"]
        legs = get(node.legs, side, Int[])
        n = length(legs)
        start_pos = sides_to_endpoints[side][1]
        end_pos = sides_to_endpoints[side][2]

        colors = Dict(
            (false, false) => PALE_GREEN,
            (true, false) => PALE_RED,
            (false, true) => PALE_RED,
            (true, true) => PALE_GREEN,
        )

        leg_bg_colors = [colors[(dual[i], flipped[i])] for i in legs]

        draw_legs_on_segment!(ax, start_pos, end_pos, n, legs, fill(true, n); fontsize=leg_fontsize, strokewidth=strokewidth, leg_background_colors=leg_bg_colors)
    end

    # Add tensor name in center
    text!(ax, center_x, center_y, text=node.name,
        align=(:center, :center), fontsize=name_fontsize, color=:black)

    # Add reflection arrows in top-right corner based on reflection type
    arrow_fontsize = name_fontsize * 0.618
    arrow_x = node_x_max
    arrow_y = node_y_max

    # Construct the label string
    label_text = ""
    if node.hor_ref
        label_text *= "↔"
    end
    if node.ver_ref
        label_text *= "↕"
    end

    # Add rotation symbol
    rot_k = mod(node.rotation, 4)
    if rot_k == 1
        label_text *= "↺"
    elseif rot_k == 2
        label_text *= "↺²"
    elseif rot_k == 3
        label_text *= "↺³"
    end

    if !isempty(label_text)
        text!(ax, arrow_x, arrow_y, text=label_text,
            align=(:right, :top), fontsize=arrow_fontsize, color=:black)
    end

    return nothing
end

"""
    plot_tensor_node(node::TensorNode)

Plot a TensorNode. Leg numbers are displayed as squares at the appropriate boundaries.

Arguments:
- node: The TensorNode to plot

Returns:
- A Makie figure containing the plot
"""
function plot_tensor_node(node::TensorNode; scale=1.0)
    aspect_ratio = node_aspect_ratio(node)
    n_vertical = max(length(get(node.legs, "top", Int[])), length(get(node.legs, "bottom", Int[])))
    n_horizontal = max(length(get(node.legs, "left", Int[])), length(get(node.legs, "right", Int[])))

    width = default_node_width(node)

    base_size = (!isnothing(node.width) || !isnothing(node.height)) ? 40 * scale : min(width / n_vertical, width * aspect_ratio / n_horizontal) * 25 * scale
    leg_fontsize = base_size * 0.7
    name_fontsize = base_size * 1.2

    fig = Figure(size=(400 * scale, aspect_ratio * 400 * scale), figure_padding=0)

    ax = Axis(fig[1, 1],
        aspect=1 / aspect_ratio,
        limits=(-0.8 * width, 0.8 * width, -aspect_ratio * width * 0.8, aspect_ratio * width * 0.8))

    draw_tensor_node!(ax, node, 0.0, 0.0, width;
        leg_fontsize=leg_fontsize,
        name_fontsize=name_fontsize,
        strokewidth=2 * scale)

    hidedecorations!(ax)
    hidespines!(ax)

    return fig
end

# Default plotting for TensorNode
function Base.show(io::IO, m::Union{MIME"image/png",MIME"image/svg+xml"}, node::TensorNode)
    fig = plot_tensor_node(node)
    show(io, m, fig)
end

"""
    plot_tensor_node_list(nodes; cols=3, scale=1.0, size=(900, 600))

Plots a list of tensor nodes in a grid layout.

# Arguments
- `nodes`: A list of tensor nodes to plot.
- `cols`: Number of columns in the grid (default: 3).
- `scale`: Overall scale factor for the nodes (default: 1.0).
- `size`: Size of the figure (default: (900, 600)).

# Returns
- A `Makie.Figure` object containing the plotted nodes.
"""
function plot_tensor_node_list(nodes; cols=3, scale=1.0, size=(900, 600))
    fig = Figure(size=size, figure_padding=0)
    for (i, node) in enumerate(nodes)
        r = div(i - 1, cols) + 1
        c = mod(i - 1, cols) + 1

        n_vertical = max(length(get(node.legs, "top", Int[])), length(get(node.legs, "bottom", Int[])))
        n_horizontal = max(length(get(node.legs, "left", Int[])), length(get(node.legs, "right", Int[])))
        aspect_ratio = node_aspect_ratio(node)

        width = default_node_width(node)

        # Consistent base_size calculation with plot_tensor_node
        base_size = (!isnothing(node.width) || !isnothing(node.height)) ? 40 * scale : min(width / n_vertical, width * aspect_ratio / n_horizontal) * 25 * scale

        leg_fontsize = base_size * 0.7
        name_fontsize = base_size * 1.2

        # Using similar limits logic as plot_tensor_node
        ax = Axis(fig[r, c], aspect=DataAspect(),
            limits=(-0.8 * width, 0.8 * width, -aspect_ratio * width * 0.8, aspect_ratio * width * 0.8))
        hidedecorations!(ax)
        hidespines!(ax)
        draw_tensor_node!(ax, node, 0.0, 0.0, width, leg_fontsize=leg_fontsize, name_fontsize=name_fontsize)
    end
    return fig
end



#################################################
# topic 3: PLOTTING TENSOR DIAGRAMS 
#################################################

"""
    calculate_node_coordinates(diagram::TensorDiagram)

If all node coordinates are provided in the diagram, returns them directly. Otherwise, computes missing node coordinates using NetworkLayout. Boundary leg positions and positions of fixed nodes are used as fixed anchor points to guide the layout.

Arguments:
- diagram: The TensorDiagram
- node_repulsion_strength: (Keyword) Strength of node repulsion used in the stress layout algorithm (default: 3.0).

Returns:
- Array of [x, y] vectors for each node position
"""
function calculate_node_coordinates(diagram::TensorDiagram; node_repulsion_strength=3.0)

    if !isnothing(diagram.node_coordinates) && !any(ismissing, diagram.node_coordinates)
        return diagram.node_coordinates
    end

    n_nodes = length(diagram.nodes)

    if n_nodes == 0
        return []
    elseif n_nodes == 1
        return [(0.0, 0.0)]
    end

    box_width = diagram.width
    box_height = diagram.height

    x_min = -box_width / 2
    x_max = box_width / 2
    y_min = -box_height / 2
    y_max = box_height / 2

    # BOUNDARY LEGS ORDER IN ADJ. MATRIX AND POSITIONS

    boundary_legs = Int[]
    boundary_legs_coordinates = Tuple{Float64,Float64}[]

    n_left = get(diagram.boundary_slots_num, "left", 0)
    n_right = get(diagram.boundary_slots_num, "right", 0)
    n_top = get(diagram.boundary_slots_num, "top", 0)
    n_bottom = get(diagram.boundary_slots_num, "bottom", 0)
    # Configuration for all boundaries: (side_name, n_total, box_dimension, fixed_coord, vary_axis, vary_min)
    boundary_configs = [
        ("left", n_left, box_height, x_min, "y", y_min),
        ("right", n_right, box_height, x_max, "y", y_min),
        ("top", n_top, box_width, y_max, "x", x_min),
        ("bottom", n_bottom, box_width, y_min, "x", x_min)
    ]
    for (side, n_total, box_size, fixed_coord, vary_axis, vary_min) in boundary_configs
        legs = get(diagram.boundary_legs, side, Int[])
        leg_posidx = get(diagram.boundary_legs_posidx, side, Int[])
        for (local_leg_idx, pos_idx) in enumerate(leg_posidx)
            varying_coord = vary_min + (pos_idx / (n_total + 1)) * box_size
            coord = vary_axis == "x" ? (varying_coord, fixed_coord) : (fixed_coord, varying_coord)
            push!(boundary_legs, legs[local_leg_idx])
            push!(boundary_legs_coordinates, coord)
        end
    end

    # ADJ MATRIX BUILDING
    n_boundary = length(boundary_legs)
    n_total = n_nodes + n_boundary
    adjacency_matrix = zeros(n_total, n_total)

    for j = 1:n_nodes
        for i = 1:n_nodes
            if any(diagram.contraction_pattern[i]' .∈ diagram.contraction_pattern[j])
                adjacency_matrix[i, j] = node_repulsion_strength
            else
                adjacency_matrix[i, j] = node_repulsion_strength * 2.0
            end
        end
    end

    boundary_leg_to_node_idx = Dict(leg => n_nodes + nd for (nd, leg) in enumerate(boundary_legs))
    for i = 1:n_nodes
        for leg in diagram.contraction_pattern[i]
            if haskey(boundary_leg_to_node_idx, leg)
                nd = boundary_leg_to_node_idx[leg]
                adjacency_matrix[i, nd] = 1.0
                adjacency_matrix[nd, i] = 1.0
            end
        end
    end

    # BOUNDARY POSITIONS PINNING
    pin_dict = Dict{Int,Point{2,Float64}}(nd + n_nodes => Point(pos[1], pos[2])
                                          for (nd, pos) in enumerate(boundary_legs_coordinates))

    if !isnothing(diagram.node_coordinates)
        for (i, coord) in enumerate(diagram.node_coordinates)
            if !ismissing(coord)
                pin_dict[i] = Point(coord[1], coord[2])
            end
        end
    end

    # CALCULATE LAYOUT
    layout_positions = stress(adjacency_matrix, pin=pin_dict)
    return layout_positions[1:n_nodes]
end

"""
    get_leg_idx_to_coordinates_dict(node, center_x, center_y, width, height, contraction_pattern)

Helper function to collect leg coordinates for a node's legs.

Arguments:
- node: The TensorNode to collect leg coordinates from
- center_x, center_y: Center coordinates of the node
- width, height: Dimensions of the node
- contraction_pattern: Vector of indices where `contraction_pattern[i]` is the  index for the i-th leg of the node.

Returns:
- leg_index_to_coordinates: Dict mapping leg indices to (x, y) coordinates
"""
function get_leg_idx_to_coordinates_dict(node, center_x, center_y, width, height, contraction_pattern)
    leg_index_to_coordinates = Dict{Int,Tuple{Float64,Float64}}()

    x_min = center_x - width / 2
    x_max = center_x + width / 2
    y_min = center_y - height / 2
    y_max = center_y + height / 2

    border_configs = Dict(
        "left" => (x_min, y_min, height, "y"),
        "right" => (x_max, y_min, height, "y"),
        "top" => (y_max, x_min, width, "x"),
        "bottom" => (y_min, x_min, width, "x")
    )

    leg_position_idx = Dict(leg => (border, pos_idx)
                            for (border, legs) in node.legs
                            for (pos_idx, leg) in enumerate(legs)) # tells which leg is at which border and at which position

    for (leg, pattern_idx) in enumerate(contraction_pattern)
        border, pos_idx = leg_position_idx[leg]
        fixed_coord, vary_min, vary_range, vary_axis = border_configs[border]
        n = length(node.legs[border])

        varying_coord = vary_min + (pos_idx / (n + 1)) * vary_range
        coord = vary_axis == "x" ? (varying_coord, fixed_coord) : (fixed_coord, varying_coord)

        leg_index_to_coordinates[leg] = coord
    end


    return leg_index_to_coordinates
end

"""
    plot_tensor_diagram(diagram::TensorDiagram; kwargs...)

Plot a TensorDiagram showing all tensor nodes and their contractions.
- Nodes are positioned using NetworkLayout.jl
- Contracted indices are shown as lines between nodes
- External indices connect to boundary markers

Arguments

Positional:
- diagram: The TensorDiagram to plot

Keyword:
- scale: Overall scale factor for the diagram (default: 1.0)
- fig_width: Figure width in pixels (default: 1200 * scale)
- fig_height: Figure height in pixels (default: 1200 * scale)
- padding_fraction: Padding as fraction of figure size (default: 0.05)
- aux_fontsize_fraction: Auxiliary text font size as fraction of padding (default: 0.4)
- strokewidth: Width of the marker stroke (default: 2 * scale)
- node_leg_fontsize: Font size for node leg labels (default: scale * 16)
- node_name_fontsize: Font size for node names (default: scale * 24)
- boundary_leg_fontsize: Font size for boundary leg labels (default: scale * 22)
- bond_label_fontsize: Font size for bond labels (default: 16 * scale)
- infinite_edge_linewidth_factor: Multiplier for infinite edge line width (default: 1.5)
- node_repulsion_strength: Strength of node repulsion in layout (default: 3.0)
- dangling_boundary_label_offset: Offset for dangling boundary leg labels from boundary (default: 0.4 * scale)
- show_grids: Whether to show grid lines (default: false)
- show_q_marks: Whether to show question marks on edges (default: false)

Returns:
- A Makie figure containing the plot
"""
function plot_tensor_diagram(diagram::TensorDiagram;
    scale=1.0,
    fig_width=800 * scale,
    fig_height=800 * scale,
    strokewidth=2 * scale,
    node_leg_fontsize=scale * 16,
    node_name_fontsize=scale * 24,
    boundary_leg_fontsize=scale * 22,
    bond_label_fontsize=16 * scale,
    infinite_edge_linewidth_factor=1.5,
    node_repulsion_strength=3.8,
    dangling_boundary_label_offset=0.4 * scale,
    show_grids=false,
    show_q_marks=false)

    # ===== 1. CONSTANTS DERIVED FROM KWARGS =====
    box_width = diagram.width
    box_height = diagram.height

    x_min = -box_width / 2
    x_max = box_width / 2
    y_min = -box_height / 2
    y_max = box_height / 2

    # ===== 2. COMPUTE NODE POSITIONS =====
    node_coordinates = calculate_node_coordinates(diagram; node_repulsion_strength=node_repulsion_strength)

    # ===== 3. COMPUTE NODE LEG POSITIONS =====
    node_and_leg_to_coordinates = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}() # (index of the node in the list of nodes, contraction pattern index) => (x,y)
    node_info = [] # (index of the node in the list of nodes, node, x, y, width, height)

    for (i, (node, node_coordinates)) in enumerate(zip(diagram.nodes, node_coordinates))
        # Calculate node dimensions
        width = default_node_width(node)
        aspect_ratio = node_aspect_ratio(node)
        height = width * aspect_ratio
        x, y = node_coordinates

        # Store node info for drawing
        push!(node_info, (i, node, x, y, width, height))

        # Calculate and store leg positions
        contraction_pattern = diagram.contraction_pattern[i]
        pattern_idx_to_coordinates = get_leg_idx_to_coordinates_dict(node, x, y, width, height, contraction_pattern)
        for (leg_idx, coordinates) in pattern_idx_to_coordinates
            node_and_leg_to_coordinates[(i, leg_idx)] = coordinates
        end
    end

    # ===== 4. COMPUTE BOUNDARY LEG POSITIONS =====
    boundary_leg_coordinates = Dict{Int,Tuple{Float64,Float64}}()

    boundary_configs = [
        ("left", x_min, y_min, box_height, "y"),
        ("right", x_max, y_min, box_height, "y"),
        ("top", y_max, x_min, box_width, "x"),
        ("bottom", y_min, x_min, box_width, "x")
    ]

    for (side, fixed_coord, vary_min, vary_range, vary_axis) in boundary_configs
        legs = get(diagram.boundary_legs, side, Int[])
        position_indices = get(diagram.boundary_legs_posidx, side, Int[])
        n_total = get(diagram.boundary_slots_num, side, 0)

        for (leg_idx, posidx) in enumerate(position_indices)
            varying_coord = vary_min + (posidx / (n_total + 1)) * vary_range
            coord = vary_axis == "x" ? (varying_coord, fixed_coord) : (fixed_coord, varying_coord)
            boundary_leg_coordinates[legs[leg_idx]] = coord
        end
    end

    # Build index mapping for connections
    pattern_index_to_nodelegs = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (node_i, pattern) in enumerate(diagram.contraction_pattern)
        for (leg_idx, pattern_idx) in enumerate(pattern)
            if !haskey(pattern_index_to_nodelegs, pattern_idx)
                pattern_index_to_nodelegs[pattern_idx] = []
            end
            push!(pattern_index_to_nodelegs[pattern_idx], (node_i, leg_idx))
        end
    end

    # ===== 5. DRAWING =====
    fig = Figure(size=(fig_width, fig_height), figure_padding=0) # TODO: adjust the width and height of the figure to the aspect ratio of the box

    # Create layout
    ax = Axis(fig[1, 1], aspect=DataAspect(),
        xticks=MultiplesTicks(5, 1, ""),
        yticks=MultiplesTicks(5, 1, ""),
        xminorticks=IntervalsBetween(5),
        yminorticks=IntervalsBetween(5),
        xminorgridvisible=true,
        yminorgridvisible=true)

    if !show_grids
        hidedecorations!(ax)
        hidespines!(ax)
    end

    # Draw boundary box
    lines!(ax, [x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min],
        color=:black, linewidth=strokewidth * 1.5)

    # Draw contraction edges (bottom layer)
    for (idx, node_list) in pattern_index_to_nodelegs
        label_text = get(diagram.labels, idx, "?")
        linewidth = label_text in INFINITE_LABELS ? strokewidth * infinite_edge_linewidth_factor : strokewidth
        line_color = label_text in INFINITE_LABELS ? :red : (label_text in FINITE_LABELS ? :blue : :black)

        if length(node_list) == 2
            # Internal contraction between two nodes
            (n1, pos1), (n2, pos2) = node_list
            pattern1 = diagram.contraction_pattern[n1]
            pattern2 = diagram.contraction_pattern[n2]

            if haskey(node_and_leg_to_coordinates, (n1, pos1)) && haskey(node_and_leg_to_coordinates, (n2, pos2))
                pos1_coords = node_and_leg_to_coordinates[(n1, pos1)]
                pos2_coords = node_and_leg_to_coordinates[(n2, pos2)]

                lines!(ax, [pos1_coords[1], pos2_coords[1]], [pos1_coords[2], pos2_coords[2]],
                    linewidth=linewidth, color=line_color)

                if label_text != "?" || show_q_marks
                    mid_x = (pos1_coords[1] + pos2_coords[1]) / 2
                    mid_y = (pos1_coords[2] + pos2_coords[2]) / 2
                    text!(ax, mid_x, mid_y, text="$idx:$label_text", align=(:center, :center),
                        fontsize=bond_label_fontsize, color=:red)
                end
            end
        elseif length(node_list) == 1
            # External leg to boundary
            (node_i, pos) = node_list[1]
            pattern = diagram.contraction_pattern[node_i]
            idx_value = pattern[pos]

            if haskey(boundary_leg_coordinates, idx_value) && haskey(node_and_leg_to_coordinates, (node_i, pos))
                node_pos = node_and_leg_to_coordinates[(node_i, pos)]
                boundary_pos = boundary_leg_coordinates[idx_value]

                lines!(ax, [node_pos[1], boundary_pos[1]], [node_pos[2], boundary_pos[2]],
                    linewidth=linewidth, color=line_color)

                if label_text != "?" || show_q_marks
                    mid_x = (node_pos[1] + boundary_pos[1]) / 2
                    mid_y = (node_pos[2] + boundary_pos[2]) / 2
                    text!(ax, mid_x, mid_y, text="$idx_value:$label_text", align=(:center, :center),
                        fontsize=bond_label_fontsize, color=:red)
                end
            end
        end
    end

    # Draw labels for unconnected boundary legs
    for (leg_idx, coord) in boundary_leg_coordinates
        if !haskey(pattern_index_to_nodelegs, leg_idx)
            label_text = get(diagram.labels, leg_idx, "?")
            if label_text == "?" && !show_q_marks
                continue
            end
            x, y = coord

            # Determine offset and alignment based on which boundary we are on
            tol = 1e-5
            text_align = (:center, :center)
            text_pos_x, text_pos_y = x, y

            if abs(x - x_min) < tol # Left boundary
                text_pos_x += dangling_boundary_label_offset
                text_align = (:left, :center)
            elseif abs(x - x_max) < tol # Right boundary
                text_pos_x -= dangling_boundary_label_offset
                text_align = (:right, :center)
            elseif abs(y - y_min) < tol # Bottom boundary
                text_pos_y += dangling_boundary_label_offset
                text_align = (:center, :bottom)
            elseif abs(y - y_max) < tol # Top boundary
                text_pos_y -= dangling_boundary_label_offset
                text_align = (:center, :top)
            end

            text!(ax, text_pos_x, text_pos_y, text="$label_text", align=text_align,
                fontsize=bond_label_fontsize, color=:red)
        end
    end

    # Draw tensor nodes (middle layer)
    for (i, node, x, y, width, height) in node_info
        draw_tensor_node!(ax, node, x, y, width;
            leg_fontsize=node_leg_fontsize, name_fontsize=node_name_fontsize, strokewidth=strokewidth)
    end


    n_left = get(diagram.boundary_slots_num, "left", 0)
    n_right = get(diagram.boundary_slots_num, "right", 0)
    n_top = get(diagram.boundary_slots_num, "top", 0)
    n_bottom = get(diagram.boundary_slots_num, "bottom", 0)

    # Draw boundary legs (top layer)
    boundary_draw_configs = [
        ("left", (x_min, y_min), (x_min, y_max), n_left),
        ("right", (x_max, y_min), (x_max, y_max), n_right),
        ("top", (x_min, y_max), (x_max, y_max), n_top),
        ("bottom", (x_min, y_min), (x_max, y_min), n_bottom)
    ]

    for (side, start_pos, end_pos, n_total) in boundary_draw_configs
        legs = get(diagram.boundary_legs, side, Int[])
        positions = get(diagram.boundary_legs_posidx, side, Int[])
        activity = [i in positions for i in 1:n_total]
        perm = sortperm(positions)
        draw_legs_on_segment!(ax, start_pos, end_pos, n_total, legs, activity; perm=invperm(perm),
            fontsize=boundary_leg_fontsize, strokewidth=strokewidth)
    end

    return fig
end

function Base.show(io::IO, m::Union{MIME"image/png",MIME"image/svg+xml"}, node::TensorDiagram)
    fig = plot_tensor_diagram(node)
    show(io, m, fig)
end

#################################################
# topic 4: DISPLAY UTILITIES
#################################################

"""
    display_scrollable_figures(figure_lists...; max_height="500px")

Display lists of figures side-by-side in a single scrollable container using a grid layout.
This is useful for comparing sequences of diagrams or plots in Jupyter notebooks.

# Arguments
- `figure_lists...`: Variable number of lists (Vectors) containing figures to be displayed. All lists must have the same length. 
- `max_height`: String specifying the maximum height of the scrollable area (default: "500px").

# Details
- All input lists must be of equal length. If you have unequal lists, pad the shorter ones with `missing`.
- The display uses an HTML grid with `overflow: auto`.
- Items that are `AbstractString` are wrapped in `<pre>` tags.
- Other items are displayed using `show(io, MIME"image/svg+xml"(), item)`.
- `missing` items result in an empty cell.

# Example
```julia
list1 = [plot1, plot2]
list2 = [plot3, plot4]
display_scrollable_figures(list1, list2; max_height="600px")
```
"""
function display_scrollable_figures(figure_lists...; max_height="500px")
    if !all(l -> length(l) == length(figure_lists[1]), figure_lists)
        error("All figure lists must have the same length. Note that you can pad shorter lists with `missing` values.")
    end

    io = IOBuffer()

    # Grid-based layout for consistent column widths and common scrolling
    n_cols = length(figure_lists)
    # Using max-content ensures columns are wide enough for the content (preventing horizontal chop).
    # The container has overflow: auto, so it will scroll if total width exceeds viewport.
    grid_style = "display: grid; grid-template-columns: repeat($(n_cols), max-content); gap: 20px; max-height: $(max_height); overflow: auto; padding: 10px;"

    print(io, """<div style="$(grid_style)">""")

    n_items = isempty(figure_lists) ? 0 : length(figure_lists[1])

    for i in 1:n_items
        for list in figure_lists
            item = list[i]
            # Cell style: border-bottom provides the row separator feel
            print(io, """<div style="padding: 10px; border-bottom: 1px solid #ccc;">""")
            if !ismissing(item)
                if item isa AbstractString
                    print(io, "<pre>", item, "</pre>")
                else
                    show(io, MIME"image/svg+xml"(), item)
                end
            end
            print(io, "</div>")
        end
    end
    print(io, "</div>")

    display(MIME("text/html"), String(take!(io)))
end
