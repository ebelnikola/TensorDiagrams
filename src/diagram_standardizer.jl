
"""
    standardize_diagram(diagram::TensorDiagram)

Renames the indices of a TensorDiagram based on its topology to a canonical form.
1. Boundary indices are renamed to -1, -2, ... based on their position (Left: Bottom->Top, Top: Left->Right, Right: Bottom->Top, Bottom: Left->Right).
2. Internal indices are renamed to 1, 2, ... based on a traversal starting from boundary indices.

Returns a tuple `(new_diag::TensorDiagram, success::Bool)`. `success` is false if there are internal indices unreachable from the boundary.
"""
function standardize_diagram(diagram::TensorDiagram)
    # Create a copy to avoid modifying the original diagram
    new_diag = copy(diagram)

    # 1. Standardize Boundary Indices
    # Canonical order of sides: Left, Top, Right, Bottom
    sides_order = ["left", "top", "right", "bottom"]

    # Collect all boundary legs with their info: (side_index, pos_idx, old_index)
    boundary_legs_info = Tuple{Int,Int,Int}[]

    for (side_i, side) in enumerate(sides_order)
        if haskey(new_diag.boundary_legs, side)
            pattern_indices = new_diag.boundary_legs[side]
            pos_indices = new_diag.boundary_legs_posidx[side]
            for (i, idx) in enumerate(pattern_indices)
                push!(boundary_legs_info, (side_i, pos_indices[i], idx))
            end
        end
    end

    # Sort boundary legs: primary key side_i, secondary key pos_idx
    sort!(boundary_legs_info, by=x -> (x[1], x[2]))

    # Create mapping for boundary indices
    canonical_boundary_order = [x[3] for x in boundary_legs_info]
    index_map = Dict(old => -i for (i, old) in enumerate(canonical_boundary_order))

    # 2. Standardize Internal Indices
    # Build Adjacency Map: index -> vector of node indicies for nodes connected by that index
    index_to_nodes = Dict{Int,Vector{Int}}()
    for (node_idx, pattern) in enumerate(new_diag.contraction_pattern)
        for idx in pattern
            if !haskey(index_to_nodes, idx)
                index_to_nodes[idx] = Int[]
            end
            push!(index_to_nodes[idx], node_idx)
        end
    end

    # Traversal to determine internal index order
    current_internal_idx = 1
    node_visited = Set{Int}()
    index_processed = Set{Int}()

    # Queue for BFS-like traversal: stores node indices
    node_queue = Int[]

    # Start traversal from boundary legs in canonical order
    for start_boundary_idx in canonical_boundary_order
        if !haskey(index_to_nodes, start_boundary_idx)
            continue # Dangling index (no tensor attached)
        end

        # Add attached nodes to queue
        let node_idx = index_to_nodes[start_boundary_idx][1]
            if !(node_idx in node_visited)
                push!(node_queue, node_idx)
                push!(node_visited, node_idx)
            end
        end

        # Process queue
        while !isempty(node_queue)
            node_idx = popfirst!(node_queue)

            # Iterate through legs of the current tensor in order

            pattern = new_diag.contraction_pattern[node_idx]
            for idx in pattern
                # We are interested in positive (internal) indices that we haven't renamed yet
                if idx > 0 && !(idx in index_processed)
                    # Assign new name
                    index_map[idx] = current_internal_idx
                    current_internal_idx += 1
                    push!(index_processed, idx)

                    # Add neighbors via this index to queue
                    for v in index_to_nodes[idx]
                        if !(v in node_visited)
                            push!(node_queue, v)
                            push!(node_visited, v)
                        end
                    end
                end
            end

        end
    end

    # Handle any remaining connected components (if any disconnected from boundary)

    all_indices = Set(Iterators.flatten(new_diag.contraction_pattern))
    unreachable_indices = filter(x -> x > 0 && !haskey(index_map, x), all_indices) |> collect

    success = true
    if !isempty(unreachable_indices)
        success = false
        sort!(unreachable_indices)
        for idx in unreachable_indices
            index_map[idx] = current_internal_idx
            current_internal_idx += 1
        end
    end

    # 3. Apply Renaming
    relabel!(new_diag, index_map)


    # 4. The last touch, we order the boundary legs in their storage such that posidx would be ordered
    for (side, legs) in new_diag.boundary_legs
        perm = sortperm(new_diag.boundary_legs_posidx[side])
        new_diag.boundary_legs[side] = legs[perm]
        new_diag.boundary_legs_posidx[side] = new_diag.boundary_legs_posidx[side][perm]
    end

    return new_diag, success
end
