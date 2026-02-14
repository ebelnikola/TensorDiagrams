"""
    standardize_diagram(diagram::TensorDiagram)

Renames the indices of a TensorDiagram based on its topology to a canonical form.
1. Boundary indices are renamed to -1, -2, ... based on their position (Left boundary: Bottom->Top, Top boundary: Left->Right, Right boundary: Bottom->Top, Bottom boundary: Left->Right).
2. Internal indices are renamed to 1, 2, ... based on a traversal starting from the new -1 boundary index.
3. Boundary indices are stored in the order of their position indices. 
4. Missing labels are populated with "?".
5. Nodes and contraction patterns are sorted based on the contraction pattern.

Returns a tuple `(new_diag::TensorDiagram, success::Bool)`. `success` is false if there are internal indices unreachable from the boundary.
"""
function standardize_diagram(diagram::TensorDiagram)
    # Create a copy to avoid modifying the original diagram
    standard_diag = copy(diagram)

    # 1. Pattern Index Update Map for Boundary Legs
    # Canonical order of sides: Left, Top, Right, Bottom
    sides_order = ["left", "top", "right", "bottom"]

    # Collect all boundary legs with their info: (side index, pos idx, old pattern index)
    boundary_legs_info = Tuple{Int,Int,Int}[]

    for (side_i, side) in enumerate(sides_order)
        if haskey(standard_diag.boundary_legs, side)
            pattern_indices = standard_diag.boundary_legs[side]
            pos_indices = standard_diag.boundary_legs_posidx[side]
            for (pos_idx, pat_idx) in zip(pos_indices, pattern_indices)
                push!(boundary_legs_info, (side_i, pos_idx, pat_idx))
            end
        end
    end

    # Sort boundary legs: primary key side_i, secondary key pos_idx
    sort!(boundary_legs_info, by=x -> (x[1], x[2]))

    # Create mapping for boundary indices
    canonical_boundary_order = [x[3] for x in boundary_legs_info]
    pat_idx_update_map = Dict(old => -i for (i, old) in enumerate(canonical_boundary_order))

    # 2. Pattern Index Update Map for Internal Legs
    # Build Adjacency Map: index -> vector of node indicies for nodes connected by that index
    index_to_nodes = Dict{Int,Vector{Int}}()
    for (node_idx, pattern) in enumerate(standard_diag.contraction_pattern)
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
        let node_idx = index_to_nodes[start_boundary_idx][1] # boundary leg should have only one attached node
            if !(node_idx in node_visited)
                push!(node_queue, node_idx)
                push!(node_visited, node_idx)
            end
        end

        # Process queue
        while !isempty(node_queue)
            node_idx = popfirst!(node_queue)

            # Iterate through legs of the current tensor in order

            pattern = standard_diag.contraction_pattern[node_idx]
            for idx in pattern
                # We are interested in positive (internal) indices that we haven't renamed yet
                if idx > 0 && !(idx in index_processed)
                    # Assign new name
                    pat_idx_update_map[idx] = current_internal_idx
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

    all_indices = Set(Iterators.flatten(standard_diag.contraction_pattern))
    unreachable_indices = filter(x -> x > 0 && !haskey(pat_idx_update_map, x), all_indices) |> collect

    success = true
    if !isempty(unreachable_indices)
        success = false
        sort!(unreachable_indices)
        for idx in unreachable_indices
            pat_idx_update_map[idx] = current_internal_idx
            current_internal_idx += 1
        end
    end

    # 3. Apply Renaming
    relabel!(standard_diag, pat_idx_update_map)


    # 4. We order the boundary legs in their storage such that posidx would be ordered
    for (side, legs) in standard_diag.boundary_legs
        perm = sortperm(standard_diag.boundary_legs_posidx[side])
        standard_diag.boundary_legs[side] = legs[perm]
        standard_diag.boundary_legs_posidx[side] = standard_diag.boundary_legs_posidx[side][perm]
    end

    # 5. Populate labels with question marks for missing keys
    all_indices = Set(vcat(standard_diag.contraction_pattern...))
    union!(all_indices, vcat(values(standard_diag.boundary_legs)...))

    for idx in all_indices
        if !haskey(standard_diag.labels, idx)
            standard_diag.labels[idx] = "?"
        end
    end

    # 6. Sort nodes and contraction patterns
    perm = sortperm(standard_diag.contraction_pattern)
    standard_diag.nodes = standard_diag.nodes[perm]
    standard_diag.contraction_pattern = standard_diag.contraction_pattern[perm]

    return standard_diag, success
end
