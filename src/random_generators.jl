#################################################
# TOC
# topic 0:           EXPORTS AND IMPORTS
# topic 1:           NODE GENERATION
# topic 2:           DIAGRAM GENERATION
#################################################

############################################
# topic 0: EXPORTS AND IMPORTS
############################################

export generate_random_tensor_node, add_random_forbidden_label_combinations!, add_random_allowed_label_combinations!, generate_random_tensor_diagram, add_random_dangling_indices!

using Random

############################################
# topic 1: NODE GENERATION
############################################

"""
    generate_random_tensor_node(num_legs::Int, name::String; min_num_allowed_rnd_labels=1, max_num_of_allowed_rnd_labels=4, necessary_labels=nothing)

Generates a random TensorNode with a specified number of legs and name.
Legs are randomly distributed among sides.
Allowed labels for each leg are a random subset of all possible labels plus any necessary labels provided.
Allowed/Forbidden label combinations are initialized as empty.

# Arguments
- `num_legs::Int`: Total number of legs.
- `name::String`: Name of the node.
- `min_num_allowed_rnd_labels::Int`: Minimum number of random allowed labels to pick per leg. Default is 1.
- `max_num_of_allowed_rnd_labels::Int`: Maximum number of random allowed labels to pick per leg. Default is 4.
- `necessary_labels`: Optional collection of labels (Vector of strings). If provided, these labels are added to `allowed_labels` for every leg. 

# Returns
- A `TensorNode` object.
"""
function generate_random_tensor_node(num_legs::Int, name::String; min_num_allowed_rnd_labels=1, max_num_of_allowed_rnd_labels=4, necessary_labels=nothing)

    sides = ["left", "right", "top", "bottom"]
    legs = Dict{String,Array{Int,1}}()

    # 1. Randomly distribute legs
    for i in 1:num_legs
        side = rand(sides)
        if !haskey(legs, side)
            legs[side] = Int[]
        end
        push!(legs[side], i)
    end
    legs = Dict(side => shuffle(legs[side]) for side in keys(legs))

    if max_num_of_allowed_rnd_labels == 0 && isnothing(necessary_labels)
        return TensorNode(
            name=name,
            color=set_alpha(rand(PALETTE), 0.6),
            legs=legs,
        )
    end


    # 2. Determine allowed labels for each leg
    all_labels = collect(INFINITE_LABELS ∪ FINITE_LABELS)

    allowed_labels = [shuffle(all_labels)[1:rand(min_num_allowed_rnd_labels:min(max_num_of_allowed_rnd_labels, length(all_labels)))] for _ in 1:num_legs]

    if !isnothing(necessary_labels)
        unique!.(append!.(allowed_labels, [necessary_labels]))
    end


    return TensorNode(
        name=name,
        color=set_alpha(rand(PALETTE), 0.6),
        legs=legs,
        allowed_labels=allowed_labels,
    )
end

function set_alpha(color::RGBA, alpha)
    return RGBA(color.r, color.g, color.b, alpha)
end

"""
    add_random_forbidden_label_combinations!(tensor_node::TensorNode; max_num_forbidden_combinations=2)

Generates random label combinations from the node's `allowed_labels` and appends them to the node's `forbidden_label_combinations`.

# Arguments
- `tensor_node::TensorNode`: The tensor node to modify in-place.
- `max_num_forbidden_combinations::Int`: The maximum number of random forbidden combinations to generate. Default is 2.
# Details
- Combinations are formed by picking one random label from `allowed_labels` for each leg.
- Duplicate combinations are removed using `unique`.
- Appends to any existing `forbidden_label_combinations` in the node.
"""
function add_random_forbidden_label_combinations!(tensor_node::TensorNode; max_num_forbidden_combinations=2)
    allowed_labels = tensor_node.allowed_labels


    forbidden_label_combinations = unique([[rand(allowed_labels[i]) for i in 1:length(allowed_labels)] for _ in 1:max_num_forbidden_combinations])
    append!(tensor_node.forbidden_label_combinations, forbidden_label_combinations)
    tensor_node.forbidden_label_combinations = unique(tensor_node.forbidden_label_combinations)
    return nothing
end

"""
    add_random_allowed_label_combinations!(tensor_node::TensorNode; max_num_allowed_combinations=2)

Generates random label combinations from the node's `allowed_labels` and appends them to the node's `allowed_label_combinations`.

# Arguments
- `tensor_node::TensorNode`: The tensor node to modify in-place.
- `max_num_allowed_combinations::Int`: The maximum number of random allowed combinations to generate. Default is 2.
# Details
- Combinations are formed by picking one random label from `allowed_labels` for each leg.
- Duplicate combinations are removed using `unique`.
- Appends to any existing `allowed_label_combinations` in the node.
"""
function add_random_allowed_label_combinations!(tensor_node::TensorNode; max_num_allowed_combinations=2)
    allowed_labels = tensor_node.allowed_labels
    allowed_label_combinations = unique([[rand(allowed_labels[i]) for i in 1:length(allowed_labels)] for _ in 1:max_num_allowed_combinations])
    append!(tensor_node.allowed_label_combinations, allowed_label_combinations)
    tensor_node.allowed_label_combinations = unique(tensor_node.allowed_label_combinations)
    return nothing
end

############################################
# topic 2: DIAGRAM GENERATION
############################################


"""
    generate_random_tensor_diagram(boundary_legs_num::Int; max_nodes_num::Int=5, max_leg_num::Int=5, max_free_slots_num_per_side::Int=2)

Generates a random TensorDiagram with a specified total number of boundary legs.
Follows a specific construction strategy to ensure connectivity and valid boundary placement.

# Arguments
- `boundary_legs_num::Int`: Total number of negative indices (boundary legs) to place.
- `max_nodes_num::Int`: Maximum number of nodes to generate (randomly picks 1 to max). Default is 5.
- `max_leg_num::Int`: Maximum number of legs per node. Default is 5.
- `max_free_slots_num_per_side::Int`: Maximum number of free slots to leave unconnected per boundary side. Default is 2.

# Construction Strategy
1. Determine `nodes_num` randomly from 1 to `max_nodes_num`.
2. Generate node leg counts such that total legs > `boundary_legs_num` and (total - boundary) is even.
3. Distribute negative indices (-1, -2, ...) into the node leg slots, preferring nodes with multiple available slots to encourage connectivity.
4. Connect remaining slots with unique internal indices (positive integers), minimizing self-loops by preferring to connect different nodes.
5. Create `TensorNode` objects for each node:
    - Adds "z" as a necessary label to all legs.
    - Adds an "all z" combination to `allowed_label_combinations` to ensure at least one consistent state.
6. Randomly distributes the boundary legs (-1, etc.) to diagram boundaries (left, right, top, bottom).
7. Returns a consistent `TensorDiagram`.

# Returns
- A `TensorDiagram` object.
"""
function generate_random_tensor_diagram(boundary_legs_num::Int; max_nodes_num::Int=5, max_leg_num::Int=5, max_free_slots_num_per_side::Int=2)

    # 1. Pick number of nodes
    # Ensure it is theoretically possible to satisfy constraints:
    # We need total_slots >= boundary_legs_num.
    min_nodes_required = max(ceil(Int, boundary_legs_num / max_leg_num), 1)
    if min_nodes_required > max_nodes_num
        error("Impossible constraints: max_nodes_num ($max_nodes_num) * max_leg_num ($max_leg_num) < boundary_legs_num ($boundary_legs_num)")
    end

    nodes_num = rand(min_nodes_required:max_nodes_num)


    # 2. Generate contraction pattern structure (leg counts)
    # We ensure that the difference between total slots and boundary legs is even
    # because internal bonds consume legs in pairs (2 legs per bond).

    node_leg_counts = [rand(1:max_leg_num) for _ in 1:nodes_num]
    total_slots = sum(node_leg_counts)
    if (total_slots - boundary_legs_num) % 2 != 0
        candidates = findall(x -> x < max_leg_num, node_leg_counts)

        if isempty(candidates)
            if max_leg_num == 1
                if nodes_num == max_nodes_num
                    node_leg_counts = node_leg_counts[1:end-1]
                    nodes_num -= 1
                else
                    push!(node_leg_counts, 1)
                    nodes_num += 1
                end
            else
                node_leg_counts[1] -= 1
            end
        else
            idx = rand(candidates)
            node_leg_counts[idx] += 1
        end
    end


    # Increment leg counts if we don't have enough slots
    while sum(node_leg_counts) < boundary_legs_num
        # Find nodes that can accept more legs
        candidates = findall(x -> x < max_leg_num, node_leg_counts)
        # Increment pair of legs
        idx1 = rand(candidates)
        node_leg_counts[idx1] += 1

        candidates = findall(x -> x < max_leg_num, node_leg_counts)

        idx2 = rand(candidates)
        for _ = 1:length(candidates)*2
            if idx2 != idx1
                break
            end
            idx2 = rand(candidates)
        end

        node_leg_counts[idx2] += 1
    end
    total_slots = sum(node_leg_counts)

    # Create empty slots for each node
    # Each slot is represented by a placeholder (0) initially
    contraction_pattern = [zeros(Int, count) for count in node_leg_counts]

    # Collect all available slot coordinates: (node_index, leg_index)
    available_slots = Tuple{Int,Int}[]
    for ni in 1:nodes_num
        for li in 1:node_leg_counts[ni]
            push!(available_slots, (ni, li))
        end
    end
    shuffle!(available_slots)

    # 3. Distribute numbers -1, -2, ... -boundary_legs_num
    # We assign them randomly, preferring slots from nodes that have at least one other slot available to allow for internal connections later.
    assigned_boundary_indices = Int[]
    for i in 1:boundary_legs_num
        # Prefer slots from nodes that have at least one other slot available (to allow for internal connection later)
        candidate_idx = findfirst(slot -> count(s -> s[1] == slot[1], available_slots) > 1, available_slots)

        if !isnothing(candidate_idx)
            (ni, li) = available_slots[candidate_idx]
            deleteat!(available_slots, candidate_idx)
        else
            (ni, li) = pop!(available_slots)
        end

        idx = -i
        contraction_pattern[ni][li] = idx
        push!(assigned_boundary_indices, idx)
    end

    # 4. Fill internal bonds
    # Connect pairs of slots with unique internal indices, avoiding self-loops if possible.
    next_internal_idx = 1
    while !isempty(available_slots)
        (n1, l1) = pop!(available_slots)

        # Try to find a slot from a different node to avoid self-contraction
        idx_different_node = findfirst(slot -> slot[1] != n1, available_slots)

        if isnothing(idx_different_node)
            (n2, l2) = pop!(available_slots)
        else
            (n2, l2) = available_slots[idx_different_node]
            deleteat!(available_slots, idx_different_node)
        end

        contraction_pattern[n1][l1] = next_internal_idx
        contraction_pattern[n2][l2] = next_internal_idx
        next_internal_idx += 1
    end

    # 5. Generate Random Nodes
    nodes = TensorNode[]
    processed_indices = Set{Int}()
    for i in 1:nodes_num
        # Create node with correct number of legs
        nl = node_leg_counts[i]

        # Pick the first finite label as the consistent one (e.g. "z")
        consistent_lbl = first(sort(collect(FINITE_LABELS))) # deterministic choice

        if rand() > 0.9
            node = generate_random_tensor_node(nl, "n$(i)"; max_num_of_allowed_rnd_labels=0)
        else
            node = generate_random_tensor_node(nl, "n$(i)"; necessary_labels=[consistent_lbl])
        end


        for (leg_idx, global_idx) in enumerate(contraction_pattern[i])
            if global_idx > 0 # internal leg
                if !(global_idx in processed_indices)
                    node.dual[leg_idx] = true
                    push!(processed_indices, global_idx)
                end
            end
        end

        push!(nodes, node)
    end

    # 6. Diagram Boundary Setup
    # Start with empty boundaries
    boundary_legs = Dict("left" => Int[], "right" => Int[], "top" => Int[], "bottom" => Int[])
    boundary_legs_posidx = Dict("left" => Int[], "right" => Int[], "top" => Int[], "bottom" => Int[])

    sides = ["left", "right", "top", "bottom"]
    shuffled_b_inds = shuffle(assigned_boundary_indices)

    # Assign each leg to a random side
    for idx in shuffled_b_inds
        side = rand(sides)
        push!(boundary_legs[side], idx)
    end

    # Finalize structure
    boundary_slots_num_dict = Dict{String,Int}()
    for side in sides
        boundary_slots_num_dict[side] = (get(boundary_legs, side, Int[]) |> length) + rand(0:max_free_slots_num_per_side)
    end

    # boundary_legs = boundary_legs |> Dict{String,Vector{Int}} what the hell was this line? 
    boundary_legs_posidx = Dict{String,Vector{Int}}()

    for side in sides
        n_legs = length(boundary_legs[side])
        n_pos = boundary_slots_num_dict[side]
        boundary_legs_posidx[side] = shuffle(collect(1:n_pos))[1:n_legs]
    end

    labels = Dict{Int,String}()



    return TensorDiagram(
        nodes=nodes,
        contraction_pattern=contraction_pattern,
        boundary_slots_num=boundary_slots_num_dict,
        boundary_legs=boundary_legs,
        boundary_legs_posidx=boundary_legs_posidx,
        labels=labels
    )
end;


"""
    add_random_dangling_indices!(diag::TensorDiagram; max_num_dangling::Int=3)

Adds random dangling (unconnected) boundary indices to a diagram in place.
Dangling indices are boundary legs that are not connected to any tensor node.

# Arguments
- `diag::TensorDiagram`: The diagram to modify in place.
- `max_num_dangling::Int`: Maximum number of dangling indices to add. Default is 3.

# Details
- Only adds indices to free slots (boundary positions not already occupied).
- Generates new negative indices that don't conflict with existing ones.
- Randomly distributes the new indices across available free slots on all sides.

# Returns
- `true` if at least one dangling index was added, `false` if no free slots were available or `max_num_dangling < 1`.
"""
function add_random_dangling_indices!(diag::TensorDiagram; max_num_dangling::Int=3)

    boundary_legs = vcat(values(diag.boundary_legs)...)
    new_neg_index_cnt = [min(minimum(boundary_legs; init=0), 0) - 1]
    nni() = generate_nfi!(new_neg_index_cnt; dir="negative")


    free_slots_total = 0
    for side in ["left", "right", "top", "bottom"]
        slots = get(diag.boundary_slots_num, side, 0)
        occupied = length(get(diag.boundary_legs, side, Int[]))
        free_slots_total += (slots - occupied)
    end

    if free_slots_total < 1 || max_num_dangling < 1
        return false
    end

    dangling_num = rand(1:min(max_num_dangling, free_slots_total))

    free_slots = Tuple{String,Int}[]
    for side in keys(diag.boundary_legs)
        positions = diag.boundary_legs_posidx[side]
        total_slots = diag.boundary_slots_num[side]

        free_slots_on_side = setdiff(1:total_slots, positions)

        append!(free_slots, [(side, slot) for slot in free_slots_on_side])
    end

    shuffle!(free_slots)
    for i in 1:dangling_num
        (side, posidx) = free_slots[i]
        new_idx = nni()
        push!(diag.boundary_legs[side], new_idx)
        push!(diag.boundary_legs_posidx[side], posidx)
    end

    return true
end
