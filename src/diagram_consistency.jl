
"""
    check_node_consistency(node::TensorNode)

Performs a series of consistency checks on a `TensorNode` and returns a dictionary of results.

# Checks performed:
- `check1`: Legs indices are consecutive integers (no gaps).
- `check2`: Number of sets of allowed labels matches the total number of legs.
- `check2'`: All allowed label sets are non-empty.
- `check3`: All allowed labels are valid (belong to `INFINITE_LABELS` or `FINITE_LABELS`).
- `check4`: All allowed label combinations have the same length as the number of legs.
- `check5`: All forbidden label combinations have the same length as the number of legs.
- `check6`: All labels used in allowed combinations are present in the node's `allowed_labels`.
- `check7`: All labels used in forbidden combinations are present in the node's `allowed_labels`.
- `check8`: The node does not define both allowed and forbidden label combinations simultaneously.
- `check9`: The `dual` vector has the same length as the total number of legs.
- `check10`: The `flipped` vector has the same length as the total number of legs.
- `check11`: `out_legs_num` is between 0 and the total number of legs.

# Returns
- A `Dict{String, Bool}` containing the result of each check.
"""
function check_node_consistency(node::TensorNode)
    results = Dict{String,Bool}()

    all_legs = sort(vcat(values(node.legs)...))
    results["check1"] = all((all_legs[2:end] - all_legs[1:end-1]) .== 1)

    num_of_legs = length(all_legs)

    results["check2"] = length(node.allowed_labels) == num_of_legs

    results["check2'"] = all([!isempty(al) for al in node.allowed_labels])

    allowed_labels_set = Set{String}(vcat(node.allowed_labels...))

    results["check3"] = allowed_labels_set ⊆ (INFINITE_LABELS ∪ FINITE_LABELS)

    results["check4"] = all(length(combination) == num_of_legs for combination in node.allowed_label_combinations)

    results["check5"] = all(length(combination) == num_of_legs for combination in node.forbidden_label_combinations)

    results["check6"] = Set(Iterators.flatten(node.allowed_label_combinations)) ⊆ allowed_labels_set

    results["check7"] = Set(Iterators.flatten(node.forbidden_label_combinations)) ⊆ allowed_labels_set

    results["check8"] = !(length(node.allowed_label_combinations) > 0 && length(node.forbidden_label_combinations) > 0)

    results["check9"] = length(node.dual) == num_of_legs

    results["check10"] = length(node.flipped) == num_of_legs

    results["check11"] = 0 <= node.out_legs_num <= num_of_legs

    return results
end


"""
    is_dual(node::TensorNode, leg_idx::Int)

Determines if a specific leg of a `TensorNode` is dual.
Condition: (leg_idx > out_legs_num) XOR dual[leg_idx] XOR flipped[leg_idx]
"""
function is_dual(node::TensorNode, leg_idx::Int)
    cond1 = leg_idx > node.out_legs_num
    cond2 = node.dual[leg_idx]
    cond3 = node.flipped[leg_idx]
    # odd number of conditions true => dual
    return (cond1 + cond2 + cond3) % 2 == 1
end

"""
    check_diagram_consistency(diag::TensorDiagram)

Performs a series of consistency checks on a `TensorDiagram` and returns a dictionary of results.

# Checks performed:
- `check0`: All nodes in the diagram are consistent (pass `check_node_consistency`).
- `check1`: Number of nodes matches number of contraction patterns.
- `check2`: Number of legs in each node matches the length of its corresponding contraction pattern.
- `check3`: All internal indices (positive) appear exactly twice.
- `check4`: All external indices (negative) appear exactly once in the contraction pattern.
- `check5`: Index 0 is not used.
- `check6`: All boundary indices are unique.
- `check7`: All negative indices in the contraction pattern are present in boundary legs.
- `check8`: All boundary legs have negative indices.
- `check9`: `boundary_slots_num` contains the "horizontal" key.
- `check10`: `boundary_slots_num` contains the "vertical" key.
- `check11`: Keys of `boundary_legs` are a subset of ["left", "right", "top", "bottom"].
- `check12`: Keys of `boundary_legs_posidx` match keys of `boundary_legs`.
- `check13`: Length of legs matches length of position indices (per side).
- `check14`: Number of legs on a side is less than or equal to the maximum allowed for that orientation (per side).
- `check15`: Minimum position index is >= 1 (per side).
- `check16`: Maximum position index is <= maximum allowed for that orientation (per side).
- `check17`: Position indices are unique (per side).
- `check18`: Keys of `labels` are a subset of all indices (contraction pattern + boundary).
- `check19`: All label values are valid. A label is valid if it is in `INFINITE_LABELS`, `FINITE_LABELS`, or "?". Additionally, a label can be of the form "target←source" where "target" is a valid label (as defined before) and "source" is a valid SpaceID string.
- `check20`: If `node_coordinates` is defined, its length matches the number of nodes.
- `check21`: Labels on legs must be in the `allowed_labels` for that leg (or "?").
- `check22`: If `allowed_label_combinations` is specified, the (partial) labels on the legs must match a subset of at least one allowed combination.
- `check23`: If `forbidden_label_combinations` is specified and all leg labels are known (no "?"), the combination must not be in the forbidden list.
- `check24`: All dangling indices (boundary indices not present in the contraction pattern) must have infinite-dimensional labels (present in `INFINITE_LABELS`).
- `check25`: All internal bonds connect a "dual" leg with a "not dual" leg.

# Returns
- A `Dict{String, Bool}` containing the result of each check.
"""
function check_diagram_consistency(diag::TensorDiagram)
    results = Dict{String,Bool}()

    results["check0"] = all([all(values(check_node_consistency(node))) for node in diag.nodes])
    results["check1"] = length(diag.nodes) == length(diag.contraction_pattern)

    results["check2"] = all([sum(length.(values(x.legs))) == length(y) for (x, y) in zip(diag.nodes, diag.contraction_pattern)])

    contraction_pattern_indices = diag.contraction_pattern |> Iterators.flatten |> collect
    freq_dict = Dict{Int,Int}()
    for idx in contraction_pattern_indices
        if haskey(freq_dict, idx)
            freq_dict[idx] += 1
        else
            freq_dict[idx] = 1
        end
    end
    contraction_pattern_indices_unique = unique(contraction_pattern_indices)

    results["check3"] = all(==(2), [freq_dict[idx] for idx in contraction_pattern_indices_unique if idx > 0])
    results["check4"] = all(==(1), [freq_dict[idx] for idx in contraction_pattern_indices_unique if idx < 0])
    results["check5"] = 0 ∉ contraction_pattern_indices_unique


    boundary_indices = vcat(values(diag.boundary_legs)...)

    results["check6"] = length(boundary_indices) == length(unique(boundary_indices))

    results["check7"] = all(idx in boundary_indices for idx in contraction_pattern_indices_unique if idx < 0)

    results["check8"] = all(boundary_indices .< 0)

    results["check9"] = haskey(diag.boundary_slots_num, "horizontal")
    results["check10"] = haskey(diag.boundary_slots_num, "vertical")
    results["check11"] = keys(diag.boundary_legs) ⊆ ["left", "right", "top", "bottom"]
    results["check12"] = keys(diag.boundary_legs_posidx) == keys(diag.boundary_legs)


    for side = ["left", "right", "top", "bottom"]
        num_key = side ∈ ["left", "right"] ? "horizontal" : "vertical"
        legs = get(diag.boundary_legs, side, Int[])
        posidx = get(diag.boundary_legs_posidx, side, Int[])

        results["check13 $side"] = length(legs) == length(posidx)
        results["check14 $side"] = length(legs) <= diag.boundary_slots_num[num_key]
        if !isempty(posidx)
            results["check15 $side"] = minimum(posidx) >= 1
            results["check16 $side"] = maximum(posidx) <= diag.boundary_slots_num[num_key]
        else
            results["check15 $side"] = true
            results["check16 $side"] = true
        end
        results["check17 $side"] = length(unique(posidx)) == length(posidx)
    end

    results["check18"] = keys(diag.labels) ⊆ Set(vcat(contraction_pattern_indices_unique, boundary_indices))

    results["check19"] = all([_valid_label(lbl) for lbl in values(diag.labels)])


    results["check20"] = isnothing(diag.node_coordinates) || length(diag.node_coordinates) == length(diag.nodes)

    for (i, node) in enumerate(diag.nodes)
        pattern = diag.contraction_pattern[i]

        labels = [get(diag.labels, idx, "?") for idx in pattern]


        results["check21 node $i"] = all([lbl == "?" || lbl in node.allowed_labels[leg] for (leg, lbl) in enumerate(labels)])

        mask = labels .!= "?"

        if length(node.allowed_label_combinations) > 0
            allowed_sub_combinations = node.allowed_label_combinations .|> x -> x[mask]
            results["check22 node $i"] = labels[mask] ∈ allowed_sub_combinations
        else
            results["check22 [x] node $i"] = true
        end

        if length(node.forbidden_label_combinations) > 0 && "?" ∉ labels
            results["check23 node $i"] = labels[mask] ∉ node.forbidden_label_combinations
        else
            results["check23 [x] node $i"] = true
        end

    end

    is_inf(lbl) = split(lbl, "←")[1] in INFINITE_LABELS
    results["check24"] = all(is_inf.([diag.labels[idx] for idx in boundary_indices if idx ∉ contraction_pattern_indices_unique]))

    bond_endpoints = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (node_idx, pattern) in enumerate(diag.contraction_pattern)
        for (leg_idx, bond_idx) in enumerate(pattern)
            if bond_idx > 0
                if !haskey(bond_endpoints, bond_idx)
                    bond_endpoints[bond_idx] = Tuple{Int,Int}[]
                end
                push!(bond_endpoints[bond_idx], (node_idx, leg_idx))
            end
        end
    end

    results["check25"] = true
    for (bond, endpoints) in bond_endpoints
        if length(endpoints) == 2
            (n1, l1) = endpoints[1]
            (n2, l2) = endpoints[2]
            if is_dual(diag.nodes[n1], l1) == is_dual(diag.nodes[n2], l2)
                results["check25"] = false
                break
            end
        end
    end

    return results
end

function _valid_label(str)
    components = split(str, "←")
    lbl = components[1]
    valid_target_space = lbl in INFINITE_LABELS ∪ FINITE_LABELS ∪ Set(["?"])
    valid_domain_space = true
    if length(components) == 2
        domain_space_id = components[2]
        valid_domain_space = is_valid_space_id(domain_space_id)
    end
    if length(components) > 2
        return false
    end
    return valid_target_space && valid_domain_space
end

"""
    print_failed_checks(results::Dict{String,Bool})

Prints the keys of the checks that failed (value is `false`).
If all checks passed, prints "all checks are passed".

# Arguments
- `results`: A dictionary mapping check names to boolean results.
"""
function print_failed_checks(results::Dict{String,Bool})
    failed_checks = [k for (k, v) in results if !v]
    if isempty(failed_checks)
        println("all checks are passed")
    else
        println("FAILED CHECKS:")
        for check in sort(failed_checks)
            println("- $check")
        end
    end
end
