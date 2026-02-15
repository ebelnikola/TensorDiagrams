#################################################
# TOC
# topic 1:           LEGS FLIPPING (NODE)
# topic 2:           ADJOINT (NODE)
# topic 3:           ADJOINT (DIAGRAM)
#################################################



#########################################################
# topic 1: LEGS FLIPPING (NODE)
#########################################################

import TensorKit: flip

"""
    flip!(node::TensorNode, idx::Int)

Flip the direction of a specific leg of a `TensorNode` in place.
Updates the `flipped` flag for the given index.

# Arguments
- `node::TensorNode`: The node to modify.
- `idx::Int`: The index of the leg to flip.

# Returns
- `node`: The modified node.
"""
function flip!(node::TensorNode, idx::Int)
    node.flipped[idx] = !node.flipped[idx]
    return node
end

"""
    flip(node::TensorNode, idx::Int)

Create a copy of a `TensorNode` with a specific leg flipped.

# Arguments
- `node::TensorNode`: The node to copy and modify.
- `idx::Int`: The index of the leg to flip.

# Returns
- A new `TensorNode` with the flipped leg.
"""
function flip(node::TensorNode, idx::Int)
    new_node = copy(node)
    flip!(new_node, idx)
    return new_node
end

"""
    flip!(node::TensorNode, indicies::Vector{Int})

Flip the direction of multiple legs of a `TensorNode` in place.

# Arguments
- `node::TensorNode`: The node to modify.
- `indicies::Vector{Int}`: A list of leg indices to flip.

# Returns
- `node`: The modified node.
"""
function flip!(node::TensorNode, indicies::Vector{Int})
    node.flipped[indicies] .= .!node.flipped[indicies]
    return node
end

"""
    flip(node::TensorNode, indicies::Vector{Int})

Create a copy of a `TensorNode` with multiple legs flipped.

# Arguments
- `node::TensorNode`: The node to copy and modify.
- `indicies::Vector{Int}`: A list of leg indices to flip.

# Returns
- A new `TensorNode` with the flipped legs.
"""
function flip(node::TensorNode, indicies::Vector{Int})
    new_node = copy(node)
    flip!(new_node, indicies)
    return new_node
end

#########################################################
# topic 2: ADJOINT OF A NODE
#########################################################


"""
    adjoint(node::TensorNode)

Compute the adjoint of a `TensorNode`.
This involves:
1. Adding or removing a prime ("'") from the name.
2. Toggling the `adjoint` flag.
3. Swapping input and output legs (reordering `legs`).
4. Reordering properties (`allowed_labels`, `dual`, `flipped`) to match the new leg order.
5. Updating `out_legs_num`.

# Arguments
- `node::TensorNode`: The node to take the adjoint of.

# Returns
- `new_node::TensorNode`: The adjoint node.
- `old_legs_to_new::Vector{Int}`: Validation permutation vector mapping old leg indices to new leg indices.
"""
function adjoint(node::TensorNode)
    new_node = copy(node)

    # 1. Add or remove a prime from the name.
    if endswith(new_node.name, "'")
        new_node.name = new_node.name[1:end-1]
    else
        new_node.name = new_node.name * "'"
    end

    # 2. Signal that this is an adjoint node
    new_node.adjoint = !node.adjoint


    # 3. Reorder legs: out legs become in legs and vice versa.
    # 3.1 Create mapping from old leg indices to new leg indices
    legs_num = maximum(maximum.(values(node.legs)))
    out_spaces = node.out_legs_num
    in_spaces = legs_num - out_spaces
    old_legs_to_new = vcat(collect(in_spaces+1:in_spaces+out_spaces), collect(1:in_spaces))


    # 3.2 Update legs by replacing old indices with new indices
    for side in keys(new_node.legs)
        new_node.legs[side] = old_legs_to_new[new_node.legs[side]]
    end

    # 3.3 For arrays whose indices are interpreted as leg indices, use the map from old to new 
    #     to reorder them
    new_legs_to_old = invperm(old_legs_to_new)
    new_node.allowed_labels = new_node.allowed_labels[new_legs_to_old]
    new_node.allowed_label_combinations = [combination[new_legs_to_old] for combination in new_node.allowed_label_combinations]
    new_node.forbidden_label_combinations = [combination[new_legs_to_old] for combination in new_node.forbidden_label_combinations]
    new_node.dual = new_node.dual[new_legs_to_old]
    new_node.flipped = new_node.flipped[new_legs_to_old]

    # 3.4 Update out_legs_num
    new_node.out_legs_num = in_spaces

    return new_node, old_legs_to_new
end




#########################################################
# topic 3: ADJOINT OF A DIAGRAM
#########################################################

"""
    adjoint(diag::TensorDiagram)

Compute the adjoint of a `TensorDiagram`.
This involves:
1. Taking the adjoint of all constituent nodes.
2. Updating the contraction pattern to account for the leg permutations incurred by node adjoints.

# Arguments
- `diag::TensorDiagram`: The diagram to take the adjoint of.

# Returns
- A new `TensorDiagram` representing the adjoint.
"""
function adjoint(diag::TensorDiagram)
    diag_a = copy(diag)

    nodes_and_perms = adjoint.(diag_a.nodes)

    new_nodes = nodes_and_perms .|> x -> x[1]
    old_legs_to_new_perms = nodes_and_perms .|> x -> x[2]
    new_legs_to_old_perms = invperm.(old_legs_to_new_perms)

    diag_a.nodes = new_nodes
    for i in eachindex(diag_a.contraction_pattern)
        diag_a.contraction_pattern[i] = diag_a.contraction_pattern[i][new_legs_to_old_perms[i]]
    end

    return diag_a
end
