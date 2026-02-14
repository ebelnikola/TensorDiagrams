# TODO: Migrate to Dictionaries.jl to handle dictionaries
##################################################################
# TOC
# topic -1:          COLOR PALETTE
# topic 0:           CONFIG LOADER
# topic 1:           TensorNode AND TensorDiagram DEFINITIONS AND PLOTTING
#    - topic 1.1:    DEFINITIONS
#    - topic 1.2:    EQUALITY AND STANDARDIZATION
#    - topic 1.3:    PLOTTING AND PRINTING
#    - topic 1.4:    CONSISTENCY CHECKS
# topic 2:           TensorNode MANIPULATIONS
#    - topic 2.1:    COPYING
#    - topic 2.2:    REFLECTION
#    - topic 2.3:    LEGS FLIPPING
#    - topic 2.4:    ADJOINT OF A NODE
# topic 3:           TensorDiagram MANIPULATIONS
#    - topic 3.1:    COPYING
#    - topic 3.2:    DANGLING INDICES
#    - topic 3.3:    RELABELING
#    - topic 3.4:    REFLECTION
#    - topic 3.5:    COMPLETION
#    - topic 3.6:    ADJOINT OF A DIAGRAM
#    - topic 3.7:    GLUING
# topic 4:           TensorDiagram PROCESSING
#    - topic 4.1:    DIAGRAM TO CODE
#    - topic 4.2:    M AND TM GENERATION
# topic 5:           UTILS
# topic 6:           RANDOM GENERATORS
##################################################################
__precompile__(false) # TODO: This should disappear eventually when we 
#                     #       wrap everything into a module. 

##################################################################
# topic -1: COLOR PALETTE 
##################################################################

include("palette.jl")

##################################################################
# topic 0: CONFIG LOADER 
##################################################################
using YAML

"""
    load_config(filename::String)

Load configuration from a YAML file and set global constants.

# Arguments
- `filename::String`: Path to the YAML configuration file.

# Details
- Reads the YAML file.
- Converts "INFINITE_LABELS" and "FINITE_LABELS" entries to `Set{String}`.
- Iterates through the configuration dictionary and defines global constants for each key-value pair using `@eval`.

# Returns
- `config`: The loaded configuration dictionary.
"""
function load_config(filename::String)
    config = YAML.load_file(filename)
    config["INFINITE_LABELS"] = Set{String}(config["INFINITE_LABELS"])
    config["FINITE_LABELS"] = Set{String}(config["FINITE_LABELS"])

    for (var, val) in config
        @eval const $(Symbol(var)) = $val
    end
    return config
end

##################################################################
# topic 1: TensorNode AND TensorDiagram DEFINITIONS AND PLOTTING   
##################################################################

##################################################################
# topic 1.1: DEFINITIONS 
##################################################################


@kwdef mutable struct TensorNode
    # Essential parameters
    name::String
    legs::Dict{String,Array{Int,1}}
    allowed_labels::Array{Array{String,1},1}
    allowed_label_combinations::Array{Array{String,1},1} = Array{Array{String,1},1}[]
    forbidden_label_combinations::Array{Array{String,1},1} = Array{Array{String,1},1}[]
    hor_ref::Bool = false
    ver_ref::Bool = false
    # TensorKit related properties
    adjoint::Bool = false
    dual::Vector{Bool} = zeros(Bool, maximum(maximum.(values(legs))))
    flipped::Vector{Bool} = zeros(Bool, maximum(maximum.(values(legs))))
    out_legs_num::Int = maximum(maximum.(values(legs)))
    # Plotting parameters
    color::RGBA = YELLOW
    width::Union{Nothing,Float64} = nothing
    height::Union{Nothing,Float64} = nothing
end

@kwdef mutable struct TensorDiagram
    # Essential parameters
    nodes::Array{TensorNode,1}
    contraction_pattern::Array{Array{Int,1},1}
    boundary_legs_num::Dict{String,Int}
    boundary_legs::Dict{String,Array{Int,1}}
    boundary_legs_posidx::Dict{String,Array{Int,1}}
    labels::Dict{Int,String} = Dict{Int,String}()
    factor::String = "1.0"
    # Plotting parameters
    node_coordinates::Union{Nothing,Array{Union{Vector{Float64},Missing},1}} = nothing
    width::Float64 = 10.0
    height::Float64 = 10.0
end

##################################################################
# topic 1.2: EQUALITY AND STANDARDIZATION
##################################################################

import Base: ==
function ==(node1::TensorNode, node2::TensorNode)
    return node1.name == node2.name &&
           node1.legs == node2.legs &&
           node1.allowed_labels == node2.allowed_labels &&
           node1.allowed_label_combinations == node2.allowed_label_combinations &&
           node1.forbidden_label_combinations == node2.forbidden_label_combinations &&
           node1.hor_ref == node2.hor_ref &&
           node1.ver_ref == node2.ver_ref &&
           node1.adjoint == node2.adjoint &&
           node1.dual == node2.dual &&
           node1.flipped == node2.flipped &&
           node1.out_legs_num == node2.out_legs_num
end

include("diagram_standardizer.jl")

"""
    ==(diag1::TensorDiagram, diag2::TensorDiagram)

Check if two `TensorDiagram`s are equivalent.
Equivalence is determined by standardizing both diagrams using `standardize_diagram` and comparing their properties.
If standardization fails for either diagram (e.g. unreachable internal indices), this function returns `false`.

# Returns
- `true` if the diagrams are equivalent.
- `false` otherwise.
"""
function ==(diag1::TensorDiagram, diag2::TensorDiagram)
    std_diag1, success1 = standardize_diagram(diag1)
    std_diag2, success2 = standardize_diagram(diag2)


    if !success1 || !success2
        return false
    end

    # Check that sets of indices are the same
    indices1 = Set(vcat(std_diag1.contraction_pattern...))
    union!(indices1, vcat(values(std_diag1.boundary_legs)...))

    indices2 = Set(vcat(std_diag2.contraction_pattern...))
    union!(indices2, vcat(values(std_diag2.boundary_legs)...))

    if indices1 != indices2
        return false
    end

    # Populate labels with question marks for missing keys
    for idx in indices1
        if !haskey(std_diag1.labels, idx)
            std_diag1.labels[idx] = "?"
        end
        if !haskey(std_diag2.labels, idx)
            std_diag2.labels[idx] = "?"
        end
    end

    perm = sortperm(std_diag1.contraction_pattern)
    std_diag1.nodes = std_diag1.nodes[perm]
    std_diag1.contraction_pattern = std_diag1.contraction_pattern[perm]

    perm = sortperm(std_diag2.contraction_pattern)
    std_diag2.nodes = std_diag2.nodes[perm]
    std_diag2.contraction_pattern = std_diag2.contraction_pattern[perm]


    return std_diag1.nodes == std_diag2.nodes &&
           std_diag1.contraction_pattern == std_diag2.contraction_pattern &&
           std_diag1.boundary_legs_num == std_diag2.boundary_legs_num &&
           std_diag1.boundary_legs == std_diag2.boundary_legs &&
           std_diag1.boundary_legs_posidx == std_diag2.boundary_legs_posidx &&
           std_diag1.labels == std_diag2.labels &&
           std_diag1.factor == std_diag2.factor
end



##################################################################
# topic 1.3: PLOTTING AND PRINTING
##################################################################

include("diagram_plotting_utils.jl")
include("diagram_printing_utils.jl")



##################################################################
# topic 1.4: CONSISTENCY CHECKS
##################################################################
include("diagram_consistency.jl")


#########################################################
# topic 2: TensorNode MANIPULATIONS
#########################################################

#########################################################
# topic 2.1: COPYING
#########################################################


import Base: copy
function copy(node::TensorNode)
    return TensorNode(
        name=node.name,
        color=node.color,
        legs=deepcopy(node.legs), # Nested structure (dict of arrays) requires deepcopy
        allowed_labels=deepcopy(node.allowed_labels), # Nested structure (array of arrays) requires deepcopy
        forbidden_label_combinations=deepcopy(node.forbidden_label_combinations), # Nested structure (array of arrays) requires deepcopy
        allowed_label_combinations=deepcopy(node.allowed_label_combinations), # Nested structure (array of arrays) requires deepcopy
        hor_ref=node.hor_ref,
        ver_ref=node.ver_ref,
        adjoint=node.adjoint,
        dual=copy(node.dual),
        flipped=copy(node.flipped),
        out_legs_num=node.out_legs_num,
        width=node.width,
        height=node.height
    )
end

#########################################################
# topic 2.2: REFLECTION
#########################################################
# TODO: rewrite reflections using copy instead of explicitly defining all fields
"""
    hor_reflection(node::TensorNode)

Create a horizontally reflected version of the TensorNode.
- Swaps left and right legs
- Reverses order of top and bottom legs
- Sets hor_ref=true

# Arguments:
- node: The TensorNode to reflect

# Returns:
- A new TensorNode with horizontal reflection applied
"""
function hor_reflection(node::TensorNode)
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

    return TensorNode(
        name=node.name,
        color=node.color,
        legs=new_legs,
        allowed_labels=node.allowed_labels,
        forbidden_label_combinations=node.forbidden_label_combinations,
        allowed_label_combinations=node.allowed_label_combinations,
        hor_ref=!node.hor_ref,
        ver_ref=node.ver_ref,
        adjoint=node.adjoint,
        dual=copy(node.dual),
        flipped=copy(node.flipped),
        out_legs_num=node.out_legs_num,
        width=node.width,
        height=node.height
    )
end
"""
    ver_reflection(node::TensorNode)

Create a vertically reflected version of the TensorNode.
- Swaps top and bottom legs
- Reverses order of left and right legs
- Sets ver_ref=true

# Arguments:
- node: The TensorNode to reflect

# Returns:
- A new TensorNode with vertical reflection applied
"""
function ver_reflection(node::TensorNode)
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

    return TensorNode(
        name=node.name,
        color=node.color,
        legs=new_legs,
        allowed_labels=node.allowed_labels,
        forbidden_label_combinations=node.forbidden_label_combinations,
        allowed_label_combinations=node.allowed_label_combinations,
        hor_ref=node.hor_ref,
        ver_ref=!node.ver_ref,
        adjoint=node.adjoint,
        dual=copy(node.dual),
        flipped=copy(node.flipped),
        out_legs_num=node.out_legs_num,
        width=node.width,
        height=node.height
    )
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

#########################################################
# topic 2.3: LEGS FLIPPING
#########################################################

import TensorKit: flip
function flip!(node::TensorNode, idx::Int)
    node.flipped[idx] = !node.flipped[idx]
    return node
end
function flip(node::TensorNode, idx::Int)
    new_node = copy(node)
    flip!(new_node, idx)
    return new_node
end
function flip!(node::TensorNode, indicies::Vector{Int})
    node.flipped[indicies] .= .!node.flipped[indicies]
    return node
end
function flip(node::TensorNode, indicies::Vector{Int})
    new_node = copy(node)
    flip!(new_node, indicies)
    return new_node
end

#########################################################
# topic 2.4: ADJOINT OF A NODE
#########################################################


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
# topic 3: TensorDiagram MANIPULATIONS
#########################################################

#########################################################
# topic 3.1: COPYING
#########################################################

function copy(diag::TensorDiagram)
    return TensorDiagram(
        nodes=deepcopy.(diag.nodes), # Nested structure (Array of TensorNodes), requires deepcopy
        contraction_pattern=deepcopy(diag.contraction_pattern), # Nested structure (Array of Arrays) requires deepcopy
        boundary_legs_num=copy(diag.boundary_legs_num), # Not nested (Dict of Ints), copy is enough
        boundary_legs=deepcopy(diag.boundary_legs), # Nested structure (Dict of Arrays), requires deepcopy
        boundary_legs_posidx=deepcopy(diag.boundary_legs_posidx), # Nested structure (Dict of Arrays), requires deepcopy
        labels=copy(diag.labels),# Not nested (Dict of Strings), copy is enough  
        factor=diag.factor,
        node_coordinates=isnothing(diag.node_coordinates) ? nothing : deepcopy(diag.node_coordinates), # Nested structure (Array of Arrays), requires deepcopy
        width=diag.width,
        height=diag.height
    )
end

#########################################################
# topic 3.2: DANGLING INDICES
#########################################################


"""
    drop_dangling_indices!(diag::TensorDiagram, side::String)

Removes boundary legs (indices) from a specific side of a diagram that are not involved in the contraction pattern.
Modifies `diag` in place.

# Arguments
- `diag`: The diagram structure to modify.
- `side`: The side to process (e.g., "left", "right", "top", "bottom").

# Returns
- `diag`: The modified diagram.
"""
function drop_dangling_indices!(diag, side)
    all_indices_in_pattern = Set(vcat(diag.contraction_pattern...))
    if haskey(diag.boundary_legs, side)
        mask = [x ∈ all_indices_in_pattern for x in diag.boundary_legs[side]]
        to_remove_from_labels = diag.boundary_legs[side][.!mask]
        for key in to_remove_from_labels
            delete!(diag.labels, key)
        end
        diag.boundary_legs[side] = diag.boundary_legs[side][mask]
        diag.boundary_legs_posidx[side] = diag.boundary_legs_posidx[side][mask]
    end

    return diag
end
"""
    drop_dangling_indices(diag_init, side::String)

Removes boundary legs (indices) from a specific side of a diagram that are not involved in the contraction pattern. This function creates a copy of the input diagram and modifies the copy.


# Arguments
- `diag_init`: The initial diagram structure.

# Returns
- `diag`: A new diagram with dangling indices removed.
"""
function drop_dangling_indices(diag_init, side)
    diag = copy(diag_init)
    drop_dangling_indices!(diag, side)
    return diag
end
"""
    drop_dangling_indices!(diag::TensorDiagram)

Removes boundary legs (indices) from all sides of the diagram that are not involved in the contraction pattern.
Modifies `diag` in place.

# Arguments
- `diag`: The diagram structure to modify.

# Returns
- `diag`: The modified diagram.
"""
function drop_dangling_indices!(diag)
    for side in keys(diag.boundary_legs)
        drop_dangling_indices!(diag, side)
    end

    return diag
end
"""
    drop_dangling_indices(diag_init)

Removes boundary legs (indices) from a diagram that are not involved in the contraction pattern.
This function creates a copy of the input diagram and modifies the copy.

# Arguments
- `diag_init`: The initial diagram structure.

# Returns
- `diag`: A new diagram with dangling indices removed.
"""
function drop_dangling_indices(diag_init)
    diag = copy(diag_init)
    drop_dangling_indices!(diag)
    return diag
end

#########################################################
# topic 3.3: RELABELING
#########################################################


"""
    relabel!(diag::TensorDiagram, index_map::Dict{Int,Int})

Relabels the indices of the diagram in place using the provided dictionary.
Updates `contraction_pattern`, `boundary_legs`, and `labels`.
Indices not present in `index_map` are left unchanged.

# Arguments
- `diag`: The diagram to relabel.
- `index_map`: A dictionary mapping old indices to new indices.

# Returns
- `diag`: The modified diagram.
"""
function relabel!(diag::TensorDiagram, index_map::Dict{Int,Int})
    if length(unique(values(index_map))) != length(keys(index_map))
        throw(ArgumentError("The relabeling dictionary contains non-unique target indices."))
    end
    diag.contraction_pattern = [[get(index_map, idx, idx) for idx in pattern] for pattern in diag.contraction_pattern]
    diag.boundary_legs = Dict(side => [get(index_map, idx, idx) for idx in legs] for (side, legs) in diag.boundary_legs)
    diag.labels = Dict(get(index_map, idx, idx) => lbl for (idx, lbl) in diag.labels)
    return diag
end
"""
    relabel(diag::TensorDiagram, index_map::Dict{Int,Int})

Creates a copy of the diagram with indices relabeled using the provided dictionary.

# Arguments
- `diag`: The diagram to relabel.
- `index_map`: A dictionary mapping old indices to new indices.

# Returns
- A new `TensorDiagram` with relabeled indices.
"""
function relabel(diag::TensorDiagram, index_map::Dict{Int,Int})
    return relabel!(copy(diag), index_map)
end

#########################################################
# topic 3.4: REFLECTION
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
    rev_sides_num = dir == "vertical" ? diag.boundary_legs_num["horizontal"] : diag.boundary_legs_num["vertical"]


    diag_r.boundary_legs[flipping_sides[1]] = get(diag.boundary_legs, flipping_sides[2], Int[])
    diag_r.boundary_legs[flipping_sides[2]] = get(diag.boundary_legs, flipping_sides[1], Int[])
    diag_r.boundary_legs_posidx[flipping_sides[1]] = get(diag.boundary_legs_posidx, flipping_sides[2], Int[])
    diag_r.boundary_legs_posidx[flipping_sides[2]] = get(diag.boundary_legs_posidx, flipping_sides[1], Int[])



    diag_r.boundary_legs_posidx[reversing_sides[1]] = (rev_sides_num + 1) .- get(diag.boundary_legs_posidx, reversing_sides[1], Int[])
    diag_r.boundary_legs_posidx[reversing_sides[2]] = (rev_sides_num + 1) .- get(diag.boundary_legs_posidx, reversing_sides[2], Int[])

    if !isnothing(diag.node_coordinates)
        diag_r.node_coordinates .= (diag.node_coordinates) .|> x -> ismissing(x) ? missing : dir == "vertical" ? [x[1], -x[2]] : [-x[1], x[2]]
    end

    return diag_r
end

#########################################################
# topic 3.5: COMPLETION
#########################################################


include("diagram_completion.jl")

#########################################################
# topic 3.6: ADJOINT OF A DIAGRAM
#########################################################

# we get an adoint of diagram by taking adjoints of all 
# nodes and accounting for the leg permutations that 
# occur when taking adjoint of individual nodes.

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



#########################################################
# topic 3.7: GLUING
#########################################################


include("diagram_gluing.jl")


#########################################################
# topic 4: TensorDiagram PROCESSING
#########################################################
import TensorOperations
using LinearAlgebra: I

#########################################################
# topic 4.1: DIAGRAM TO CODE
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


#########################################################
# topic 4.2: M AND TM GENERATION
#########################################################

include("space_id.jl")

#########################################################
# topic 5: UTILS
#########################################################

include("utils.jl")


#########################################################
# topic 6: RANDOM GENERATORS
#########################################################

include("random_generators.jl")

