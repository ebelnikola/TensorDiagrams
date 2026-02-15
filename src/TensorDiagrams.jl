##################################################################
# TOC
# topic -1:          COLOR PALETTE
# topic 0:           CONFIG LOADER
# topic 1:           BASICS
#    - topic 1.1:    DEFINITIONS
#    - topic 1.2:    COPYING
#    - topic 1.3:    EQUALITY AND STANDARDIZATION
#    - topic 1.4:    CONSISTENCY CHECKS
#    - topic 1.5:    PLOTTING AND PRINTING
# topic 2:           GEOMETRICAL TRANSFORMATIONS
#    - topic 2.1:    REFLECTION
#    - topic 2.2:    ROTATION
# topic 3:           TENSOR OPERATIONS
#    - topic 3.1:    ADJOINT AND LEGS FLIPPING
# topic 4:           DIAGRAM COMPOSITION AND MODIFICATION
#    - topic 4.1:    GLUING
#    - topic 4.2:    RELABELING AND DANGLING INDICES
#    - topic 4.3:    COMPLETION
# topic 5:           PROCESSING
#    - topic 5.1:    DIAGRAM TO CODE
#    - topic 5.2:    M AND TM GENERATION
# topic 6:           UTILS
# topic 7:           RANDOM GENERATORS
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
# topic 1: BASICS
##################################################################

##################################################################
# topic 1.1: DEFINITIONS 
##################################################################
"""
    TensorNode

A structure representing a single tensor in a tensor network diagram.

# Fields

### Essential Parameters
- `name::String`: The name of the tensor.
- `legs::Dict{String,Array{Int,1}}`: A dictionary mapping sides to leg indices.
    - Allowed keys: "left", "right", "top", "bottom".
    - Values: Arrays of integers numbering the legs of the tensor, i.e. leg ordinals. 

### Diagram Completion Parameters
- `allowed_labels::Array{Array{String,1},1}`: Constraints on the vector space labels for each leg.
    - An array of arrays, where the i-th element corresponds to the i-th leg.
    - Each inner array contains the valid space labels (strings) that the leg can take.
    - Empty by default. 

- `allowed_label_combinations::Array{Array{String,1},1}`: A list of specific combinations of labels that are valid for the entire tensor.
    - If non-empty, the tensor's legs can only have labels that form one of these combinations.
    - Empty by default.
- `forbidden_label_combinations::Array{Array{String,1},1}`: A list of specific combinations of labels that are forbidden.
    - The tensor cannot have labels that form any of these combinations.
    - Empty by default.

### TensorKit Related Properties
- `adjoint::Bool`: Whether this node is the adjoint (conjugate transpose) of some original tensor. Default is `false`.
- `dual::Vector{Bool}`: A boolean vector indicating if each leg is in the dual space.
    - Length must match the total number of legs.
- `flipped::Vector{Bool}`: A boolean vector indicating if the arrow direction of each leg was flipped using `flip` from TensorKit.
    - Length must match the total number of legs.
- `out_legs_num::Int`: The number of outgoing legs.
    - Default is the total number of legs, i.e, all legs are outgoing.

### Plotting Parameters
- `hor_ref::Bool`: Indicates if the node has been horizontally reflected. Default `false`.
- `ver_ref::Bool`: Indicates if the node has been vertically reflected. Default `false`.
- `rotation::Int`: Rotation of the node in degrees (multiple of 90). Default 0.
- `color::RGBA`: The fill color of the node in plots. Default `YELLOW`, see "palette.jl".
- `width::Union{Nothing,Float64}`: Explicit width of the node for plotting. If `nothing`, width is calculated automatically.
- `height::Union{Nothing,Float64}`: Explicit height of the node for plotting. If `nothing`, height is calculated automatically.
"""
@kwdef mutable struct TensorNode
    # Required parameters
    name::String
    legs::Dict{String,Array{Int,1}}
    # Paramters for diagram completions
    allowed_labels::Array{Array{String,1},1} = map(_ -> String[], 1:maximum(maximum.(values(legs))))
    allowed_label_combinations::Array{Array{String,1},1} = Array{Array{String,1},1}[]
    forbidden_label_combinations::Array{Array{String,1},1} = Array{Array{String,1},1}[]
    # TensorKit related properties
    adjoint::Bool = false
    dual::Vector{Bool} = zeros(Bool, maximum(maximum.(values(legs))))
    flipped::Vector{Bool} = zeros(Bool, maximum(maximum.(values(legs))))
    out_legs_num::Int = maximum(maximum.(values(legs)))
    # Plotting parameters
    hor_ref::Bool = false
    ver_ref::Bool = false
    rotation::Int = 0
    color::RGBA = YELLOW
    width::Union{Nothing,Float64} = nothing
    height::Union{Nothing,Float64} = nothing
end

"""
    TensorDiagram

A structure representing a tensor network diagram, consisting of connected tensor nodes with legs decorated by labels indicating vector subspaces. Also contains a factor string, which should be a valid Julia expression that can be evaluated to a number. 

# Fields

### Required Parameters
- `nodes::Array{TensorNode,1}`: The list of `TensorNode` objects comprising the diagram.
- `contraction_pattern::Array{Array{Int,1},1}`: Specifies how the nodes are connected in the ncon fashion.
    - Each element corresponds to a node (in the same order as `nodes`).
    - The i-th element labels the legs of i-th node by pattern indices: positive indices for legs that are contracted (internal edges) and negative indices for legs that connect to the boundary.
- `boundary_slots_num::Dict{String,Int}`: Specifies the number of slots on each boundary side. A slot is a position where a leg can be attached.
    - Allowed keys: "left", "right", "top", "bottom".
- `boundary_legs::Dict{String,Array{Int,1}}`: Maps boundary sides to the pattern indices of the legs ending there.
    - Allowed keys: "left", "right", "top", "bottom".
    - Values: Arrays of pattern indices present on that boundary.
- `boundary_legs_posidx::Dict{String,Array{Int,1}}`: Specifies the geometric ordering of legs on the boundary.
    - Allowed keys: "left", "right", "top", "bottom".
    - Values: Arrays of integers from 1 to the number of slots on that side, indicating the position index, i.e. the slot ordinal, of each leg on its respective side. The slots are always ordered from bottom to top and from left to right. 

### Decorations
- `labels::Dict{Int,String}`: A dictionary mapping pattern indices to string labels.
    - These labels are meant to indicate particular subspaces of the total leg space. 

### Factor
- `factor::String`: A symbolic scalar factor multiplying the diagram. Expressed as a string that should be a valid Julia expression. Default "1.0".

### Plotting Parameters
- `node_coordinates::Union{Nothing,Array{Union{Vector{Float64},Missing},1}}`: Explicit 2D coordinates for the centers of the nodes.
    - If `nothing`, layout is computed automatically.
    - Supports `missing` for individual nodes to let the layout engine place them.
    - If not `nothing`, should contain cooridnates for each tensor. Otherwise the diagram will be rendered incorrectly. 
- `width::Float64`: The total width of the canvas/box representing the diagram boundaries. Default 10.0.
- `height::Float64`: The total height of the canvas/box representing the diagram boundaries. Default 10.0.
"""
@kwdef mutable struct TensorDiagram
    # Required parameters
    nodes::Array{TensorNode,1}
    contraction_pattern::Array{Array{Int,1},1}
    boundary_slots_num::Dict{String,Int}
    boundary_legs::Dict{String,Array{Int,1}}
    boundary_legs_posidx::Dict{String,Array{Int,1}}
    # Decorations
    labels::Dict{Int,String} = Dict{Int,String}()
    # Factor
    factor::String = "1.0"
    # Plotting parameters
    node_coordinates::Union{Nothing,Array{Union{Vector{Float64},Missing},1}} = nothing
    width::Float64 = 10.0
    height::Float64 = 10.0
end

##################################################################
# topic 1.2: COPYING
##################################################################
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
        rotation=node.rotation,
        adjoint=node.adjoint,
        dual=copy(node.dual),
        flipped=copy(node.flipped),
        out_legs_num=node.out_legs_num,
        width=node.width,
        height=node.height
    )
end

function copy(diag::TensorDiagram)
    return TensorDiagram(
        nodes=deepcopy.(diag.nodes), # Nested structure (Array of TensorNodes), requires deepcopy
        contraction_pattern=deepcopy(diag.contraction_pattern), # Nested structure (Array of Arrays) requires deepcopy
        boundary_slots_num=copy(diag.boundary_slots_num), # Not nested (Dict of Ints), copy is enough
        boundary_legs=deepcopy(diag.boundary_legs), # Nested structure (Dict of Arrays), requires deepcopy
        boundary_legs_posidx=deepcopy(diag.boundary_legs_posidx), # Nested structure (Dict of Arrays), requires deepcopy
        labels=copy(diag.labels),# Not nested (Dict of Strings), copy is enough  
        factor=diag.factor,
        node_coordinates=isnothing(diag.node_coordinates) ? nothing : deepcopy(diag.node_coordinates), # Nested structure (Array of Arrays), requires deepcopy
        width=diag.width,
        height=diag.height
    )
end

##################################################################
# topic 1.3: EQUALITY AND STANDARDIZATION
##################################################################

import Base: ==
function ==(node1::TensorNode, node2::TensorNode)
    node1 = unreflect(unrotate(node1))
    node2 = unreflect(unrotate(node2))
    return node1.name == node2.name &&
           node1.legs == node2.legs &&
           node1.allowed_labels == node2.allowed_labels &&
           node1.allowed_label_combinations == node2.allowed_label_combinations &&
           node1.forbidden_label_combinations == node2.forbidden_label_combinations &&
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
- `true` if the diagrams are equivalent and all internal indices are reachable from the boundary.
- `false` otherwise. 
"""
function ==(diag1::TensorDiagram, diag2::TensorDiagram)
    std_diag1, success1 = standardize_diagram(diag1)
    std_diag2, success2 = standardize_diagram(diag2)


    if !success1 || !success2
        return false
    end

    return std_diag1.nodes == std_diag2.nodes &&
           std_diag1.contraction_pattern == std_diag2.contraction_pattern &&
           std_diag1.boundary_slots_num == std_diag2.boundary_slots_num &&
           std_diag1.boundary_legs == std_diag2.boundary_legs &&
           std_diag1.boundary_legs_posidx == std_diag2.boundary_legs_posidx &&
           std_diag1.labels == std_diag2.labels &&
           std_diag1.factor == std_diag2.factor
end

##################################################################
# topic 1.4: CONSISTENCY CHECKS
##################################################################
include("diagram_consistency.jl")

##################################################################
# topic 1.5: PLOTTING AND PRINTING
##################################################################

include("diagram_plotting_utils.jl")
include("diagram_printing_utils.jl")


#########################################################
# topic 2: GEOMETRICAL TRANSFORMATIONS
#########################################################

#########################################################
# topic 2.1: REFLECTION
#########################################################
include("diagram_reflection.jl")

#########################################################
# topic 2.2: ROTATION
#########################################################
include("diagram_rotation.jl")

#########################################################
# topic 3: TENSOR OPERATIONS
#########################################################

#########################################################
# topic 3.1: ADJOINT AND LEGS FLIPPING
#########################################################
include("diagram_adjoints_and_flips.jl")

#########################################################
# topic 4: DIAGRAM COMPOSITION AND MODIFICATION
#########################################################

#########################################################
# topic 4.1: GLUING
#########################################################
include("diagram_gluing.jl")

#########################################################
# topic 4.2: RELABELING AND DANGLING INDICES
#########################################################
# Moved to diagram_standardizer.jl 

#########################################################
# topic 4.3: COMPLETION
#########################################################
include("diagram_completion.jl")

#########################################################
# topic 5: PROCESSING
#########################################################

#########################################################
# topic 5.1: DIAGRAM TO CODE
#########################################################
include("diagram_to_code.jl")

#########################################################
# topic 5.2: M AND TM GENERATION
#########################################################
include("space_id.jl")

#########################################################
# topic 6: UTILS
#########################################################
include("utils.jl")

#########################################################
# topic 7: RANDOM GENERATORS
#########################################################
include("random_generators.jl")
