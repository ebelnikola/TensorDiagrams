#################################################
# TOC
# topic 0:         EXPORTS
# topic 1:         SPACE ID DEFINITION AND UTILS
# topic 2:         SPACE ID TO LAYOUT AND SPACE
#################################################

##################################################
# topic 0: EXPORTS
##################################################

export SpaceID, string_to_space_id, is_valid_space_id, get_space_id, @sid_str
export space_ids_to_layout, space_id_to_space, space_id_to_dim

##################################################
# topic 1: SPACE ID DEFINITION AND UTILS
##################################################

"""
    SpaceID(labels, posidx, slots_number)

Identifies a vector space on a boundary of a tensor diagram.

# Fields
- `labels::Vector{String}`: The labels occupying the slots (e.g. `"z"`, `"q"`).
- `posidx::Vector{Int}`: Position indices of each label within the boundary slots.
- `slots_number::Int`: Total number of slots on the boundary (including empty ones).

Use [`string_to_space_id`](@ref) or the `sid"..."` string macro to construct from a string representation,
where labels occupy positions and `'-'` denotes an empty slot (e.g. `sid"z-q"`).
"""
struct SpaceID
    labels::Vector{String}
    posidx::Vector{Int}
    slots_number::Int
end

import Base: ==
function ==(sid1::SpaceID, sid2::SpaceID)
    return sid1.labels == sid2.labels && sid1.posidx == sid2.posidx && sid1.slots_number == sid2.slots_number
end

import Base: hash
function hash(sid::SpaceID, h::UInt)
    return hash(sid.slots_number, hash(sid.posidx, hash(sid.labels, h)))
end

function Base.string(sid::SpaceID)
    slots = fill("-", sid.slots_number)
    for (label, pos) in zip(sid.labels, sid.posidx)
        if 1 <= pos <= sid.slots_number
            slots[pos] = label
        else
            @warn "Label $label at position $pos is out of bounds for SpaceID with $(sid.slots_number) slots."
        end
    end
    return join(slots)
end

function Base.show(io::IO, sid::SpaceID)
    print(io, "sid:" * string(sid))
end

"""
    string_to_space_id(s::String) -> SpaceID

Parse a string representation into a `SpaceID`.
Labels are matched greedily (longest first) from the set of valid labels (`FINITE_LABELS ∪ INFINITE_LABELS`).
A `'-'` character denotes an empty slot.

# Examples
```julia
string_to_space_id("z-q")   # labels=["z","q"], posidx=[1,3], slots_number=3
string_to_space_id("zz")    # labels=["z","z"], posidx=[1,2], slots_number=2
```
"""
function string_to_space_id(s)
    labels = String[]
    posidx = Int[]

    valid_labels = union(FINITE_LABELS, INFINITE_LABELS)
    # Sort valid_labels by length descending to match longest first (e.g. z' before z)
    sorted_valid_labels = sort(collect(valid_labels), by=length, rev=true)

    current_pos = 1
    i = 1 # char index in string
    str_len = length(s)

    while i <= str_len
        if s[i] == '-'
            # Empty slot
            current_pos += 1
            i += 1
            continue
        end

        # Try to match a label
        matched = false
        for lbl in sorted_valid_labels
            len = length(lbl)
            if i + len - 1 <= str_len && s[i:i+len-1] == lbl
                push!(labels, lbl)
                push!(posidx, current_pos)
                current_pos += 1
                i += len
                matched = true
                break
            end
        end

        if !matched
            error("Could not parse SpaceID from string '$s': invalid character/label at index $i")
        end
    end

    return SpaceID(labels, posidx, current_pos - 1)
end

"""
    is_valid_space_id(s::String)

Check if a given string represents a valid SpaceID.
"""
function is_valid_space_id(s)
    try
        string_to_space_id(s)
        return true
    catch
        return false
    end
end

macro sid_str(str)
    return :(string_to_space_id($str))
end

"""
    get_space_id(diag; side="top")

Identify the vector space on a specific boundary of a diagram.
The identity is defined by the sorted sequence of labels and their position indices.

# Arguments
- `diag`: The diagram to inspect.
- `side`: The boundary side to check ("top", "left", "right", "bottom"). Default is "top".

# Returns
- A `SpaceID` object.
"""
function get_space_id(diag; side="top")
    legs_idx = diag.boundary_legs[side]
    legs_posidx = diag.boundary_legs_posidx[side]
    legs_labels = [split(get(diag.labels, idx, "?"), "←")[1] for idx in legs_idx]
    perm = sortperm(legs_posidx)
    legs_labels = legs_labels[perm]
    legs_posidx = legs_posidx[perm]

    slots_number = diag.boundary_slots_num[side]

    return SpaceID(legs_labels, legs_posidx, slots_number)
end


##################################################
# topic 2: SPACE ID TO LAYOUT AND SPACE
##################################################

"""
    space_ids_to_layout(space_ids; space_dims_dict::Dict{String,Int})

Generate a layout mapping for a set of space IDs. Each space ID is assigned a unique range of indices.

# Arguments
- `space_ids`: A collection of space identifiers (`SpaceID` objects).
- `space_dims_dict`: A dictionary mapping space labels to their dimensions.

# Returns
- A dictionary mapping each space ID to a `UnitRange{Int}` representing its location in the layout.
"""
function space_ids_to_layout(space_ids; space_dims_dict::Dict{String,Int})
    space_layout = Dict{SpaceID,UnitRange{Int}}()
    cnt = 1
    for space_id in space_ids
        local_finite_labels = Set(keys(space_dims_dict) ∩ space_id.labels)

        powers = Dict(space => count(==(space), space_id.labels) for space in local_finite_labels)

        dim = prod(space_dims_dict[space]^powers[space] for space in local_finite_labels; init=1)

        space_layout[space_id] = cnt:(cnt+dim-1)
        cnt += dim
    end
    space_layout
end

import TensorKit: sectors, dim, Sector, fuse
"""
    space_ids_to_layout(space_ids; spaces_dict::Dict{String,ElementarySpace})

Generate a layout mapping for a set of space IDs using TensorKit spaces. 
Each space ID is assigned a dictionary mapping each irrep (Sector) in the fused space to its range of indices.

# Arguments
- `space_ids`: A collection of space identifiers (`SpaceID` objects).
- `spaces_dict`: A dictionary mapping space labels to TensorKit elementary spaces.

# Returns
- A dictionary mapping each `SpaceID` to a `Dict{Sector, UnitRange{Int}}`.
"""
function space_ids_to_layout(space_ids; spaces_dict::Dict{String,ElementarySpace})

    S = sectortype(first(values(spaces_dict)))

    space_layout = Dict{SpaceID,Dict{S,UnitRange{Int}}}()
    cnt = Dict{S,Int}()



    for space_id in space_ids
        # 1. Construct the full product space for this SpaceID
        #    space_id_to_space already handles the logic of looking up labels in spaces_dict
        #    and tensor-producting them together.
        product_space = fuse(space_id_to_space(space_id; spaces_dict=spaces_dict))

        # 2. Iterate over sectors present in this product space
        #    We want to assign a range for each sector based on its dimension in the fused space.
        #    Note: 'sectors(product_space)' iterates over the coupled sectors.

        current_id_layout = Dict{S,UnitRange{Int}}()

        for s in sectors(product_space)
            d = dim(product_space, s)
            cnt_local = get(cnt, s, 1)
            current_id_layout[s] = cnt_local:(cnt_local+d-1)
            cnt[s] = cnt_local + d
        end

        space_layout[space_id] = current_id_layout
    end

    return space_layout
end

import TensorKit: ⊗, one

"""
    space_id_to_space(space_id::SpaceID; spaces_dict)

Construct the tensor product space for the labels in `space_id` that are present in `spaces_dict`.
The labels are ordered by their position index (`posidx`) and mapped to spaces using `spaces_dict`.
Labels not found in `spaces_dict` (e.g. infinite-dimensional labels) are skipped.

# Arguments
- `space_id::SpaceID`: The space ID to convert.
- `spaces_dict`: A dictionary mapping space labels to TensorKit space objects.

# Returns
- A TensorKit space object representing the product space.
"""
function space_id_to_space(space_id::SpaceID; spaces_dict)
    # Sort labels by position index
    perm = sortperm(space_id.posidx)
    sorted_labels = space_id.labels[perm]

    # Map to spaces using spaces_dict (only those present in the dict)
    spaces = [spaces_dict[lbl] for lbl in sorted_labels if haskey(spaces_dict, lbl)]

    if isempty(spaces)
        if isempty(spaces_dict)
            error("Cannot determine unit space: spaces_dict is empty.")
        end
        return one(first(values(spaces_dict)))
    end

    return reduce(⊗, spaces)
end

"""
    space_id_to_dim(space_id::SpaceID; dims_dict) -> Number

Compute the total dimension of the space identified by `space_id`.
Only labels present as keys in `dims_dict` are included; other labels
(e.g. infinite-dimensional ones) are skipped.
The dimension is the product of `dims_dict[lbl]` for each matching label,
accounting for multiplicity (repeated labels contribute multiple factors).

# Arguments
- `space_id::SpaceID`: The space identifier.
- `dims_dict`: A dictionary mapping label strings to their integer dimensions (e.g. `Dict("z" => 2, "q" => 8)`).

# Returns
- The product of dimensions of matching labels. Returns `1` if no labels match.

# Examples
```julia
sid = string_to_space_id("zqz")
space_id_to_dim(sid; dims_dict=Dict("z" => 2, "q" => 8))  # 2 * 8 * 2 = 32
```
"""
function space_id_to_dim(space_id::SpaceID; dims_dict)
    prod([dims_dict[lbl] for lbl in space_id.labels if haskey(dims_dict, lbl)])
end