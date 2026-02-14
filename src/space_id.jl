
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

    slots_key = (side == "left" || side == "right") ? "horizontal" : "vertical"
    slots_number = diag.boundary_legs_num[slots_key]

    return SpaceID(legs_labels, legs_posidx, slots_number)
end