#################################################
# TOC
# topic 0:           GENERAL UTILITIES
# topic 1:           TENSORS MANIPULATION UTILITIES
# topic 2:           DIAGRAM PROCESSING RELATED UTILITIES
# topic 3:           PRETTY VISUALS
#################################################

#################################################
# topic 0: GENERAL UTILITIES
#################################################
"""
    walk_through_dict_tree(dict_tree, keys)

Traverse a nested dictionary structure following a sequence of keys.

# Arguments
- dict_tree: The nested dictionary to traverse
- keys: An iterable of keys to follow

# Returns:
- The value found at the end of the path.
"""
function walk_through_dict_tree(dict_tree, keys)
    for key in keys
        dict_tree = dict_tree[key]
        if !(dict_tree isa Dict)
            return dict_tree
        end
    end
    return dict_tree
end

"""
    generate_nfi!(free_index_array)

Generate the next free index and update the counter.

# Arguments
- free_index_array: A 1-element array containing the current free index value. This array is modified in-place.

# Returns:
- Int: The current free index value before decrementing.
"""
function generate_nfi!(free_index_array; dir="negative")
    c = dir == "negative" ? -1 : 1
    nfi = free_index_array[1]
    free_index_array[1] = free_index_array[1] + c
    return nfi
end

"""
    extract_variables(expr_str::String)

Parse a Julia expression string and extract all variable names used in it.
Ignores function names in function calls (e.g., in "f(x)", extracts "x" but not "f").
Handles field access (e.g., in "obj.field", extracts "obj").

# Arguments
- `expr_str::String`: A string containing a valid Julia expression.

# Returns
- `Set{String}`: A set of unique variable names found in the expression.
"""
function extract_variables(expr_str::String)
    expr = Meta.parse(expr_str)
    vars = Set{String}()
    _extract_vars_from_expr!(vars, expr)
    return vars
end

function _extract_vars_from_expr!(vars, x)
    # Fallback for literals and other types (Int, Float, String, QuoteNode, nothing, etc.)
    return
end

function _extract_vars_from_expr!(vars, s::Symbol)
    push!(vars, string(s))
end

function _extract_vars_from_expr!(vars, ex::Expr)
    if ex.head == :call
        # Check if the function part is complex (e.g. d.f(o))
        # If it is a Symbol, we treat it as a function name and skip it.
        # If it is an Expr (like field access), we extract variables from it.
        if !(ex.args[1] isa Symbol)
            _extract_vars_from_expr!(vars, ex.args[1])
        end

        # Gather variables from remaining arguments
        for arg in ex.args[2:end]
            _extract_vars_from_expr!(vars, arg)
        end
    elseif ex.head == :.
        # Handle field access: obj.field -> extract from obj.
        # ex.args[1] is the object. ex.args[2] is the field (QuoteNode), which falls back to nothing.
        _extract_vars_from_expr!(vars, ex.args[1])
    elseif ex.head == :macrocall
        # Skip macro name (1st) and LineNumberNode (2nd)
        for arg in ex.args[3:end]
            _extract_vars_from_expr!(vars, arg)
        end
    else
        # Generic traversal for other expression types (including :ref, :=, etc.)
        for arg in ex.args
            _extract_vars_from_expr!(vars, arg)
        end
    end
end


"""
    remove_last_lines(str, n::Int)

Removes the last `n` lines from the input string `str`.

# Arguments
- `str`: The input string, potentially multi-line.
- `n::Int`: The number of lines to remove from the end.

# Returns
- A string with the last `n` lines removed. If `n` is greater than or equal to the total number of lines, returns an empty string.
"""
function remove_last_lines(str, n::Int)
    lines = split(str, '\n')
    if n >= length(lines)
        return ""
    else
        return join(lines[1:end-n], '\n')
    end
end


#################################################
# topic 1: TENSORS MANIPULATION UTILITIES
#################################################


using ChainRulesCore

"""
    combine_tensors_by_layout(tensors, layout)

Combines a list of tensors into a single larger tensor based on the specified layout.

# Arguments
- `tensors`: A collection of tensors (arrays) to be combined.
- `layout`: A collection of tuples of ranges, where each tuple specifies the location of the corresponding tensor in the result.
          Each tuple must have the same length as the number of dimensions of the resulting tensor.

# Returns
- `res`: The combined tensor with zeros in unoccupied regions.

# Safety
This function checks for overlapping regions in the layout unless `@inbounds` is used.
Overlapping regions will cause an error.
"""
function combine_tensors_by_layout(tensors, layout)
    @boundscheck begin
        for i in 1:length(layout)
            for j in i+1:length(layout)
                loc1 = layout[i]
                loc2 = layout[j]
                overlap = true
                for k in eachindex(loc1)
                    if isempty(intersect(loc1[k], loc2[k]))
                        overlap = false
                        break
                    end
                end
                if overlap
                    error("Layouts for tensor $i and tensor $j overlap.")
                end
            end
        end
    end

    max_dims = [maximum(last.(layout .|> x -> x[i])) for i in eachindex(layout[1])]
    res = zeros(max_dims...)
    for (tensor, loc) in zip(tensors, layout)
        res[loc...] .= tensor
    end
    return res
end

function ChainRulesCore.rrule(::typeof(combine_tensors_by_layout), tensors, layout)
    y = combine_tensors_by_layout(tensors, layout)
    function combine_tensors_by_layout_pullback(ȳ)
        y_bar = unthunk(ȳ)
        tensors_bar = [y_bar[loc...] for loc in layout]
        return NoTangent(), tensors_bar, NoTangent()
    end
    return y, combine_tensors_by_layout_pullback
end



#################################################
# topic 2: DIAGRAM PROCESSING RELATED UTILITIES
#################################################

"""
    space_ids_to_layout(space_ids; space_dims_dict::Dict{String,Int}=Dict("z" => 2))

Generate a layout mapping for a set of space IDs. Each space ID is assigned a unique range of indices.

# Arguments
- `space_ids`: A collection of space identifiers (`SpaceID` objects).
- `space_dims_dict`: A dictionary mapping space labels to their dimensions. Default is `Dict("z" => 2)`.

# Returns
- A dictionary mapping each space ID to a `UnitRange{Int}` representing its location in the layout.

#SR - DONE 10.2.2026 Can we correct this function such that it would not use the global array FINITE_LABELS ? 
let it just use dimensions from the provided dictionary and if some key is missing ignore this label, 
    i.e. assign dimension 1
"""
function space_ids_to_layout(space_ids; space_dims_dict::Dict{String,Int}=Dict("z" => 2))
    space_layout = Dict{SpaceID,UnitRange{Int}}()
    space_ids = space_ids
    cnt = 1
    for space_id in space_ids
        #local_finite_labels = Set(FINITE_LABELS ∩ space_id.labels)
        local_finite_labels = Set(keys(space_dims_dict) ∩ space_id.labels)

        powers = Dict(space => count(==(space), space_id.labels) for space in local_finite_labels)

        dim = prod(space_dims_dict[space]^powers[space] for space in local_finite_labels; init=1)

        space_layout[space_id] = cnt:(cnt+dim-1)
        cnt += dim
    end
    space_layout
end

"""
    space_id_to_space(space_id::SpaceID; spaces_dict)

Construct the product space corresponding to the finite-dimensional labels in `space_id`.
The labels are ordered by their position index (`posidx`) and mapped to spaces using `spaces_dict`.

# Arguments
- `space_id::SpaceID`: The space ID to convert.
- `spaces_dict`: A dictionary mapping space labels to TensorKit space objects (default: Dict("z" => Z2Space(1, 1))).

# Returns
- A TensorKit space object representing the product space.
"""
function space_id_to_space(space_id::SpaceID; spaces_dict=Dict("z" => Z2Space(1, 1)))
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

#################################################
# topic 3: PRETTY VISUALS
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

