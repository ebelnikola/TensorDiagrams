#################################################
# TOC
# topic 0:           EXPORTS
# topic 1:           GENERAL UTILITIES
# topic 2:           CODE PROCESSING UTILITIES
#################################################

#################################################
# topic 0: EXPORTS
#################################################

export walk_through_dict_tree, generate_nfi!, extract_variables, remove_last_lines

#################################################
# topic 1: GENERAL UTILITIES
#################################################
"""
    walk_through_dict_tree(dict_tree, keys)

Traverse a nested dictionary structure following a sequence of keys.

# Arguments
- `dict_tree`: The nested dictionary to traverse
- `keys`: An iterable of keys to follow

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
- `free_index_array`: A 1-element array containing the current free index value. This array is modified in-place.

# Returns:
- `Int`: The current free index value before decrementing.
"""
function generate_nfi!(free_index_array; dir="negative")
    c = dir == "negative" ? -1 : 1
    nfi = free_index_array[1]
    free_index_array[1] = free_index_array[1] + c
    return nfi
end


################################################
# topic 2: CODE PROCESSING UTILITIES
################################################

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

