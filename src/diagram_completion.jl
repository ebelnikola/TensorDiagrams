#################################################
# TOC
# topic 1:           DIAGRAM COMPLETION
#################################################

"""
    fill_labels(diag::TensorDiagram)

Fill in all question mark labels ("?") in a tensor diagram with all valid label combinations.

This function generates all possible complete diagrams by replacing "?" labels with valid
labels from each node's `allowed_labels`, while respecting each node's `forbidden_label_combinations`.

Arguments:
- diag: A TensorDiagram with potentially incomplete labels (containing "?")

Returns:
- Vector{TensorDiagram}: Array of complete diagrams with all valid label assignments

"""
function fill_labels(diag::TensorDiagram)
    # 1. Check if there are any question marks to be filled
    result = TensorDiagram[diag]

    pattern_labels = Set(vcat(diag.contraction_pattern...))
    boundary_labels = Set(vcat(values(diag.boundary_legs)...))
    all_labels = union(pattern_labels, boundary_labels)
    has_question_marks = any(lbl == "?" for lbl in values(diag.labels)) || keys(diag.labels) != all_labels

    if !has_question_marks
        return result
    end

    # 2. Iteratively fill in the question marks for each node
    for (i, node) in enumerate(diag.nodes)
        # for each node we fill the new_result array with diagrams whose indices surrounding the node are filled

        new_result = TensorDiagram[]
        for diag_loc in result

            # 2.1. Identify indices with "?" labels for the current node and extract allowed labels

            pattern_to_process = diag_loc.contraction_pattern[i]
            labels_of_pattern = [get(diag_loc.labels, ind, "?") for ind in pattern_to_process]
            mask = labels_of_pattern .== "?"
            allowed_labels = node.allowed_labels[mask]

            # 2.2. Run over all allowed labels, 
            #      check if the resulting pattern is consistent with "allowed_label_combinations" and "forbidden_label_combinations",
            #      and push the diagram into new_result if it is.  
            for sub_pattern in Iterators.product(allowed_labels...)
                labels_of_pattern[mask] .= sub_pattern # In case when there are no question marks, sub_pattern is empty, but mask is all false too, 
                # so this loop works appropriately even in this case. 

                pattern_is_allowed = all(labels_of_pattern[i] in lbls for (i, lbls) in enumerate(node.allowed_labels)) && # each label must be allowed on its own 
                                     (isempty(node.allowed_label_combinations) || # either we don't have a list of allowed combinations at all
                                      (labels_of_pattern in node.allowed_label_combinations)) # or the pattern belongs to this list 

                pattern_is_forbidden = labels_of_pattern in node.forbidden_label_combinations

                if !pattern_is_allowed || pattern_is_forbidden
                    continue
                end

                new_diag = copy(diag_loc)
                for j in eachindex(pattern_to_process)
                    new_diag.labels[pattern_to_process[j]] = labels_of_pattern[j]
                end
                push!(new_result, new_diag)
            end
        end
        result = new_result
    end
    return result
end;


