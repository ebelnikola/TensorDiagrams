using Test
using Aqua

include("../src/TensorDiagrams.jl")
load_config((@__DIR__) * "/tests_config.yml")


###################################################################
# topic 1: LOADER AND RANDOM GENERATORS
###################################################################

@testset "YAML LOADER                     " begin
    # Tests the configuration loader to ensure  sets of labels
    # (finite vs infinite) are correctly parsed from the YAML file.
    @test INFINITE_LABELS == Set(["x", "u", "d", "r"])
    @test FINITE_LABELS == Set(["z", "q", "g", "z'", "t", "w", "c"])
end

@testset "GENERATION OF RANDOM NODES      " begin
    # Checks basic consistency of generated node.
    node_consistency = true
    # Checks if requested necessary_labels are indeed present in allowed_labels.
    necessary_labels_added = true
    # Checks consistency after adding forbidden combinations (should still be consistent if allowed ones are empty).
    node_consistency_after_adding_forbidden = true
    # Checks consistency after adding allowed combinations (should still be consistent if forbidden ones are empty).
    node_consistency_after_adding_allowed = true
    # Checks that having both allowed and forbidden combinations is flagged as inconsistent.
    node_inconsistency_after_adding_both = true

    N = 1000
    for _ = 1:N

        necessary_labels = rand(INFINITE_LABELS ∪ FINITE_LABELS, rand(0:2))

        T = generate_random_tensor_node(rand(1:20), "T"; max_num_of_allowed_rnd_labels=rand(4:1000), necessary_labels=necessary_labels)

        check_results = check_node_consistency(T)
        if !all(values(check_results))
            node_consistency = false
        end

        if !all([necessary_labels ⊆ leg_al for leg_al in T.allowed_labels])
            necessary_labels_added = false
        end

        Tp = copy(T)
        add_random_forbidden_label_combinations!(Tp; max_num_forbidden_combinations=rand(0:1000))

        check_results = check_node_consistency(Tp)
        if !all(values(check_results))
            node_consistency_after_adding_forbidden = false
        end

        Tp = copy(T)
        add_random_allowed_label_combinations!(Tp; max_num_allowed_combinations=rand(0:1000))

        check_results = check_node_consistency(Tp)
        if !all(values(check_results))
            node_consistency_after_adding_allowed = false
        end

        add_random_forbidden_label_combinations!(Tp; max_num_forbidden_combinations=rand(1:1000))
        add_random_allowed_label_combinations!(Tp; max_num_allowed_combinations=rand(1:1000))

        check_results = check_node_consistency(Tp)
        if all(values(check_results))
            node_inconsistency_after_adding_both = false
        end
    end
    @test node_consistency
    @test necessary_labels_added
    @test node_consistency_after_adding_forbidden
    @test node_consistency_after_adding_allowed
    @test node_inconsistency_after_adding_both
end

@testset "GENERATION OF RANDOM DIAGRAMS   " begin
    # Checks that generated diagrams pass all consistency checks.
    diagram_consistency = true
    # Checks that the generator throws an error when constraints are impossible.
    error_handling = true

    N = 1000
    for _ = 1:N
        boundary_slots_num = rand(0:20)
        max_nodes_num = rand(1:10)
        max_leg_num = rand(1:8)

        # Expected error condition: not enough slots to fit boundary legs
        # Same logic as in generate_random_tensor_diagram
        min_nodes_required = max(ceil(Int, boundary_slots_num / max_leg_num), 1)
        expect_error = min_nodes_required > max_nodes_num

        if expect_error
            try
                generate_random_tensor_diagram(boundary_slots_num; max_nodes_num=max_nodes_num, max_leg_num=max_leg_num)
                # If we are here, no error was thrown.
                error_handling = false
            catch
                # Expected error caught
            end
        else
            d = generate_random_tensor_diagram(boundary_slots_num; max_nodes_num=max_nodes_num, max_leg_num=max_leg_num)
            check_results = check_diagram_consistency(d)
            if !all(values(check_results))
                diagram_consistency = false
            end
        end
    end

    @test diagram_consistency
    @test error_handling
end

###################################################################
# topic 2: EQUALITIES AND COPYING
###################################################################

@testset "TENSOR NODE EQUALITY AND COPYING" begin
    # Tests reflexivity: a node should be equal to itself.
    equality_reflexivity = true
    # Tests inequality: different nodes should not be equal.
    equality_inequality = true
    # Tests that a copy is equal to the original.
    copy_equality = true
    # Tests that metadata (like 'width') is ignored in equality checks.
    metadata_ignorance = true


    N = 1000
    for _ = 1:N
        necessary_labels = rand(INFINITE_LABELS ∪ FINITE_LABELS, rand(0:2))

        T1 = generate_random_tensor_node(rand(1:20), "T"; max_num_of_allowed_rnd_labels=rand(4:1000), necessary_labels=necessary_labels)


        T2 = copy(T1)
        add_random_allowed_label_combinations!(T2; max_num_allowed_combinations=3)

        if T1 != T1
            equality_reflexivity = false
        end
        if T1 == T2
            equality_inequality = false
        end

        T1p = copy(T1)
        if T1 != T1p
            copy_equality = false
        end

        T1p.width = 2.0
        T1p.height = 3.0
        T1p.color = BLUE
        if T1 != T1p
            metadata_ignorance = false
        end
    end

    @test equality_reflexivity
    @test equality_inequality
    @test copy_equality
    @test metadata_ignorance
end


@testset "DIAGRAM EQUALITY AND COPYING    " begin
    # Tests reflexivity: a diagram should be equal to itself.
    reflexivity = true
    # Tests that a copy is equal to the original.
    copy_equality = true
    # Tests that setting a label to "?" (unknown/default) does not break equality.
    unknown_label_equality = true
    # Tests that modifying the structure (reflecting a node) breaks equality.
    structure_inequality = true
    # Test that setting a specific label (different from "?") breaks equality.
    label_inequality = true
    # Test that relabeling indices of the pattern preserves equality.
    relabeling_equality = true
    # Test that permuting nodes preserves equality.
    permutation_equality = true

    # Test that adding dangling indices breaks equality.
    dangling_indices_break_equality = true


    bubble_counter = 0
    N = 1000
    for _ = 1:N
        boundary_slots_num = rand(0:20)
        max_nodes_num = rand(1:10)
        max_leg_num = rand(1:8)

        min_nodes_required = max(ceil(Int, boundary_slots_num / max_leg_num), 1)
        expect_error = min_nodes_required > max_nodes_num
        if expect_error
            continue
        end

        d1 = generate_random_tensor_diagram(boundary_slots_num; max_nodes_num=max_nodes_num, max_leg_num=max_leg_num)
        _, no_bubbles = standardize_diagram(d1)

        diagram_is_empty = length(d1.nodes) == 0


        if !no_bubbles
            bubble_counter += 1
        end

        # we are fine with tests fail if there standardization fails 
        # which happens when there is a bubble in the diagram that cannot be reached 
        # from the boundary

        if d1 != d1 && no_bubbles
            reflexivity = false
        end

        d2 = copy(d1)
        if d1 != d2 && no_bubbles
            copy_equality = false
        end

        if !diagram_is_empty
            all_ind = unique(collect(Iterators.flatten(d2.contraction_pattern)))
            idx = rand(all_ind)
            d2.labels[idx] = "?"
            if d1 != d2 && no_bubbles
                unknown_label_equality = false
            end

            d3 = copy(d2)
            d3.nodes[1] = hor_reflection(d3.nodes[1])
            if d1 == d3 && no_bubbles
                structure_inequality = false
            end

            d3 = copy(d2)
            idx = rand(all_ind)
            d3.labels[idx] = "x"
            if d2 == d3 && no_bubbles
                label_inequality = false
            end

            inner_ind = all_ind[all_ind.>0]
            boundary_ind = all_ind[all_ind.<0]
            inner_ind_shuffled = shuffle(inner_ind)
            boundary_ind_shuffled = shuffle(boundary_ind)
            all_ind = vcat(inner_ind, boundary_ind)
            all_ind_shuffled = vcat(inner_ind_shuffled, boundary_ind_shuffled)
            relabeing_dict = Dict{Int,Int}(ind_before => ind_after for (ind_before, ind_after) in zip(all_ind, all_ind_shuffled))

            d2 = relabel(d1, relabeing_dict)

            relabeling_test = d1 == d2


            if !relabeling_test && no_bubbles
                relabeling_equality = false
            end

            if relabeling_test
                perm = randperm(length(d2.nodes))
                d2.nodes = d2.nodes[perm]
                d2.contraction_pattern = d2.contraction_pattern[perm]

                if d1 != d2 && no_bubbles
                    permutation_equality = false
                end
            end
        end

        d2 = copy(d1)
        success = add_random_dangling_indices!(d2; max_num_dangling=3)
        if success
            if d1 == d2
                dangling_indices_break_equality = false
            end
        end
    end

    @test reflexivity
    @test copy_equality
    @test unknown_label_equality
    @test structure_inequality
    @test label_inequality
    @test relabeling_equality
    @test permutation_equality
    @test dangling_indices_break_equality

    println("\nTest encountered $bubble_counter of bubbles in $N experiments\n")

end

###################################################################
# topic 3: REFLECTIONS
###################################################################


@testset "NODE REFLECTIONS                " begin
    # Tests that horizontal reflection is an involution (f(f(x)) == x).
    involution_h = true
    # Tests that vertical reflection is an involution (f(f(x)) == x).
    involution_v = true
    # Checks that legs vary correctly under horizontal reflection.
    legs_permutation_h = true
    # Checks that legs vary correctly under vertical reflection.
    legs_permutation_v = true

    # Checks that reflection with incorrect direction throws an error.
    reflect_throws_for_incorrect_dir = true

    N = 1000
    for _ = 1:N
        necessary_labels = rand(INFINITE_LABELS ∪ FINITE_LABELS, rand(0:2))

        T1 = generate_random_tensor_node(rand(1:20), "T"; max_num_of_allowed_rnd_labels=rand(4:1000), necessary_labels=necessary_labels)

        add_random_allowed_label_combinations!(T1; max_num_allowed_combinations=rand(3:10))
        add_random_forbidden_label_combinations!(T1; max_num_forbidden_combinations=rand(3:10))

        if T1 != hor_reflection(hor_reflection(T1))
            involution_h = false
        end
        if T1 != ver_reflection(ver_reflection(T1))
            involution_v = false
        end

        T1_r = reflect(T1; dir="horizontal")

        for leg_name in keys(T1.legs)
            if leg_name in ["top", "bottom"]
                if T1_r.legs[leg_name] != reverse(T1.legs[leg_name])
                    legs_permutation_h = false
                end
            elseif leg_name == "left"
                if T1_r.legs["right"] != T1.legs[leg_name]
                    legs_permutation_h = false
                end
            elseif leg_name == "right"
                if T1_r.legs["left"] != T1.legs[leg_name]
                    legs_permutation_h = false
                end
            end
        end

        T1_r2 = reflect(T1; dir="vertical")

        for leg_name in keys(T1.legs)
            if leg_name in ["left", "right"]
                if T1_r2.legs[leg_name] != reverse(T1.legs[leg_name])
                    legs_permutation_v = false
                end
            elseif leg_name == "top"
                if T1_r2.legs["bottom"] != T1.legs[leg_name]
                    legs_permutation_v = false
                end
            elseif leg_name == "bottom"
                if T1_r2.legs["top"] != T1.legs[leg_name]
                    legs_permutation_v = false
                end
            end
        end

        try
            reflect(T1; dir="diagonal")
            reflect_throws_for_incorrect_dir = false
        catch
            # Expected error caught
        end
    end

    @test involution_h
    @test involution_v
    @test legs_permutation_h
    @test legs_permutation_v
    @test reflect_throws_for_incorrect_dir
end


@testset "COMPUTING M MATRICES                " begin end