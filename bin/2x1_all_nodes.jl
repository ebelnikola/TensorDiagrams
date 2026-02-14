# space labels
#INFINITE_LABELS:
#  - "x" # infinite dimensional space
#  - "u" # infinite dimensional space
#  - "d" # infinite dimensional space
#  - "r" # infinite dimensional space
#FINITE_LABELS:
#  - "z" # space where A_n lives
#  - "z'" # space into which J_n_zp maps
#  - "q" # the space of the side leg of Q
#  - "g" # the wiggly space of T
#  - "z'"# the part of z⊗z that will go to the infinite dimensional sectors  
#  - "t" # trivial part of auxilary space
#  - "w" # weights part of auxilary space
#  - "c" # cancellation part of auxilary space

##############################
# FINITE DIMENSIONAL TENSORS #
##############################

An_nd = TensorNode(
    name="An",
    color=BLUE,
    legs=Dict("left" => [1], "top" => [2], "right" => [3], "bottom" => [4]),
    allowed_labels=[["z"], ["z"], ["z"], ["z"]],
    width=1.2,
    height=1.2
)
Jn_z_nd = TensorNode(
    name="Jnz",
    color=BLUE,
    legs=Dict("top" => [1], "bottom" => [2, 3]),
    allowed_labels=[["z"], ["z"], ["z"]],
    out_legs_num=1,
    width=3.0,
    height=1.4
)
Jn_zp_nd = TensorNode(
    name="Jnzp",
    color=BLUE,
    legs=Dict("top" => [1], "bottom" => [2, 3]),
    allowed_labels=[["z'"], ["z"], ["z"]],
    out_legs_num=1,
    width=3.0,
    height=1.4
)
Q_nd = TensorNode(
    name="Q",
    color=BLUE,
    legs=Dict("right" => [1], "top" => [2, 3], "bottom" => [4]),
    allowed_labels=[["q"], ["z"], ["z"], ["z"]],
    out_legs_num=1,
    width=1.4,
    height=1.2
)
T_nd = TensorNode(
    name="T",
    color=BLUE,
    legs=Dict("right" => [1], "left" => [2], "top" => [3, 4], "bottom" => [5]),
    allowed_labels=[["g"], ["z"], ["z"], ["z"], ["z"]],
    dual=Bool[false, false, false, false, true],
    out_legs_num=1,
    width=1.4,
    height=1.2
)
U_nd = TensorNode(
    name="U",
    color=BLUE,
    legs=Dict("top" => [1, 2, 3, 4], "bottom" => [5, 6]),
    allowed_labels=[["z"], ["z"], ["z"], ["z"], ["z"], ["z"]],
    out_legs_num=0,
    width=4.0,
    height=1.2
)

################################
# INFINITE DIMENSIONAL TENSORS #
################################

b_nd = TensorNode(
    name="b",
    color=PINK,
    legs=Dict("left" => [1], "top" => [2], "right" => [3], "bottom" => [4]),
    allowed_labels=[["z", "u", "d", "r"], ["z", "x"], ["z", "u", "d", "r"], ["z", "x"]],
    forbidden_label_combinations=[["z", "z", "z", "z"]],
    width=1.2,
    height=1.2
)

b_1x_nd = TensorNode(
    name="b₁ₓ",
    color=PINK,
    legs=Dict("left" => [1], "top" => [2], "right" => [3], "bottom" => [4]),
    allowed_labels=[["z", "x"], ["z", "x"], ["z", "x"], ["z", "x"]],
    allowed_label_combinations=[
        ["x", "z", "z", "z"],
        ["z", "x", "z", "z"],
        ["z", "z", "x", "z"],
        ["z", "z", "z", "x"],
    ],
    width=1.2,
    height=1.2
)

bnc_nd = TensorNode( # this is b without corner contributions, we will need it in some diagrams
    name="bnc",
    color=PINK,
    legs=Dict("left" => [1], "top" => [2], "right" => [3], "bottom" => [4]),
    allowed_labels=[["z", "u", "d", "r"], ["z", "x"], ["z", "u", "d", "r"], ["z", "x"]],
    forbidden_label_combinations=[
        ["z", "z", "z", "z"],
        ["z", "x", "u", "z"],
        ["z", "z", "d", "x"],
        ["u", "x", "z", "z"],
        ["d", "z", "z", "x"]
    ],
    width=1.2,
    height=1.2
)

bc_nd = TensorNode( # this is only corner contributions to b when the left leg of b is in "z". This is used in Bw, see discussion after 4.27
    name="bc",
    color=PINK,
    legs=Dict("left" => [1], "top" => [2], "right" => [3], "bottom" => [4]),
    allowed_labels=[["z", "u", "d", "r"], ["z", "x"], ["z", "u", "d", "r"], ["z", "x"]],
    allowed_label_combinations=[
        ["z", "x", "u", "z"],
        ["z", "z", "d", "x"],
    ],
    width=1.2,
    height=1.2
)

###############
# AUXILIARIES #
###############

aux_width = 1.2

et_nd = TensorNode(
    name="et",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["t"]],
    width=aux_width,
    height=0.6
)

ew1_nd = TensorNode(
    name="ew1",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["w"]],
    width=aux_width,
    height=0.6
)

ew2_nd = TensorNode(
    name="ew2",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["w"]],
    width=aux_width,
    height=0.6
)

# this node will represent axuiliary vector w1 * ew2 + ew3
ew23_nd = TensorNode(
    name="ew23",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["w"]],
    width=aux_width,
    height=0.6
)

# this node will represent axuiliary vector w1 * ew1 + ew3
ew13_nd = TensorNode(
    name="ew13",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["w"]],
    width=aux_width,
    height=0.6
)

ec_nd = TensorNode(
    name="ec",
    color=BLUE,
    legs=Dict("right" => [1]),
    allowed_labels=[["c"]],
    width=aux_width,
    height=0.6
)

all_nodes = [An_nd, Jn_z_nd, Jn_zp_nd, Q_nd, T_nd, U_nd, b_nd, bnc_nd, bc_nd, b_1x_nd, et_nd, ew1_nd, ew2_nd, ew23_nd, ew13_nd, ec_nd]
nodes_consistency = check_node_consistency.(all_nodes) .|> x -> all(values(x))

if !all(nodes_consistency)
    inconsistent_indices = findall(!, nodes_consistency)
    inconsistent_nodes_names = all_nodes[inconsistent_indices] .|> x -> x.name


    @warn "Nodes $(inconsistent_nodes_names) are inconsistent! Please use check_node_consistency to see details."
end

# Plot all nodes
plot_tensor_node_list(all_nodes; cols=5, size=(1600, 900))