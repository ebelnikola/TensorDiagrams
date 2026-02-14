include("src/TensorDiagrams.jl");
load_config("config_dev.yml")
include("bin/2x1_all_nodes.jl")
boundary_legs_3 = Dict("left" => [-1, -2], "right" => [-5, -6, -7, -8], "top" => [-3, -4], "bottom" => [-9, -10])

An_1f_nd = flip(An_nd, 1)
tst = TensorNode(
    name="test node",
    legs=Dict("top" => [1, 2], "bottom" => [3, 4])
)
An_nd |> pretty_print
diag = TensorDiagram(
    nodes=[hor_reflection(An_1f_nd), tst, Q_nd, b_nd, hor_reflection(et_nd), ew1_nd],
    contraction_pattern=[[1, -3, -2, -9], [1, 2, -6, -10], [-7, 3, 4, 2], [3, -4, -8, 4], [-1], [-5]],
    boundary_slots_num=Dict("horizontal" => 6, "vertical" => 3),
    boundary_legs=boundary_legs_3,
    boundary_legs_posidx=Dict("left" => [1, 4], "right" => [1, 4, 5, 6], "top" => [1, 3], "bottom" => [1, 3]),
    labels=Dict(-4 => "x", -8 => "u",),
    node_coordinates=[
        [-2, 0],
        [2, 0],
        [2, 1.8],
        [2, 3.5],
        [-3, -3],
        [3, -3]
    ],
    factor="1/w2"
)

plot_tensor_diagram(diag; scale=1.0, strokewidth=2.0)

plot_tensor_node(
    TensorNode(
        name="test node",
        legs=Dict("top" => [1, 2], "bottom" => [3, 4]),
    ), scale=1)

TensorNode(
    name="test node",
    legs=Dict("top" => [1, 2, 51], "bottom" => [3, 4]),
)

include("src/TensorDiagrams.jl");
plot_tensor_node_list(all_nodes; cols=5, size=(2000, 1600))

