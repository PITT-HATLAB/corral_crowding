import networkx as nx

# !pip install mqt.bench
from qiskit import transpile
from qiskit.transpiler import CouplingMap

# from rustworkx.visualization import graphviz_draw, mpl_draw

# from mqt.bench import CompilerSettings, QiskitSettings, TKETSettings, get_benchmark


def build_graphs(snails, qubits, edges):
    # Create the snail-qubit graph
    snail_qubit_graph = nx.Graph()
    node_to_index = {
        node: idx for idx, node in enumerate(snails + qubits)
    }  # Map node values to graph indices
    snail_qubit_graph.add_nodes_from(snails + qubits)  # Add all nodes
    snail_qubit_graph.add_edges_from(
        [(node_to_index[u], node_to_index[v], 0) for u, v in edges]
    )

    # Create the qubit connectivity graph
    qubit_connectivity = nx.Graph()
    qubit_to_index = {
        qubit: idx for idx, qubit in enumerate(qubits)
    }  # Map qubit values to graph indices
    qubit_connectivity.add_nodes_from(qubits)  # Add only qubit nodes

    # Map each snail to its connected qubits
    snail_to_qubits = {snail: [] for snail in snails}
    for u, v in edges:
        if u in snails and v in qubits:
            snail_to_qubits[u].append(v)
        elif v in snails and u in qubits:
            snail_to_qubits[v].append(u)

    # Add edges between qubits sharing the same snail
    for connected_qubits in snail_to_qubits.values():
        for i, q1 in enumerate(connected_qubits):
            for q2 in connected_qubits[i + 1 :]:
                # Add edge between the indices of the qubits
                if not qubit_connectivity.has_edge(
                    qubit_to_index[q1], qubit_to_index[q2]
                ):
                    qubit_connectivity.add_edge(
                        qubit_to_index[q1], qubit_to_index[q2], True
                    )
    return snail_qubit_graph, qubit_connectivity


########################################################################
# 2-qubit module, ring topology with 16 qubits
snails_ring = [i for i in range(1, 33, 2)]  # Odd indices for snails
qubits_ring = [i for i in range(0, 32, 2)]  # Even indices for qubits
edges_ring = [(qubits_ring[i], snails_ring[i]) for i in range(16)] + [
    (snails_ring[i], qubits_ring[(i + 1) % 16]) for i in range(16)
]  # Ring connections

ring = [snails_ring, qubits_ring, edges_ring]


########################################################################
# 2-qubit module, square-lattice topology with 16 qubits
# Define snails and qubits
qubits_square = [0, 2, 4, 6, 11, 13, 15, 17, 22, 24, 26, 28, 33, 35, 37, 39]
snails_square = [
    1,
    3,
    5,
    7,
    8,
    9,
    10,
    12,
    14,
    16,
    18,
    19,
    20,
    21,
    23,
    25,
    27,
    29,
    30,
    31,
    32,
    34,
    36,
    38,
]
edges_square = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (0, 7),
    (2, 8),
    (4, 9),
    (6, 10),
    (7, 11),
    (8, 13),
    (9, 15),
    (10, 17),
    (11, 12),
    (12, 13),
    (13, 14),
    (14, 15),
    (15, 16),
    (16, 17),
    (11, 18),
    (13, 19),
    (15, 20),
    (17, 21),
    (18, 22),
    (19, 24),
    (20, 26),
    (21, 28),
    (22, 23),
    (23, 24),
    (24, 25),
    (25, 26),
    (26, 27),
    (27, 28),
    (22, 29),
    (24, 30),
    (26, 31),
    (28, 32),
    (29, 33),
    (30, 35),
    (31, 37),
    (32, 39),
    (33, 34),
    (34, 35),
    (35, 36),
    (36, 37),
    (37, 38),
    (38, 39),
]

# Combine the data for the lattice
square = [snails_square, qubits_square, edges_square]


########################################################################
# Define snails and qubits for the double ring
qubits_tworing = [1, 3, 5, 7, 9, 11, 13, 15, 16, 18, 20, 22, 24, 26, 28, 30]
snails_tworing = [0, 2, 4, 6, 8, 10, 12, 14, 17, 19, 21, 23, 25, 27, 29, 31]
edges_tworing = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 8),
    (8, 9),
    (9, 10),
    (10, 11),
    (11, 12),
    (12, 13),
    (13, 14),
    (14, 15),
    (15, 0),
    (16, 17),
    (17, 18),
    (18, 19),
    (19, 20),
    (20, 21),
    (21, 22),
    (22, 23),
    (23, 24),
    (24, 25),
    (25, 26),
    (26, 27),
    (27, 28),
    (28, 29),
    (29, 30),
    (30, 31),
    (31, 16),
    (0, 16),
    (1, 17),
    (2, 18),
    (3, 19),
    (4, 20),
    (5, 21),
    (6, 22),
    (7, 23),
    (8, 24),
    (9, 25),
    (10, 26),
    (11, 27),
    (12, 28),
    (13, 29),
    (14, 30),
    (15, 31),
]

# Combine the data for the double ring
tworing = [snails_tworing, qubits_tworing, edges_tworing]

########################################################################
qubits_hex = [0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27]
snail_hex = [4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23]
edges_hex = [
    (0, 4),
    (0, 5),
    (1, 5),
    (1, 6),
    (2, 6),
    (2, 7),
    (3, 7),
    (4, 8),
    (5, 9),
    (6, 10),
    (7, 11),
    (8, 12),
    (9, 12),
    (9, 13),
    (10, 13),
    (10, 14),
    (11, 14),
    (11, 15),
    (12, 16),
    (13, 17),
    (14, 18),
    (15, 19),
    (16, 20),
    (16, 21),
    (17, 21),
    (17, 22),
    (18, 22),
    (18, 23),
    (19, 23),
    (20, 24),
    (21, 25),
    (22, 26),
    (23, 27),
]

hex_topo = [snail_hex, qubits_hex, edges_hex]

########################################################################
qubits_corral = [0, 1, 2, 3, 4, 5, 6, 7, 15, 16, 17, 18, 19, 20, 21, 22]
snail_corral = [8, 9, 10, 11, 12, 13, 14, 23]
edges_corral = [
    (0, 8),
    (1, 8),
    (1, 9),
    (2, 9),
    (2, 10),
    (3, 10),
    (3, 11),
    (4, 11),
    (4, 12),
    (5, 12),
    (5, 13),
    (6, 13),
    (6, 14),
    (7, 14),
    (8, 15),
    (8, 16),
    (9, 16),
    (9, 17),
    (10, 17),
    (10, 18),
    (11, 18),
    (11, 19),
    (12, 19),
    (12, 20),
    (13, 20),
    (13, 21),
    (14, 21),
    (14, 22),
    (23, 7),
    (23, 22),
    (23, 0),
    (23, 15),
]

corral = [snail_corral, qubits_corral, edges_corral]

########################################################################
qubits_denselattice = [
    0,
    2,
    4,
    # 6,
    8,
    10,
    12,
    # 14,
    15,
    17,
    19,
    21,
    23,
    25,
    27,
    # 29,
    # 30,
    32,
    34,
    36,
]
snail_denselattice = [1, 3, 5, 7, 9, 11, 13, 16, 18, 20, 22, 24, 26, 28, 31, 33, 35]
edges_denselattice = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    # (5, 6),
    (7, 8),
    (8, 9),
    (9, 10),
    (10, 11),
    (11, 12),
    (12, 13),
    # (13, 14),
    (15, 16),
    (16, 17),
    (17, 18),
    (18, 19),
    (19, 20),
    (20, 21),
    (22, 23),
    (23, 24),
    (24, 25),
    (25, 26),
    (26, 27),
    (27, 28),
    # (28, 29),
    # (30, 31),
    (31, 32),
    (32, 33),
    (33, 34),
    (34, 35),
    (35, 36),
    (0, 7),
    (1, 8),
    (2, 9),
    (3, 10),
    (4, 11),
    (5, 12),
    # (6, 13),
    (7, 15),
    (8, 16),
    (9, 17),
    (10, 18),
    (11, 19),
    (12, 20),
    (13, 21),
    (15, 22),
    (16, 23),
    (17, 24),
    (18, 25),
    (19, 26),
    (20, 27),
    (21, 28),
    # (22, 30),
    (23, 31),
    (24, 32),
    (25, 33),
    (26, 34),
    (27, 35),
    (28, 36),
]
denselattice = [snail_denselattice, qubits_denselattice, edges_denselattice]

best_snails = [0, 1, 2, 3, 4, 5, 6, 7]
best_qubits = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
best_edges = [
    (1, 21),
    (1, 22),
    (1, 23),
    (1, 24),
    (5, 18),
    (5, 19),
    (5, 20),
    (5, 24),
    (6, 15),
    (6, 16),
    (6, 17),
    (6, 23),
    (7, 12),
    (7, 13),
    (7, 14),
    (7, 22),
    (8, 9),
    (8, 10),
    (8, 11),
    (8, 21),
    (2, 11),
    (2, 14),
    (2, 17),
    (2, 20),
    (3, 10),
    (3, 13),
    (3, 16),
    (3, 19),
    (4, 9),
    (4, 12),
    (4, 15),
    (4, 18),
]
best = [best_snails, best_qubits, best_edges]
