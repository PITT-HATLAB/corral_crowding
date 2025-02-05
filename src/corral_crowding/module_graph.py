import matplotlib.pyplot as plt
import networkx as nx


class QuantumModuleGraph:
    def __init__(self, qubit_frequencies, snail_frequency):
        """Initializes a quantum module graph with only qubit-qubit and snail-qubit interactions.

        Args:
            qubit_frequencies (list): List of qubit bare frequencies in GHz.
            snail_frequency (float): SNAIL bare frequency in GHz.
        """
        self.G = nx.Graph()
        self.qubit_frequencies = qubit_frequencies
        self.snail_frequency = snail_frequency

        # Build the graph structure
        self._add_nodes()
        self._add_edges()

    def _add_nodes(self):
        """Adds qubits and the SNAIL as nodes in the graph."""
        for i, freq in enumerate(self.qubit_frequencies):
            self.G.add_node(f"Q{i}", type="qubit", bare_frequency=freq, color="green")

        self.G.add_node(
            "SNAIL", type="snail", bare_frequency=self.snail_frequency, color="red"
        )

    def _add_edges(self):
        """Adds only qubit-qubit and snail-qubit interactions as edges."""
        num_qubits = len(self.qubit_frequencies)

        # Qubit-qubit interactions (blue edges)
        for i in range(num_qubits):
            for j in range(i + 1, num_qubits):
                self.G.add_edge(
                    f"Q{i}", f"Q{j}", interaction="qubit-qubit", color="blue"
                )

        # Qubit-SNAIL interactions (orange edges)
        for i in range(num_qubits):
            self.G.add_edge(f"Q{i}", "SNAIL", interaction="snail-qubit", color="orange")

    def set_frequencies(self, qubit_frequencies, snail_frequency):
        """Sets the bare frequencies of qubits and the SNAIL."""
        self.qubit_frequencies = qubit_frequencies
        self.snail_frequency = snail_frequency

        # Update the bare frequencies of qubits and the SNAIL
        for i, freq in enumerate(self.qubit_frequencies):
            self.G.nodes[f"Q{i}"]["bare_frequency"] = freq

        self.G.nodes["SNAIL"]["bare_frequency"] = snail_frequency

    def get_interaction_frequencies(self):
        interaction_freqs = {
            "qubit-resonance": {},  # { Q0: freq }
            "snail-resonance": {},  # { SNAIL: freq }
            ################################
            "qubit-qubit": {},  # { (Q0, Q1): freq_diff }
            "snail-qubit": {},  # { (Q2, SNAIL): freq_diff }
            "qubit-sub": {},  # { Q0: freq/2 }
            # Intermodule placeholders for future updates
            "snail-qubit (inter)": {},
            "qubit-sub (inter)": {},
            # placeholder for snail speedlimit
            "snail-sub": {},  # { SNAIL: freq/2 }
        }

        # Store interaction frequencies for qubit-qubit and snail-qubit edges
        for u, v in self.G.edges:
            f_u = self.G.nodes[u]["bare_frequency"]
            f_v = self.G.nodes[v]["bare_frequency"]
            freq_diff = abs(f_u - f_v)

            interaction_type = self.G.edges[u, v]["interaction"]

            if interaction_type in ["qubit-qubit", "snail-qubit"]:
                interaction_freqs[interaction_type][(u, v)] = freq_diff

        # Store individual qubit resonances and subharmonics
        for node in self.G.nodes:
            node_type = self.G.nodes[node]["type"]
            freq = self.G.nodes[node]["bare_frequency"]

            if node_type == "qubit":
                interaction_freqs["qubit-resonance"][(node,)] = freq
                interaction_freqs["qubit-sub"][(node,)] = freq / 2

            elif node_type == "snail":
                interaction_freqs["snail-resonance"][(node,)] = freq
                interaction_freqs["snail-sub"][(node,)] = freq / 2

        return interaction_freqs

    def plot_graph(self):
        """Plots the quantum module graph with specified node and edge colors."""
        pos = nx.spring_layout(self.G, seed=42)  # Layout for better visualization

        # Node colors
        node_colors = [self.G.nodes[node]["color"] for node in self.G.nodes]

        # Edge colors
        edge_colors = [self.G.edges[u, v]["color"] for u, v in self.G.edges]

        # Labels for nodes
        labels = {
            node: f"{node}\n{self.G.nodes[node]['bare_frequency']:.2f} GHz"
            for node in self.G.nodes
        }

        plt.figure(figsize=(5, 5))

        # Draw nodes
        nx.draw(
            self.G,
            pos,
            with_labels=True,
            labels=labels,
            node_size=2000,
            node_color=node_colors,
            edgecolors="black",
        )

        # Draw edges with specific colors
        nx.draw_networkx_edges(self.G, pos, edge_color=edge_colors, width=2)

        # Add edge labels
        edge_labels = {
            (u, v): self.G.edges[u, v]["interaction"] for u, v in self.G.edges
        }
        nx.draw_networkx_edge_labels(self.G, pos, edge_labels=edge_labels, font_size=8)

        plt.title("Quantum Module (SNAIL + Qubits)")
        plt.show()

    def plot_interaction_frequencies(self):
        """Plots the bare qubit frequencies, SNAIL frequency, and calculated interaction frequencies."""
        interaction_freqs = self.get_interaction_frequencies()
        all_freqs = list(interaction_freqs["qubit-resonance"].values()) + list(
            interaction_freqs["snail-resonance"].values()
        )

        with plt.style.context(["ieee", "use_mathtext", "science"]):
            fig, ax = plt.subplots(figsize=(3.5, 1))

            max_freq = max(all_freqs) * 1.1
            ax.set_xlim(0, max_freq)
            ax.get_yaxis().set_visible(False)

            # Prepare to track which labels have been added to avoid duplicates
            added_labels = set()

            # Updated Color Mapping
            color_map = {
                "qubit-qubit": "blue",
                "qubit-resonance": "green",
                "snail-qubit": "orange",
                "qubit-sub": "lightgreen",  # Qubit subharmonic in light green
                "snail-sub": "lightcoral",  # Snail subharmonic in light red
                "snail-resonance": "red",  # Snail resonance in red
            }

            # Plot interactions with color-coding
            for interaction_type, freqs in interaction_freqs.items():
                if not freqs:
                    continue

                color = color_map.get(interaction_type, "black")

                label = (
                    interaction_type.replace("-", " ").title()
                    if interaction_type not in added_labels
                    else ""
                )

                for freq in freqs.values():
                    ax.axvline(
                        freq,
                        color=color,
                        linestyle="-",
                        linewidth=2,
                        alpha=0.8,
                        label=label,
                    )
                    added_labels.add(interaction_type)

            ax.set_xlabel("Frequency (GHz)")
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(
                by_label.values(),
                by_label.keys(),
                loc="upper center",
                bbox_to_anchor=(0.5, -0.25),
                ncol=2,
                fontsize=10,
            )

        plt.show()

    def get_graph(self):
        """Returns the NetworkX graph object."""
        return self.G
