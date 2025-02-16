import matplotlib.pyplot as plt
import networkx as nx


class QuantumModuleGraph:
    def __init__(self, num_qubits):
        self.G = nx.Graph()
        self.num_qubits = num_qubits  # Default for topology setup
        self._add_edges()

    def _add_edges(self):
        for i in range(self.num_qubits):
            for j in range(i + 1, self.num_qubits):
                self.G.add_edge(
                    f"Q{i}", f"Q{j}", interaction="qubit-qubit", color="blue"
                )
            self.G.add_edge(f"Q{i}", "SNAIL", interaction="snail-qubit", color="orange")

    def get_interaction_frequencies(self, qubit_frequencies, snail_frequency):
        interaction_freqs = {
            "qubit-qubit": {},
            "snail-qubit": {},
            "qubit-resonance": {},
            "snail-resonance": {},
            "qubit-sub": {},
            "snail-sub": {},
        }
        for i, freq in enumerate(qubit_frequencies):
            interaction_freqs["qubit-resonance"][f"Q{i}"] = freq
            interaction_freqs["qubit-sub"][f"Q{i}"] = freq / 2
        interaction_freqs["snail-resonance"]["SNAIL"] = snail_frequency
        interaction_freqs["snail-sub"]["SNAIL"] = snail_frequency / 2
        for u, v in self.G.edges:
            if u.startswith("Q") and v.startswith("Q"):
                interaction_freqs["qubit-qubit"][(u, v)] = abs(
                    qubit_frequencies[int(u[1:])] - qubit_frequencies[int(v[1:])]
                )
            elif v == "SNAIL":
                interaction_freqs["snail-qubit"][(u, v)] = abs(
                    qubit_frequencies[int(u[1:])] - snail_frequency
                )
        return interaction_freqs

    def plot_graph(self, qubit_frequencies, snail_frequency):
        pos = nx.spring_layout(self.G, seed=42)
        labels = {
            node: (
                f"{node}\n{qubit_frequencies[int(node[1:])]:.2f} GHz"
                if node.startswith("Q")
                else f"SNAIL\n{snail_frequency:.2f} GHz"
            )
            for node in self.G.nodes
        }
        node_colors = [
            "green" if node.startswith("Q") else "red" for node in self.G.nodes
        ]
        plt.figure(figsize=(2, 2))
        nx.draw(
            self.G,
            pos,
            with_labels=True,
            node_color=node_colors,
            edgecolors="black",
            width=2,
        )
        nx.draw_networkx_edges(
            self.G,
            pos,
            edge_color=[self.G.edges[e]["color"] for e in self.G.edges],
            width=2,
        )
        plt.show()

    def plot_interaction_frequencies(self, qubit_frequencies, snail_frequency):
        all_freqs = list(qubit_frequencies) + [snail_frequency]
        interaction_freqs = self.get_interaction_frequencies(
            qubit_frequencies, snail_frequency
        )
        with plt.style.context(["ieee", "use_mathtext", "science"]):
            fig, ax = plt.subplots(figsize=(3.5, 1))
            max_freq = max(all_freqs) * 1.05
            ax.set_xlim(0, max_freq)
            ax.get_yaxis().set_visible(False)
            added_labels = set()
            color_map = {
                "qubit-qubit": "blue",
                "qubit-resonance": "green",
                "snail-qubit": "orange",
                "qubit-sub": "gray",
                "snail-sub": "magenta",
                "snail-resonance": "red",
            }
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
                        linewidth=1.5,
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
