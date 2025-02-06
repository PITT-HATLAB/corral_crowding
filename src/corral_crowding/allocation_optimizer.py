# import networkx
import lovelyplots
import networkx as nx
import numpy as np
import rustworkx as rx
import scienceplots
from rustworkx.visualization import graphviz_draw, mpl_draw
from scipy.optimize import minimize
from scipy.stats import gmean
from tqdm import tqdm

from corral_crowding.detuning_fit import compute_infidelity_parameters, decay_fit
from corral_crowding.module_graph import QuantumModuleGraph


class GateFidelityOptimizer:
    def __init__(
        self,
        module,
        lambdaq,
        eta,
        g3,
        alpha=0.12,  # 120Mhz
        min_bare_space_ghz=0.2,  # 200mhz
        qubit_bounds=(3.3, 5.7),
        snail_bounds=(4.2, 4.7),
        drop_k=0,  # 0 for best, 1 to drop worst, 2 to drop 2 worst, etc
    ):
        self.lambdaq = lambdaq
        self.eta = eta
        self.g3 = g3
        self.alpha = alpha
        self.min_bare_space_ghz = min_bare_space_ghz
        self.qubit_bounds = qubit_bounds
        self.snail_bounds = snail_bounds
        self.module_graph = module
        self.best_frequencies = None
        self.best_cost = np.inf
        self.drop_k = drop_k

        detuning_list = np.linspace(50, 1000, 64)
        self.infidelity_params, _ = compute_infidelity_parameters(
            detuning_list, lambdaq=lambdaq, eta=eta, alpha=120e6, g3=g3
        )

    def _unit_crosstalk(self, intended_freq, spectator_key, spectator_freq):
        distance = np.abs(intended_freq - spectator_freq)
        units_distance = distance * 1e3  # Convert GHz → MHz

        # Strong crosstalk penalty
        if distance < 0.05:
            return 0.5

        # anharmonicity higher level transitions penalty
        if spectator_key == "qubit-qubit" and distance < self.alpha:
            return 0.25

        params = self.infidelity_params.get(spectator_key)
        if params is None:
            raise KeyError(f"Unknown interaction type: {spectator_key}")

        return decay_fit(units_distance, *params)

    def _compute_gate_infidelity(self, edge, interaction_data):
        driven_freq = interaction_data["qubit-qubit"][edge]
        gate_infidelities = sum(
            self._unit_crosstalk(driven_freq, interaction_type, spectator_freq)
            for interaction_type in ["qubit-qubit", "snail-qubit", "qubit-sub"]
            for spectator_edge, spectator_freq in interaction_data[
                interaction_type
            ].items()
            if spectator_edge != edge
        )
        return gate_infidelities

    def _compute_bare_infidelity(self, edge, interaction_data):
        # for each qubit-resonance, add penality of 1 if closer than min_bare_space_ghz to
        # any of the other bare resonances (qubit-resonance or snail-resonance)
        driven_freq = interaction_data["qubit-resonance"][edge]
        cost = lambda freq: 1.0 - freq / self.min_bare_space_ghz
        gate_infidelity = sum(
            (
                cost(np.abs(driven_freq - spec_freq))
                if np.abs(driven_freq - spec_freq) < self.min_bare_space_ghz
                else 0
            )
            for spectator_edge, spec_freq in interaction_data["qubit-resonance"].items()
            if spectator_edge != edge
        )
        gate_infidelity += sum(
            (
                cost(np.abs(driven_freq - spec_freq))
                if np.abs(driven_freq - spec_freq) < self.min_bare_space_ghz
                else 0
            )
            for spec_freq in interaction_data["snail-resonance"].values()
        )
        return gate_infidelity

    def compute_total_infidelity(self, frequencies):
        qubit_frequencies, snail_frequency = frequencies[:-1], frequencies[-1]
        interaction_data = self.module_graph.get_interaction_frequencies(
            qubit_frequencies, snail_frequency
        )
        two_qubit_crowding = [
            self._compute_gate_infidelity(edge, interaction_data)
            for edge in interaction_data["qubit-qubit"]
        ]
        two_qubit_crowding = sum(two_qubit_crowding[self.drop_k :])
        one_qubit_crowding = sum(
            self._compute_bare_infidelity(edge, interaction_data)
            for edge in interaction_data["qubit-resonance"]
        )
        return two_qubit_crowding + one_qubit_crowding

    def optimize_frequencies(self, attempts=64):
        qubit_count = self.module_graph.num_qubits
        self.best_cost = np.inf
        for _ in tqdm(range(attempts)):
            initial_guess = np.append(
                np.random.uniform(
                    self.qubit_bounds[0], self.qubit_bounds[1], qubit_count
                ),
                np.random.uniform(self.snail_bounds[0], self.snail_bounds[1]),
            )
            result = minimize(
                self.compute_total_infidelity,
                initial_guess,
                bounds=[self.qubit_bounds] * qubit_count + [self.snail_bounds],
                method="Nelder-Mead",
            )
            # print(result)
            if result.fun < self.best_cost:
                self.best_cost = result.fun
                self.best_frequencies = result.x
                best_result = result
        print(best_result.message)
        return self.best_frequencies, self.best_cost

    def get_final_average(self):
        if self.best_frequencies is None:
            print("No optimized frequencies available.")
            return

        qubit_frequencies, snail_frequency = (
            self.best_frequencies[:-1],
            self.best_frequencies[-1],
        )
        interaction_data = self.module_graph.get_interaction_frequencies(
            qubit_frequencies, snail_frequency
        )
        gate_infidelities = {
            edge: self._compute_gate_infidelity(edge, interaction_data)
            for edge in list(interaction_data["qubit-qubit"])[self.drop_k :]
        }

        avg_gate_infidelity = gmean(list(gate_infidelities.values()))
        return avg_gate_infidelity

    def report_results(self):
        if self.best_frequencies is None:
            print("No optimized frequencies available.")
            return

        qubit_frequencies, snail_frequency = (
            self.best_frequencies[:-1],
            self.best_frequencies[-1],
        )
        interaction_data = self.module_graph.get_interaction_frequencies(
            qubit_frequencies, snail_frequency
        )

        print("Qubit Frequencies:", qubit_frequencies, "GHz")
        print(f"SNAIL Frequency: {snail_frequency} GHz")
        print("Gate Infidelities:")
        gate_infidelities = {
            edge: self._compute_gate_infidelity(edge, interaction_data)
            for edge in list(interaction_data["qubit-qubit"])[self.drop_k :]
        }
        for edge, infidelity in gate_infidelities.items():
            freq = interaction_data["qubit-qubit"][edge]
            print(f"  Gate {edge}:{freq:.6f} GHz → Infidelity: {infidelity:.6f}")

        avg_gate_infidelity = gmean(list(gate_infidelities.values()))
        print(f"\nAverage Gate Infidelity: {avg_gate_infidelity:.6f}")

        self.module_graph.plot_graph(qubit_frequencies, snail_frequency)
        self.module_graph.plot_interaction_frequencies(
            qubit_frequencies, snail_frequency
        )
