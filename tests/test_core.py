import unittest

import numpy as np

from biology.cell import Cell
from core.mechanics import mechanical_forces_indexed, mechanical_forces_reference, update_positions
from core.scheduler import advance_cell_cycles
from core.spatial_index import SpatialIndex
from core.simulation import Simulation, SimulationConfig
from core.world import World
from biology.regulatory_toggle import (
    analyze_fixed_points, signal_pulse_experiment, signal_sequence_experiment,
)
from analysis.experiments import describe


def cell(identifier, x, y=0.0, radius=1.0, fate=1):
    return Cell(identifier, None, np.array([x, y], dtype=float), radius,
                0.0, 0.0, np.zeros(2), fate=fate)


class SchedulerTests(unittest.TestCase):
    def test_age_and_cycle_continue_at_capacity(self):
        a, b = cell(1, 0), cell(2, 50)
        world = World(cells=[a, b])
        ready = advance_cell_cycles(world, 0.5, 1.0, 2)
        self.assertEqual(ready, [])
        self.assertEqual((a.age, b.age), (0.5, 0.5))
        self.assertEqual((a.cycle_progress, b.cycle_progress), (0.5, 0.5))

    def test_disabled_proliferation_ages_without_cycle_progress(self):
        a = cell(1, 0)
        a.cycle_progress = 0.4
        self.assertEqual(advance_cell_cycles(World(cells=[a]), 0.25, 1.0, 5, False), [])
        self.assertEqual(a.age, 0.25)
        self.assertEqual(a.cycle_progress, 0.4)

    def test_capacity_returns_only_available_divisions(self):
        cells = [cell(i, i * 10.0) for i in range(3)]
        for item in cells:
            item.cycle_progress = 1.0
        ready = advance_cell_cycles(World(cells=cells), 0.1, 1.0, 4)
        self.assertEqual(ready, [cells[0]])
        self.assertTrue(all(item.cycle_progress == 1.0 for item in cells))


class FieldTests(unittest.TestCase):
    def test_no_flux_diffusion_conserves_integrated_mass(self):
        world = World(size=10.0, grid_size=11)
        rng = np.random.default_rng(2)
        world.morphogen[:] = rng.random(world.morphogen.shape)
        before = world.morphogen.sum() * world.spacing**2
        world.update_morphogen(0.01, 0.2, 0.0, 0.0)
        after = world.morphogen.sum() * world.spacing**2
        self.assertAlmostEqual(before, after, places=12)
        self.assertTrue(np.all(world.morphogen >= 0))

    def test_decay_and_bilinear_secretion_mass(self):
        world = World(size=10.0, grid_size=11)
        world.morphogen.fill(1.0)
        world.update_morphogen(0.1, 0.0, 0.2, 0.0)
        self.assertAlmostEqual(world.morphogen.mean(), 0.98)
        world.morphogen.fill(0.0)
        world.cells.append(Cell(1, None, np.array([0.13, -0.27]), 1.0,
                                0.0, 0.0, np.zeros(2), organizer=True))
        row = world.update_morphogen(0.1, 0.0, 0.0, 2.0)
        self.assertAlmostEqual(world.morphogen.sum() * world.spacing**2, 0.2)
        self.assertAlmostEqual(row["secretion_added"], 0.2)
        self.assertEqual(row["clipping_removed"], 0.0)

        for grid_size in (10, 20):
            edge_world = World(size=10.0, grid_size=grid_size)
            edge_world.cells.append(Cell(1, None, np.array([4.999, -4.999]), 1.0,
                                         0.0, 0.0, np.zeros(2), organizer=True))
            edge_world.update_morphogen(0.1, 0.0, 0.0, 2.0)
            self.assertAlmostEqual(edge_world.morphogen.sum() * edge_world.spacing**2, 0.2)

    def test_explicit_instability_is_rejected(self):
        with self.assertRaises(ValueError):
            World(size=10, grid_size=11).update_morphogen(1, 10, 0, 0)

    def test_disabled_field_clearing_is_accounted_as_removal(self):
        world = World(size=10, grid_size=10)
        world.morphogen.fill(0.5)
        row = world.update_morphogen(0.1, 0.2, 0.1, 1.0, enabled=False)
        self.assertEqual(row["final_mass"], 0.0)
        self.assertAlmostEqual(row["field_reset_removed"], row["initial_mass"])
        self.assertEqual(row["uptake_removed"], 0.0)
        self.assertEqual(row["residual"], 0.0)

    def test_diffusion_grid_refinement_reduces_analytic_mode_error(self):
        errors = []
        length, diffusion, duration = 10.0, 0.2, 0.5
        for count in (20, 40, 80):
            world = World(size=length, grid_size=count)
            dx = world.spacing
            x = -length / 2 + (np.arange(count) + 0.5) * dx
            mode = np.cos(np.pi * (x + length / 2) / length)
            world.morphogen[:] = (1.0 + mode)[None, :]
            dt = 0.15 * dx**2 / diffusion
            steps = int(np.ceil(duration / dt))
            dt = duration / steps
            for _ in range(steps):
                world.update_morphogen(dt, diffusion, 0.0, 0.0, record=False)
            expected = 1.0 + mode * np.exp(-diffusion * (np.pi / length) ** 2 * duration)
            errors.append(float(np.sqrt(np.mean((world.morphogen.mean(axis=0) - expected) ** 2))))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_diffusion_time_refinement_reduces_euler_error(self):
        count, length, diffusion, duration = 32, 10.0, 0.2, 0.5
        world = World(size=length, grid_size=count)
        dx = world.spacing
        x = -length / 2 + (np.arange(count) + 0.5) * dx
        mode = np.cos(np.pi * (x + length / 2) / length)
        eigenvalue = -4 * np.sin(np.pi / (2 * count)) ** 2 / dx**2
        exact = 1.0 + mode * np.exp(diffusion * eigenvalue * duration)
        errors = []
        for fraction in (0.24, 0.12, 0.06):
            world.morphogen[:] = (1.0 + mode)[None, :]
            candidate_dt = fraction * dx**2 / diffusion
            steps = int(np.ceil(duration / candidate_dt))
            dt = duration / steps
            for _ in range(steps):
                world.update_morphogen(dt, diffusion, 0.0, 0.0, record=False)
            errors.append(float(np.sqrt(np.mean((world.morphogen.mean(axis=0) - exact) ** 2))))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])


class MechanicsTests(unittest.TestCase):
    def test_pair_force_cutoff_overlap_and_coincident(self):
        cells = [cell(0, 0), cell(1, 2.5), cell(2, 20), cell(3, 0)]
        forces = mechanical_forces_reference(cells)
        self.assertTrue(np.all(np.isfinite(forces)))
        self.assertTrue(np.allclose(forces.sum(axis=0), 0.0))
        self.assertGreater(np.linalg.norm(forces[0]), 0.0)
        self.assertTrue(np.allclose(forces[2], 0.0))

    def test_spatial_hash_matches_all_pairs(self):
        rng = np.random.default_rng(7)
        cells = [cell(i, *rng.uniform(-10, 10, 2), radius=rng.uniform(0.5, 1.5))
                 for i in range(120)]
        cells.extend([cell(120, 1, 1), cell(121, 1, 1)])
        ref = mechanical_forces_reference(cells)
        fast = mechanical_forces_indexed(cells)
        np.testing.assert_allclose(fast, ref, rtol=0, atol=1e-12)

    def test_indexed_neighbors_match_brute_force_counts(self):
        rng = np.random.default_rng(9)
        cells = [cell(i, *rng.uniform(-20, 20, 2), radius=rng.uniform(0.5, 2.0))
                 for i in range(150)]
        index = SpatialIndex(cells, 5.0)
        for i, current in enumerate(cells):
            actual = sum(
                j != i and np.linalg.norm(cells[j].position - current.position)
                < current.radius + cells[j].radius + 1.0
                for j in index.query(current.position, current.radius + max(c.radius for c in cells) + 1.0)
            )
            expected = sum(
                j != i and np.linalg.norm(other.position - current.position)
                < current.radius + other.radius + 1.0
                for j, other in enumerate(cells)
            )
            self.assertEqual(actual, expected)

    def test_hard_wall_keeps_cell_centers_inside_radius_adjusted_domain(self):
        item = cell(0, 20.0, -20.0, radius=1.0)
        world = World(size=10.0, grid_size=10, cells=[item])
        update_positions(world, 0.1, 0.1, 0.035, chemotaxis=False)
        np.testing.assert_allclose(item.position, [4.0, -4.0])


class RecordingTests(unittest.TestCase):
    def test_recording_modes_do_not_change_trajectory(self):
        summaries = []
        for mode in ("minimal", "standard", "trace"):
            sim = Simulation(4, SimulationConfig(max_cells=24, recording_mode=mode,
                                                  recording_interval=3))
            sim.run(20)
            summaries.append((sim.summary(), sim.events,
                              [(c.id, c.position.tolist(), c.gene_expression.tolist())
                               for c in sim.world.cells]))
            if mode == "trace":
                self.assertTrue(sim.world.mass_history)
                self.assertTrue(sim.world.cells[0].history)
            else:
                self.assertFalse(sim.world.mass_history)
                self.assertFalse(sim.world.cells[0].history)
        self.assertEqual(summaries[0], summaries[1])
        self.assertEqual(summaries[1], summaries[2])

    def test_experiment_summary_preserves_categorical_fates(self):
        simulation = Simulation(2, SimulationConfig(max_cells=1, recording_mode="minimal"))
        simulation.run(1)
        result = describe(simulation)
        self.assertNotIn("fate_signal_correlation", result)
        self.assertEqual(result["fate_proportions"]["Stem"], 1.0)
        self.assertIn("mean", result["mean_signal_by_fate"]["Stem"])


class RegulatoryToggleTests(unittest.TestCase):
    def test_default_parameters_have_two_stable_states_and_a_saddle(self):
        points = analyze_fixed_points()
        self.assertEqual(len(points), 3)
        self.assertEqual(sum(point["stable"] for point in points), 2)
        self.assertTrue(all(point["residual_norm"] < 1e-5 for point in points))

    def test_signal_pulse_can_leave_a_maintained_induced_state(self):
        result = signal_pulse_experiment()
        after = np.array(result["after_withdrawal"])
        self.assertGreater(after[0], 1.5)
        self.assertLess(after[1], 0.3)

    def test_initial_conditions_select_distinct_unforced_attractors(self):
        low_a = signal_pulse_experiment(duration=100, pulse_duration=0,
                                        initial=(2.0, 0.05))
        low_b = signal_pulse_experiment(duration=100, pulse_duration=0,
                                        initial=(0.05, 2.0))
        self.assertGreater(low_a["after_withdrawal"][0], low_a["after_withdrawal"][1])
        self.assertGreater(low_b["after_withdrawal"][1], low_b["after_withdrawal"][0])

    def test_induction_strength_and_integration_step_are_explicit(self):
        mild = signal_pulse_experiment(signal=0.2, pulse_duration=2.0)
        strong = signal_pulse_experiment(signal=1.0, pulse_duration=2.0)
        self.assertGreater(strong["after_pulse"][0], mild["after_pulse"][0])
        coarse = signal_pulse_experiment(dt=0.05)
        fine = signal_pulse_experiment(dt=0.01)
        np.testing.assert_allclose(coarse["after_withdrawal"],
                                   fine["after_withdrawal"], atol=2e-4, rtol=0)

    def test_regulatory_state_persists_through_repeated_signal_changes(self):
        result = signal_sequence_experiment(
            (0.05, 2.0), [(1.0, 10.0), (0.0, 30.0), (0.2, 5.0), (0.0, 30.0)]
        )
        states = np.asarray(result["states"])
        self.assertGreater(states[1, 0], states[1, 1])
        self.assertGreater(states[-1, 0], 1.5)
        self.assertLess(states[-1, 1], 0.3)


if __name__ == "__main__":
    unittest.main()
