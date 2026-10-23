"""Run controlled mechanism and robustness comparisons for the default model."""

from dataclasses import replace
import json

import numpy as np

from biology.cell import FATE_NAMES
from core.simulation import Simulation, SimulationConfig


def describe(simulation: Simulation) -> dict:
    cells = [cell for cell in simulation.world.cells if cell.alive]
    counts = {name: 0 for name in FATE_NAMES}
    centroids = {}
    signal_by_fate = {}
    positions = np.array([cell.position for cell in cells])
    center = positions.mean(axis=0) if len(positions) else np.zeros(2)
    for cell in cells:
        name = FATE_NAMES[cell.fate]
        counts[name] += 1
    for fate, name in enumerate(FATE_NAMES):
        members = [cell for cell in cells if cell.fate == fate]
        if members:
            centroids[name] = np.mean([cell.position for cell in members], axis=0).tolist()
            signal_by_fate[name] = round(float(np.mean([
                simulation.world.sample_morphogen(cell.position) for cell in members
            ])), 5)

    centroid_points = np.array(list(centroids.values()))
    centroid_separation = None
    if len(centroid_points) > 1:
        pairwise = centroid_points[:, None, :] - centroid_points[None, :, :]
        distances = np.linalg.norm(pairwise, axis=2)
        centroid_separation = round(float(distances[np.triu_indices(len(distances), 1)].mean()), 5)

    fate_values = np.array([cell.fate for cell in cells], dtype=float)
    local_signals = np.array([
        simulation.world.sample_morphogen(cell.position) for cell in cells
    ])
    correlation = None
    if len(cells) > 1 and np.std(fate_values) > 0 and np.std(local_signals) > 0:
        correlation = round(float(np.corrcoef(fate_values, local_signals)[0, 1]), 5)

    mass_history = simulation.world.mass_history
    return {
        "seed": simulation.seed,
        "steps": simulation.tick,
        "time": round(simulation.tick * simulation.config.dt, 4),
        "dt": simulation.config.dt,
        "cells": len(cells),
        "fates": counts,
        "divisions": sum(event["type"] == "division" for event in simulation.events),
        "fate_transitions": sum(event["type"] == "fate_transition" for event in simulation.events),
        "mean_distance_from_colony_center": round(float(np.linalg.norm(positions - center, axis=1).mean()), 5) if len(cells) else 0,
        "mean_fate_centroid_separation": centroid_separation,
        "fate_centroids": centroids,
        "mean_signal_by_fate": signal_by_fate,
        "fate_signal_correlation": correlation,
        "morphogen_mass": round(float(
            simulation.world.morphogen.sum() * simulation.world.spacing**2
        ), 6),
        "mass_budget": {
            key: round(sum(float(row[key]) for row in mass_history), 8)
            for key in (
                "diffusion_delta",
                "degradation_delta",
                "secretion_added",
                "uptake_removed",
                "clipping_removed",
                "boundary_removed",
            )
        },
        "max_mass_balance_residual": max(
            (abs(float(row["residual"])) for row in simulation.world.mass_history),
            default=0.0,
        ),
    }


def run_case(name: str, config: SimulationConfig, seed: int, duration: float) -> dict:
    steps = round(duration / config.dt)
    simulation = Simulation(seed=seed, config=config)
    simulation.run(steps)
    return {"case": name, **describe(simulation)}


def main() -> None:
    baseline = SimulationConfig()
    duration = 240 * baseline.dt
    cases = [
        run_case("baseline", baseline, 7, duration),
        run_case("signaling_off", replace(baseline, signaling_enabled=False), 7, duration),
        run_case("adhesion_off", replace(baseline, adhesion=0.0), 7, duration),
        run_case("chemotaxis_off", replace(baseline, chemotaxis_enabled=False), 7, duration),
        run_case("proliferation_off", replace(baseline, proliferation_enabled=False), 7, duration),
        run_case("seed_1", baseline, 1, duration),
        run_case("seed_17", baseline, 17, duration),
        run_case("half_dt", replace(baseline, dt=baseline.dt / 2), 7, duration),
    ]
    print(json.dumps(cases, indent=2))


if __name__ == "__main__":
    main()
