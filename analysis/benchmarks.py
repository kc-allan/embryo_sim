"""Fixed-population component benchmarks for local algorithm comparisons."""

import json
import tempfile
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from biology.cell import Cell
from biology.gene_regulation import update_expression
from core.mechanics import mechanical_forces_indexed, mechanical_forces_reference, update_positions
from core.scheduler import advance_cell_cycles
from core.simulation import SimulationConfig
from core.spatial_index import SpatialIndex
from core.world import World
from visualization.renderer import render_svg


def timed(function: Callable[[], Any], repeats: int = 3) -> float:
    function()  # warm-up
    start = time.perf_counter()
    for _ in range(repeats):
        function()
    return (time.perf_counter() - start) / repeats


def benchmark(sizes: Sequence[int] = (100, 300, 1000, 3000),
              seed: int = 123) -> dict[str, Any]:
    report: dict[str, Any] = {"seed": seed, "units": "seconds per operation", "populations": []}
    rng = np.random.default_rng(seed)
    for count in sizes:
        positions = rng.uniform(-400, 400, (count, 2))
        cells = [Cell(i, None, positions[i].copy(), 1.0, 0.0, 0.0,
                      np.zeros(2), fate=1) for i in range(count)]
        world = World(size=1000, grid_size=128, cells=cells)
        mechanics = timed(lambda: mechanical_forces_indexed(cells), repeats=3)
        position_update = timed(lambda: update_positions(world, 0.01, 0.1, 0.035, False), repeats=3)
        reference = None
        if count <= 300:
            reference = timed(lambda: mechanical_forces_reference(cells), repeats=2)
            np.testing.assert_allclose(mechanical_forces_indexed(cells),
                                       mechanical_forces_reference(cells), atol=1e-12, rtol=0)
        index = SpatialIndex(cells, 3.0)
        neighbors = timed(lambda: [index.query(cell.position, 3.0) for cell in cells], repeats=3)
        scheduler = timed(lambda: advance_cell_cycles(world, 0.01, 1.7, count + 1), repeats=3)
        field = timed(lambda: world.update_morphogen(0.01, 0.2, 0.01, 0.0,
                                                      record=False), repeats=3)
        expression = np.zeros(2)
        gene = timed(lambda: update_expression(expression, 0.4, 0.01), repeats=1000)
        event = {"type": "division", "parent_id": 1, "daughter_id": 2,
                 "time": 1.0, "state": [0.1, 0.2]}
        event_recording = timed(lambda: [dict(event) for _ in range(count)], repeats=10)
        payload = {"config": SimulationConfig().__dict__,
                   "cells": [{"id": cell.id, "position": cell.position.tolist()}
                             for cell in cells], "events": [event] * count}
        serialization = timed(lambda: json.dumps(payload), repeats=3)
        with tempfile.TemporaryDirectory() as directory:
            rendering = timed(lambda: render_svg(world, 1, Path(directory) / "bench.svg"), repeats=2)
        report["populations"].append({
            "cells": count, "mechanics_indexed": mechanics,
            "mechanics_position_update": position_update,
            "mechanics_reference": reference, "neighbor_queries": neighbors,
            "cell_cycle_scheduler": scheduler, "chemical_field": field,
            "gene_update_per_cell": gene, "event_recording": event_recording,
            "serialization": serialization, "rendering": rendering,
        })
    return report


if __name__ == "__main__":
    print(json.dumps(benchmark(), indent=2))
