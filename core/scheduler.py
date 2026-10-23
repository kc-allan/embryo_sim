"""Cell-cycle progression and capacity-aware division eligibility."""

import numpy as np

from core.spatial_index import SpatialIndex
from core.world import World


def advance_cell_cycles(world: World, dt: float, division_period: float,
                        capacity: int, enabled: bool = True) -> list:
    """Age all living cells; return ready cells that fit current capacity.

    Capacity limits new cells only. At capacity, cycles can reach (and remain
    at) 1.0, indicating division readiness, while chronological age continues.
    Disabling proliferation advances age but pauses cell-cycle progress.
    """
    if dt < 0 or division_period <= 0 or capacity < 0:
        raise ValueError("dt/capacity must be nonnegative and division_period positive")
    living = [cell for cell in world.cells if cell.alive]
    for cell in living:
        cell.age += dt
    if not enabled or not living:
        return []
    slots = max(0, capacity - len(living))
    max_radius = max(cell.radius for cell in living)
    search_radius = 2.0 * max_radius + 1.0
    index = SpatialIndex(living, max(search_radius, 1e-9))
    ready = []
    for i, cell in enumerate(living):
        nearby = 0
        query_radius = cell.radius + max_radius + 1.0
        for j in index.query(cell.position, query_radius):
            other = living[j]
            if j != i and np.linalg.norm(other.position - cell.position) < cell.radius + other.radius + 1.0:
                nearby += 1
        crowding = min(nearby / 8.0, 0.85)
        cell.cycle_progress = min(
            1.0,
            cell.cycle_progress + dt * max(0.12, 1.0 - crowding) / division_period,
        )
        if cell.cycle_progress >= 1.0 and len(ready) < slots:
            ready.append(cell)
    return ready
