import numpy as np

from core.world import World


def advance_cell_cycles(world: World, dt: float, division_period: float, capacity: int) -> list:
    """Advance cell cycles and return ready-to-divide cells in stable order."""
    living = [cell for cell in world.cells if cell.alive]
    if len(living) >= capacity:
        return []
    ready = []
    for cell in living:
        cell.age += dt
        nearby = sum(
            1 for other in living
            if other.id != cell.id
            and np.linalg.norm(other.position - cell.position) < cell.radius + other.radius + 1.0
        )
        crowding = min(nearby / 8.0, 0.85)
        cell.cycle_progress += dt * max(0.12, 1.0 - crowding) / division_period
        if cell.cycle_progress >= 1.0:
            ready.append(cell)
            if len(living) + len(ready) >= capacity:
                break
    return ready
