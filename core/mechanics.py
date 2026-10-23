import numpy as np

from core.world import World


def update_positions(
    world: World,
    dt: float,
    mobility: float,
    adhesion: float,
    chemotaxis: bool = True,
) -> None:
    """Resolve soft cell overlap, weak cohesion, and stem-cell chemotaxis."""
    cells = [cell for cell in world.cells if cell.alive]
    if not cells:
        return
    positions = np.stack([cell.position for cell in cells])
    radii = np.array([cell.radius for cell in cells])
    forces = np.zeros_like(positions)
    displacement = positions[None, :, :] - positions[:, None, :]
    distance = np.linalg.norm(displacement, axis=2)
    np.fill_diagonal(distance, np.inf)
    minimum = radii[:, None] + radii[None, :]
    direction = displacement / np.maximum(distance[:, :, None], 1e-9)
    overlap = np.maximum(minimum - distance, 0.0)
    near = np.maximum(distance - minimum, 0.0)
    pair_force = (-overlap * 1.8 + np.minimum(near, 1.2) * adhesion)[:, :, None] * direction
    forces += pair_force.sum(axis=1)
    for index, cell in enumerate(cells):
        if chemotaxis and cell.fate == 0:
            forces[index] += world.morphogen_gradient(cell.position) * 0.12
    next_positions = positions + mobility * dt * forces
    limit = world.size / 2.0 - radii
    next_positions = np.clip(next_positions, -limit[:, None], limit[:, None])
    for cell, position in zip(cells, next_positions):
        cell.position = position
