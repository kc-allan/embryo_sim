"""Versioned, finite-range mechanics for the center-based 2D model."""

import numpy as np

from core.spatial_index import SpatialIndex
from core.world import World

MECHANICS_MODEL_VERSION = "finite-range-v2"


def _pair_direction(a, b, id_a: int, id_b: int) -> np.ndarray:
    delta = b - a
    distance = float(np.linalg.norm(delta))
    if distance > 1e-12:
        return delta / distance, distance
    # Deterministic regularization for a physically underdetermined overlap.
    return (np.array([1.0, 0.0]) if id_a < id_b else np.array([-1.0, 0.0])), 0.0


def _pair_force(cell_i, cell_j, repulsion: float, adhesion: float,
                adhesion_range: float) -> np.ndarray:
    direction, distance = _pair_direction(
        cell_i.position, cell_j.position, cell_i.id, cell_j.id
    )
    contact = cell_i.radius + cell_j.radius
    overlap = max(contact - distance, 0.0)
    gap = max(distance - contact, 0.0)
    scalar = -repulsion * overlap
    if 0.0 < gap < adhesion_range:
        # Smooth attraction that vanishes at contact and at the finite cutoff.
        scalar += adhesion * np.sin(np.pi * gap / adhesion_range)
    return scalar * direction


def mechanical_forces_reference(cells, repulsion=1.8, adhesion=0.035,
                                adhesion_range=1.2) -> np.ndarray:
    """Simple O(N²) reference force calculation, including each pair once."""
    if repulsion < 0 or adhesion < 0 or adhesion_range <= 0:
        raise ValueError("mechanics parameters must be nonnegative and range positive")
    cells = list(cells)
    forces = np.zeros((len(cells), 2), dtype=float)
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            force = _pair_force(cells[i], cells[j], repulsion, adhesion, adhesion_range)
            forces[i] += force
            forces[j] -= force
    return forces


def mechanical_forces_indexed(cells, repulsion=1.8, adhesion=0.035,
                               adhesion_range=1.2) -> np.ndarray:
    """Spatial-hash implementation equivalent to the reference force law."""
    if repulsion < 0 or adhesion < 0 or adhesion_range <= 0:
        raise ValueError("mechanics parameters must be nonnegative and range positive")
    cells = list(cells)
    if not cells:
        return np.zeros((0, 2), dtype=float)
    max_radius = max(cell.radius for cell in cells)
    cutoff = 2.0 * max_radius + adhesion_range
    index = SpatialIndex(cells, max(cutoff, 1e-9))
    forces = np.zeros((len(cells), 2), dtype=float)
    for i, j in index.candidate_pairs(cutoff):
        force = _pair_force(cells[i], cells[j], repulsion, adhesion, adhesion_range)
        forces[i] += force
        forces[j] -= force
    return forces


def update_positions(world: World, dt: float, mobility: float, adhesion: float,
                     chemotaxis: bool = True, repulsion: float = 1.8,
                     adhesion_range: float = 1.2) -> None:
    """Apply overdamped finite-range forces, then project cells at hard walls.

    ``adhesion`` and ``repulsion`` are normalized force-per-distance
    coefficients in the model's arbitrary units. Pair forces vanish beyond
    the sum of the radii and ``adhesion_range``. Exact coincident centers use a
    deterministic x-axis separation convention.
    """
    if dt < 0 or mobility < 0:
        raise ValueError("dt and mobility must be nonnegative")
    cells = [cell for cell in world.cells if cell.alive]
    if not cells:
        return
    forces = mechanical_forces_indexed(cells, repulsion, adhesion, adhesion_range)
    for index, cell in enumerate(cells):
        if chemotaxis and cell.fate == 0:
            forces[index] += world.morphogen_gradient(cell.position) * 0.12
    positions = np.stack([cell.position for cell in cells]) + mobility * dt * forces
    radii = np.array([cell.radius for cell in cells])
    limit = world.size / 2.0 - radii
    positions = np.clip(positions, -limit[:, None], limit[:, None])
    for cell, position in zip(cells, positions):
        cell.position = position


def update_positions_legacy(world: World, dt: float, mobility: float,
                            adhesion: float, chemotaxis: bool = True) -> None:
    """Historical unbounded-attraction law retained for comparisons only."""
    cells = [cell for cell in world.cells if cell.alive]
    if not cells:
        return
    positions = np.stack([cell.position for cell in cells])
    radii = np.array([cell.radius for cell in cells])
    displacement = positions[None, :, :] - positions[:, None, :]
    distance = np.linalg.norm(displacement, axis=2)
    np.fill_diagonal(distance, np.inf)
    minimum = radii[:, None] + radii[None, :]
    direction = displacement / np.maximum(distance[:, :, None], 1e-9)
    overlap = np.maximum(minimum - distance, 0.0)
    near = np.maximum(distance - minimum, 0.0)
    forces = ((-overlap * 1.8 + np.minimum(near, 1.2) * adhesion)[:, :, None]
              * direction).sum(axis=1)
    for index, cell in enumerate(cells):
        if chemotaxis and cell.fate == 0:
            forces[index] += world.morphogen_gradient(cell.position) * 0.12
    next_positions = positions + mobility * dt * forces
    limit = world.size / 2.0 - radii
    next_positions = np.clip(next_positions, -limit[:, None], limit[:, None])
    for cell, position in zip(cells, next_positions):
        cell.position = position
