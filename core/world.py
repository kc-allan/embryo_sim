from dataclasses import dataclass, field

import numpy as np

from biology.cell import Cell


@dataclass
class World:
    size: float = 80.0
    grid_size: int = 128
    cells: list[Cell] = field(default_factory=list)
    morphogen: np.ndarray = field(init=False)
    mass_history: list[dict[str, float | int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.morphogen = np.zeros((self.grid_size, self.grid_size), dtype=float)

    @property
    def spacing(self) -> float:
        return self.size / (self.grid_size - 1)

    def _grid_position(self, position: np.ndarray) -> tuple[float, float]:
        x = (position[0] + self.size / 2.0) / self.spacing
        y = (position[1] + self.size / 2.0) / self.spacing
        return float(np.clip(x, 0, self.grid_size - 1)), float(np.clip(y, 0, self.grid_size - 1))

    def sample_morphogen(self, position: np.ndarray) -> float:
        x, y = self._grid_position(position)
        x0, y0 = int(x), int(y)
        x1, y1 = min(x0 + 1, self.grid_size - 1), min(y0 + 1, self.grid_size - 1)
        fx, fy = x - x0, y - y0
        top = self.morphogen[y0, x0] * (1 - fx) + self.morphogen[y0, x1] * fx
        bottom = self.morphogen[y1, x0] * (1 - fx) + self.morphogen[y1, x1] * fx
        return float(top * (1 - fy) + bottom * fy)

    def morphogen_gradient(self, position: np.ndarray) -> np.ndarray:
        x, y = self._grid_position(position)
        ix, iy = int(round(x)), int(round(y))
        left, right = max(ix - 1, 0), min(ix + 1, self.grid_size - 1)
        up, down = max(iy - 1, 0), min(iy + 1, self.grid_size - 1)
        dx = max(right - left, 1) * self.spacing
        dy = max(down - up, 1) * self.spacing
        return np.array([
            (self.morphogen[iy, right] - self.morphogen[iy, left]) / dx,
            (self.morphogen[down, ix] - self.morphogen[up, ix]) / dy,
        ])

    def update_morphogen(
        self,
        dt: float,
        diffusion: float,
        decay: float,
        secretion: float,
        tick: int = 0,
        enabled: bool = True,
    ) -> dict[str, float | int]:
        voxel_area = self.spacing**2
        initial_mass = float(self.morphogen.sum() * voxel_area)
        if not enabled:
            self.morphogen.fill(0.0)
            diagnostic = {
                "tick": tick,
                "initial_mass": initial_mass,
                "diffusion_delta": 0.0,
                "degradation_delta": 0.0,
                "secretion_added": 0.0,
                "uptake_removed": 0.0,
                "clipping_removed": 0.0,
                "boundary_removed": 0.0,
                "final_mass": 0.0,
                "residual": 0.0,
            }
            self.mass_history.append(diagnostic)
            return diagnostic

        old_field = self.morphogen.copy()
        padded = np.pad(self.morphogen, 1, mode="edge")
        laplacian = (
            padded[1:-1, :-2] + padded[1:-1, 2:]
            + padded[:-2, 1:-1] + padded[2:, 1:-1]
            - 4.0 * self.morphogen
        ) / self.spacing**2
        diffusion_delta = float(dt * diffusion * laplacian.sum() * voxel_area)
        degradation_delta = float(-dt * decay * old_field.sum() * voxel_area)
        self.morphogen += dt * (diffusion * laplacian - decay * self.morphogen)
        secretion_added = 0.0
        for cell in self.cells:
            if cell.alive and cell.organizer:
                x, y = self._grid_position(cell.position)
                self.morphogen[int(round(y)), int(round(x))] += secretion * dt
                secretion_added += secretion * dt * voxel_area
        unclipped = self.morphogen.copy()
        np.clip(self.morphogen, 0.0, 1.0, out=self.morphogen)
        clipping_removed = float((unclipped - self.morphogen).sum() * voxel_area)
        boundary_removed = float(
            self.morphogen[0, :].sum() + self.morphogen[-1, :].sum()
            + self.morphogen[1:-1, 0].sum() + self.morphogen[1:-1, -1].sum()
        ) * voxel_area
        self.morphogen[0, :] = 0.0
        self.morphogen[-1, :] = 0.0
        self.morphogen[:, 0] = 0.0
        self.morphogen[:, -1] = 0.0
        final_mass = float(self.morphogen.sum() * voxel_area)
        expected_mass = (
            initial_mass + diffusion_delta + degradation_delta + secretion_added
            - clipping_removed - boundary_removed
        )
        diagnostic = {
            "tick": tick,
            "initial_mass": initial_mass,
            "diffusion_delta": diffusion_delta,
            "degradation_delta": degradation_delta,
            "secretion_added": float(secretion_added),
            "uptake_removed": 0.0,
            "clipping_removed": clipping_removed,
            "boundary_removed": boundary_removed,
            "final_mass": final_mass,
            "residual": final_mass - expected_mass,
        }
        self.mass_history.append(diagnostic)
        return diagnostic
