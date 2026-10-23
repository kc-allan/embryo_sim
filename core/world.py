from dataclasses import dataclass, field

import numpy as np

from biology.cell import Cell


@dataclass
class World:
    size: float = 80.0
    grid_size: int = 128
    cells: list[Cell] = field(default_factory=list)
    morphogen: np.ndarray = field(init=False)

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

    def update_morphogen(self, dt: float, diffusion: float, decay: float, secretion: float) -> None:
        padded = np.pad(self.morphogen, 1, mode="edge")
        laplacian = (
            padded[1:-1, :-2] + padded[1:-1, 2:]
            + padded[:-2, 1:-1] + padded[2:, 1:-1]
            - 4.0 * self.morphogen
        ) / self.spacing**2
        self.morphogen += dt * (diffusion * laplacian - decay * self.morphogen)
        for cell in self.cells:
            if cell.alive and cell.organizer:
                x, y = self._grid_position(cell.position)
                self.morphogen[int(round(y)), int(round(x))] += secretion * dt
        np.clip(self.morphogen, 0.0, 1.0, out=self.morphogen)
        self.morphogen[0, :] = 0.0
        self.morphogen[-1, :] = 0.0
        self.morphogen[:, 0] = 0.0
        self.morphogen[:, -1] = 0.0
