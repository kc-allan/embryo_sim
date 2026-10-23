from dataclasses import dataclass

import numpy as np


@dataclass
class Cell:
    id: int
    parent_id: int | None
    position: np.ndarray
    radius: float
    age: float
    cycle_progress: float
    gene_expression: np.ndarray
    fate: int = 0
    alive: bool = True
    organizer: bool = False


FATE_NAMES = ("Stem", "Type A", "Type B", "Type C")
FATE_COLORS = ("#e8b84a", "#f16f63", "#63b6a5", "#7489d8")
