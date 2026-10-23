from dataclasses import dataclass, field
from typing import Literal, TypeAlias, TypedDict

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.float64]
Fate: TypeAlias = Literal[0, 1, 2, 3]


class CellHistoryEntry(TypedDict):
    tick: int
    time: float
    age: float
    position: list[float]
    fate: str
    signal: float
    gene_expression: list[float]
    cycle_progress: float


@dataclass
class Cell:
    id: int
    parent_id: int | None
    position: FloatArray
    radius: float
    age: float
    cycle_progress: float
    gene_expression: FloatArray
    fate: Fate = 0
    alive: bool = True
    organizer: bool = False
    history: list[CellHistoryEntry] = field(default_factory=list, repr=False)


FATE_NAMES = ("Stem", "Type A", "Type B", "Type C")
FATE_COLORS = ("#e8b84a", "#f16f63", "#63b6a5", "#7489d8")
