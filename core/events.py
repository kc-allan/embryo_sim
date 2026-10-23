"""Typed records for observable simulation events."""

from typing import Literal, TypeAlias, TypedDict


class DivisionEvent(TypedDict):
    type: Literal["division"]
    tick: int
    time: float
    parent_id: int
    daughter_id: int
    inherited_fate: str
    parent_age: float
    parent_position: list[float]
    daughter_position: list[float]
    parent_gene_expression: list[float]
    daughter_gene_expression: list[float]


FateTransitionEvent = TypedDict(
    "FateTransitionEvent",
    {
        "type": Literal["fate_transition"],
        "tick": int,
        "time": float,
        "cell_id": int,
        "from": str,
        "to": str,
        "age": float,
        "position": list[float],
        "signal": float,
        "gene_expression_before": list[float],
        "gene_expression_after": list[float],
        "parent_id": int | None,
    },
)

SimulationEvent: TypeAlias = DivisionEvent | FateTransitionEvent
