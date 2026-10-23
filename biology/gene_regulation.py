from typing import TypeAlias, cast

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.float64]


def update_expression(expression: FloatArray, morphogen: float, dt: float) -> FloatArray:
    """A small, illustrative two-gene network with smooth signal response."""
    signal = np.clip(morphogen, 0.0, 1.0)
    target = np.array(
        [1.0 / (1.0 + np.exp(-18.0 * (signal - 0.62))),
         1.0 / (1.0 + np.exp(18.0 * (signal - 0.28)))],
        dtype=float,
    )
    return cast(FloatArray, np.clip(expression + (target - expression) * min(dt * 1.8, 1.0), 0.0, 1.0))
