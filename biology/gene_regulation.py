import numpy as np


def update_expression(expression: np.ndarray, morphogen: float, dt: float) -> np.ndarray:
    """A small, illustrative two-gene network with smooth signal response."""
    signal = np.clip(morphogen, 0.0, 1.0)
    target = np.array(
        [1.0 / (1.0 + np.exp(-18.0 * (signal - 0.62))),
         1.0 / (1.0 + np.exp(18.0 * (signal - 0.28)))],
        dtype=float,
    )
    return np.clip(expression + (target - expression) * min(dt * 1.8, 1.0), 0.0, 1.0)
