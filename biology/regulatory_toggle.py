"""Independent two-factor toggle model for regulatory commitment experiments.

The model is an experimental hypothesis, not the default cell rule. In
dimensionless model time its Hill-repression equations are

    dA/dt = (beta_A * (1 + input_gain * u)) / (1 + B**n) - gamma*A
    dB/dt = beta_B / (1 + A**n) - gamma*B

where ``u`` is an externally imposed signal in [0, 1]. Parameters must be
measured or deliberately declared for each experiment; this model has no
biological calibration.
"""

from dataclasses import dataclass
from collections.abc import Sequence
from typing import Any, TypeAlias, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray: TypeAlias = NDArray[np.float64]


@dataclass(frozen=True)
class ToggleParameters:
    beta_a: float = 2.0
    beta_b: float = 2.0
    gamma: float = 1.0
    hill_n: int = 4
    input_gain: float = 4.0

    def __post_init__(self) -> None:
        if min(self.beta_a, self.beta_b, self.gamma, self.input_gain) < 0:
            raise ValueError("production, decay, and input parameters must be nonnegative")
        if self.gamma == 0 or self.hill_n < 1:
            raise ValueError("gamma and hill exponent must be positive")


def derivative(state: ArrayLike, signal: float,
               parameters: ToggleParameters = ToggleParameters()) -> FloatArray:
    if not 0.0 <= signal <= 1.0:
        raise ValueError("signal must be between 0 and 1")
    a, b = np.maximum(np.asarray(state, dtype=float), 0.0)
    n = parameters.hill_n
    da = parameters.beta_a * (1 + parameters.input_gain * signal) / (1 + b**n) - parameters.gamma * a
    db = parameters.beta_b / (1 + a**n) - parameters.gamma * b
    return cast(FloatArray, np.array([da, db], dtype=np.float64))


def integrate(state: ArrayLike, duration: float, dt: float, signal: float,
              parameters: ToggleParameters = ToggleParameters()) -> FloatArray:
    """Integrate with fixed-step RK4, shortening the last step to duration."""
    if duration < 0 or dt <= 0:
        raise ValueError("duration must be nonnegative and dt positive")
    result = cast(FloatArray, np.maximum(np.asarray(state, dtype=np.float64).copy(), 0.0))
    elapsed = 0.0
    while elapsed < duration - 1e-14:
        h = min(dt, duration - elapsed)
        k1 = derivative(result, signal, parameters)
        k2 = derivative(result + h * k1 / 2, signal, parameters)
        k3 = derivative(result + h * k2 / 2, signal, parameters)
        k4 = derivative(result + h * k3, signal, parameters)
        result += h * (k1 + 2*k2 + 2*k3 + k4) / 6
        result = np.maximum(result, 0.0)
        elapsed += h
    return result


def equilibria(initial_states: Sequence[ArrayLike] | None = None,
               duration: float = 100.0, dt: float = 0.01,
               signal: float = 0.0,
               parameters: ToggleParameters = ToggleParameters(),
               tolerance: float = 1e-5) -> list[FloatArray]:
    """Find attractors by relaxation from distinct initial states.

    This numerical basin survey supplements, but does not replace, solving
    fixed-point equations and evaluating the Jacobian at each equilibrium.
    """
    if initial_states is None:
        initial_states = ([0.05, 0.05], [2.0, 0.05], [0.05, 2.0], [2.0, 2.0])
    points: list[FloatArray] = []
    for initial in initial_states:
        point = integrate(initial, duration, dt, signal, parameters)
        if not any(np.linalg.norm(point - prior) < tolerance for prior in points):
            points.append(point)
    return sorted(points, key=lambda point: (point[0], point[1]))


def analyze_fixed_points(signal: float = 0.0,
                         parameters: ToggleParameters = ToggleParameters()) -> list[dict[str, Any]]:
    """Return relaxed fixed points, residuals, and local Jacobian eigenvalues."""
    results: list[dict[str, Any]] = []
    for point in equilibria(signal=signal, parameters=parameters):
        eps = 1e-5
        jacobian = np.column_stack([
            (derivative(point + np.eye(2)[i] * eps, signal, parameters)
             - derivative(point - np.eye(2)[i] * eps, signal, parameters)) / (2 * eps)
            for i in range(2)
        ])
        eigenvalues = np.linalg.eigvals(jacobian)
        results.append({
            "state": point.tolist(),
            "residual_norm": float(np.linalg.norm(derivative(point, signal, parameters))),
            "jacobian": jacobian.tolist(),
            "eigenvalues": [[float(value.real), float(value.imag)] for value in eigenvalues],
            "stable": bool(np.all(eigenvalues.real < 0)),
        })
    return results


def signal_pulse_experiment(duration: float = 40.0, pulse_duration: float = 10.0,
                            dt: float = 0.01, signal: float = 1.0,
                            initial: Sequence[float] = (0.05, 2.0),
                            parameters: ToggleParameters = ToggleParameters()) -> dict[str, Any]:
    """Apply a transient input then withdraw it, returning phase endpoints."""
    induced = integrate(initial, pulse_duration, dt, signal, parameters)
    final = integrate(induced, max(0.0, duration - pulse_duration), dt, 0.0, parameters)
    return {"initial": list(initial), "after_pulse": induced.tolist(),
            "after_withdrawal": final.tolist(), "signal": signal,
            "pulse_duration": pulse_duration, "duration": duration,
            "dt": dt, "parameters": parameters.__dict__.copy()}


def signal_sequence_experiment(initial: Sequence[float],
                               phases: Sequence[tuple[float, float]],
                               dt: float = 0.01,
                               parameters: ToggleParameters = ToggleParameters()) -> dict[str, Any]:
    """Return the regulatory state after each ``(signal, duration)`` phase."""
    state = cast(FloatArray, np.asarray(initial, dtype=np.float64))
    states = [state.tolist()]
    for signal, duration in phases:
        state = integrate(state, duration, dt, signal, parameters)
        states.append(state.tolist())
    return {"states": states, "phases": [list(phase) for phase in phases],
            "dt": dt, "parameters": parameters.__dict__.copy()}
