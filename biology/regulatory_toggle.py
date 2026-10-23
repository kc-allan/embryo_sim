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

import numpy as np


@dataclass(frozen=True)
class ToggleParameters:
    beta_a: float = 2.0
    beta_b: float = 2.0
    gamma: float = 1.0
    hill_n: int = 4
    input_gain: float = 4.0

    def __post_init__(self):
        if min(self.beta_a, self.beta_b, self.gamma, self.input_gain) < 0:
            raise ValueError("production, decay, and input parameters must be nonnegative")
        if self.gamma == 0 or self.hill_n < 1:
            raise ValueError("gamma and hill exponent must be positive")


def derivative(state, signal: float, parameters: ToggleParameters = ToggleParameters()):
    if not 0.0 <= signal <= 1.0:
        raise ValueError("signal must be between 0 and 1")
    a, b = np.maximum(np.asarray(state, dtype=float), 0.0)
    n = parameters.hill_n
    da = parameters.beta_a * (1 + parameters.input_gain * signal) / (1 + b**n) - parameters.gamma * a
    db = parameters.beta_b / (1 + a**n) - parameters.gamma * b
    return np.array([da, db])


def integrate(state, duration: float, dt: float, signal: float,
              parameters: ToggleParameters = ToggleParameters()) -> np.ndarray:
    """Integrate with fixed-step RK4, shortening the last step to duration."""
    if duration < 0 or dt <= 0:
        raise ValueError("duration must be nonnegative and dt positive")
    result = np.maximum(np.asarray(state, dtype=float).copy(), 0.0)
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


def equilibria(initial_states=None, duration=100.0, dt=0.01,
               signal=0.0, parameters=ToggleParameters(), tolerance=1e-5):
    """Find attractors by relaxation from distinct initial states.

    This numerical basin survey supplements, but does not replace, solving
    fixed-point equations and evaluating the Jacobian at each equilibrium.
    """
    if initial_states is None:
        initial_states = ([0.05, 0.05], [2.0, 0.05], [0.05, 2.0], [2.0, 2.0])
    points = []
    for initial in initial_states:
        point = integrate(initial, duration, dt, signal, parameters)
        if not any(np.linalg.norm(point - prior) < tolerance for prior in points):
            points.append(point)
    return sorted(points, key=lambda point: (point[0], point[1]))


def analyze_fixed_points(signal=0.0, parameters=ToggleParameters()):
    """Return relaxed fixed points, residuals, and local Jacobian eigenvalues."""
    results = []
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


def signal_pulse_experiment(duration=40.0, pulse_duration=10.0, dt=0.01,
                            signal=1.0, initial=(0.05, 2.0),
                            parameters=ToggleParameters()):
    """Apply a transient input then withdraw it, returning phase endpoints."""
    induced = integrate(initial, pulse_duration, dt, signal, parameters)
    final = integrate(induced, max(0.0, duration - pulse_duration), dt, 0.0, parameters)
    return {"initial": list(initial), "after_pulse": induced.tolist(),
            "after_withdrawal": final.tolist(), "signal": signal,
            "pulse_duration": pulse_duration, "duration": duration,
            "dt": dt, "parameters": parameters.__dict__.copy()}


def signal_sequence_experiment(initial, phases, dt=0.01,
                               parameters=ToggleParameters()):
    """Return the regulatory state after each ``(signal, duration)`` phase."""
    state = np.asarray(initial, dtype=float)
    states = [state.tolist()]
    for signal, duration in phases:
        state = integrate(state, duration, dt, signal, parameters)
        states.append(state.tolist())
    return {"states": states, "phases": [list(phase) for phase in phases],
            "dt": dt, "parameters": parameters.__dict__.copy()}
