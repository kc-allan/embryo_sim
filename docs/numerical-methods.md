# Numerical methods

## Update order

One call to `Simulation.step()` performs: (1) morphogen diffusion, degradation,
and organizer secretion; (2) default regulatory activity and fate updates;
(3) chronological age and cycle progression, followed by divisions;
(4) mechanics and chemotaxis; (5) aggregate metrics and configured recordings.
New daughter cells therefore do not update regulatory state until the next
tick. Division randomness uses the simulation's seeded NumPy generator.

## Reaction-diffusion field

The field is cell-centered, square, with `dx = world_size/grid_size`. The
five-point Laplacian uses edge-value replication, representing zero normal
flux at the boundary. The explicit Euler update is

\[
C^{n+1}=C^n+dt(DL_h C^n-\lambda C^n)+dt S^n.
\]

For a nonnegative field with nonnegative source, the implementation enforces
the conservative positivity bound
`dt * (4*D/dx**2 + decay) <= 1`. This is a sufficient stability/positivity
condition for this split linear operator, not a convergence guarantee. Pure
diffusion conserves `sum(C)*dx**2` to floating-point precision. Degradation is
first-order explicit and therefore does not exactly match continuous-time
exponential decay at finite `dt`. The source is applied after field transport
using bilinear weights whose sum is one, including clamped behavior near the
domain edge. There is no concentration saturation or clipping. Disabling
signaling clears the existing field as an experiment intervention; the
removed amount is recorded separately from biological uptake, which is not
implemented.

## Mechanics and spatial search

Mechanics is overdamped, with no inertial state. The all-pairs reference and
spatial-hash implementation use the same pair-force equation. The hash bucket
width is the largest possible pair cutoff (twice the largest cell radius plus
adhesion range); each query visits all buckets intersecting its search radius,
then exact Euclidean distances are evaluated. Pairs are ordered and unique.
The index is rebuilt from current positions each update, so divisions and
movement cannot leave stale entries. For roughly uniform populations and
bounded local density, work follows the local candidate count; dense clusters
can still require quadratic work.

## Recording

`minimal` records periodic aggregate metrics and all biological events and
lineage IDs. `standard` additionally stores periodic cell snapshots. `trace`
stores full per-cell per-step histories and per-step mass diagnostics. The
recording choice does not consume random values or alter update rules. `-vvv`
and `--trace-cell` select trace mode in the CLI.

## Numerical limitations

The position grid interpolation and cell mechanics are not a coupled
continuum solver. Source deposition, edge handling, and grid conventions are
explicitly specified, but no physical calibration is available. Stability
tests answer whether an implementation remains bounded under the stated
conditions; convergence tests answer whether a numerical answer approaches a
reference as resolution changes. Neither establishes biological validity.
