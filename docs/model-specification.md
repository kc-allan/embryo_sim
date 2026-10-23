# Model specification

## Scope and units

The current implementation is a simplified 2D center-based model. It is an
experimental mechanism sandbox, not a calibrated embryo model. Position,
radius, force, chemical concentration, secretion, and time use declared model
units. `dt` is model time per update; it has no conversion to minutes or hours.

## Cell state and identity

Each `Cell` stores position, radius, chronological age, dimensionless cycle
progress, two regulatory activities, a fate label, parent ID, and organizer
status. These quantities have distinct meanings:

- **Regulatory state:** the current values in `gene_expression`.
- **Fate label:** the current abstract category Stem/Type A/Type B/Type C.
- **Commitment:** in the default model, Type A and Type C labels persist and
  are inherited at division; this is a rule for label persistence, not a
  self-maintaining molecular mechanism. Type B can be reassigned.
- **Phenotype:** no independent phenotype model is implemented.
- **Lineage:** parent-child identifiers and division events.

The simple default expression rule relaxes activities toward logistic targets
of the local morphogen signal:

\[
\dot g_k = \alpha(T_k(C)-g_k),\qquad
T_A=\sigma(18(C-0.62)),\quad T_B=\sigma(-18(C-0.28)),
\]

with `alpha=1.8`, Euler integration, and activities clipped to [0,1]. This is
a signal-response rule; it is not a demonstrated bistable gene network.
Fate rules are threshold-based and live in `biology/differentiation.py`.

## Cell cycle and division

Chronological age increases by `dt` for every living cell, including at the
capacity limit and with proliferation disabled. Disabled proliferation pauses
cycle progress. Otherwise, the configured division period is the uncrowded
time to add one normalized cycle unit. At each update, progress advances by

\[
\Delta q = \frac{dt}{T}\max(0.12,1-\min(n/8,0.85)),
\]

where `n` counts neighbors within `r_i+r_j+1` model-length units. Progress is
capped at 1 while waiting for capacity. The capacity is a computational cap,
not an environmental carrying-capacity law. Division resets the parent's
progress by subtracting one and initializes the daughter at zero. Daughter
position is offset from the parent; regulatory activities receive independent
small Gaussian perturbations and are clipped to [0,1]. Fate and lineage are
inherited. Cell age and cycle progress are reset.

## Chemical field

The morphogen follows an explicit reaction-diffusion approximation:

\[
\partial_t C = D\nabla^2 C - \lambda C + S.
\]

`World` uses a cell-centered square grid, spacing `size/grid_size`, and edge
replication for a no-flux boundary. The concentration field is nonnegative;
there is no upper concentration clipping. Organizer secretion is bilinearly
distributed over the surrounding grid cells. The `secretion` parameter is
integrated amount per model-time unit; each deposited contribution is divided
by the grid-cell area so its total is grid-resolution independent. The default
value preserves the previous default grid's total source amount per step.

The mass integral is `sum(C) * spacing**2`. Its ledger reports diffusion,
degradation, secretion, actual uptake (currently zero), signaling-ablation
field removal, and residual. A closed ledger is bookkeeping,
not proof of biological validity.

## Mechanics model revision

Current mechanics version: `finite-range-v2`. For a pair with center distance
`d`, contact distance `R=r_i+r_j`, overlap `o=max(R-d,0)`, and gap
`g=max(d-R,0)`, the force on cell i toward j is

\[
\mathbf F_{ij}=\left[-k_r o + \mathbf 1_{0<g<\ell}
k_a\sin(\pi g/\ell)\right]\hat{\mathbf r}_{ij}.
\]

This gives linear short-range overlap repulsion and smooth attraction only
between contact and finite range `ell`; pair forces are zero outside `R+ell`.
At coincident centers, an ID-ordered x-axis direction regularizes the
otherwise undefined unit vector. Movement uses overdamped Euler updates
`x += mobility * dt * sum(F)`. Cells are projected to remain inside hard square
walls by clamping their centers to the radius-adjusted domain. Parameters are
normalized model quantities, not calibrated physical constants. The previous
unbounded-attraction rule remains available as `update_positions_legacy` for
historical comparison; it is not the current baseline.

## Separate regulatory toggle experiment

`biology/regulatory_toggle.py` is not coupled to the multicellular simulation.
It implements

\[
\dot A=\frac{\beta_A(1+g u)}{1+B^n}-\gamma A,\qquad
\dot B=\frac{\beta_B}{1+A^n}-\gamma B.
\]

The defaults `(beta_A,beta_B,gamma,n,g)=(2,2,1,4,4)` are dimensionless,
illustrative parameters. Numerical basin relaxation and a finite-difference
Jacobian find two stable unforced states and one saddle. A signal pulse from
the B-high initial state moves the system to the A-high state, which remains
after withdrawal for this parameter set. The sequence experiment can apply
repeated signal pulses and withdrawals to inspect state persistence. This is
a model demonstration, not evidence that real cell identities behave this
way.
