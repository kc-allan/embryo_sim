# Architecture

## Runtime flow

`core/simulation.py` owns configuration, seeded randomness, event creation,
recording policy, and update ordering. It composes independent modules rather
than defining an intended tissue geometry.

- `biology/cell.py` defines cell state and abstract fate labels.
- `biology/gene_regulation.py` contains the default signal-response model.
- `biology/differentiation.py` contains rule-based fate assignment.
- `biology/regulatory_toggle.py` is an isolated alternative-network
  experiment; it is not part of default multicellular runs.
- `core/world.py` owns the chemical grid, interpolation, source deposition,
  diffusion, and mass diagnostics.
- `core/spatial_index.py` provides deterministic local cell queries.
- `core/mechanics.py` contains the force law, all-pairs reference, indexed
  force solver, and historical force update.
- `core/scheduler.py` handles chronological aging, crowding, and cycle
  progression independently of simulation capacity.
- `analysis/experiments.py` defines ablations and summaries; `analysis/benchmarks.py`
  measures component runtime.
- `visualization/renderer.py` creates a static SVG of a completed state.

## Shared spatial index

The uniform spatial hash is rebuilt from the living-cell list when mechanics
or cycle scheduling runs. A caller supplies a conservative query radius and
performs exact geometric filtering. This keeps neighbor discovery semantics
shared while allowing mechanics and crowding to use different exact distance
rules. Results are sorted by stable input index. Any movement, division, or
removal requires rebuilding; no index persists across state changes.

## Introducing a biological hypothesis

Add a mechanism as an isolated, documented model first. State its equations,
units or dimensionless conventions, assumptions, parameters, and numerical
method. Add unit tests against analytic behavior or a trusted reference. Keep
the existing model available as a baseline, then expose the alternative
through explicit configuration or a standalone experiment. Couple it into
the multicellular loop only after its independent behavior is understood.

Keep four concerns distinguishable in code and reports: biological
hypothesis, mathematical model, numerical implementation, and simulation
infrastructure. A successful output need not resemble a predefined anatomy;
it must be interpretable and reproducible.
