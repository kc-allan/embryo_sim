# Mission

## Purpose

Capability checklist for a trustworthy, extensible computational morphogenesis framework.

## Now

- [x] Cell age advances independently of division eligibility and population capacity.
- [x] Cell-cycle progression and computational capacity have explicit semantics.
- [x] Chemical boundaries and source deposition match documented mathematical assumptions.
- [x] Chemical mass accounting uses a defined spatial integration measure.
- [x] Fate and signal statistics preserve categorical meaning.
- [x] Cell mechanics have a documented finite interaction range and reference solver.
- [x] Shared spatial queries accelerate mechanics and crowding without missing neighbors.
- [x] Recording modes preserve reproducibility while making full traces optional.
- [x] Regulatory state, fate commitment, phenotype, and lineage have distinct definitions.
- [x] Numerical methods and experiment evidence are documented and reproducible.

## Next

- [ ] Cell division partitioning is validated against explicit intracellular inheritance assumptions.
- [ ] Cell growth and volume changes are represented independently of cell-cycle timing.
- [ ] Cell death and removal rules can be studied as isolated hypotheses.
- [ ] Multiple signaling species can be coupled through declared reaction terms.
- [ ] Uptake and secretion kinetics can be configured per cell state.
- [ ] Cell shape and polarity can be represented when a supported hypothesis requires them.
- [ ] Mechanical parameters have sensitivity and identifiability analyses.
- [ ] Experimental observations can be imported for parameter calibration.
- [ ] Model versions and run metadata support archival reproduction.
- [ ] Quantitative morphology metrics are validated against synthetic ground truths.

## Later

- [ ] Three-dimensional domains support verified boundary and neighborhood operations.
- [ ] Tissue topology can change through explicit rearrangement rules.
- [ ] Multiscale coupling can combine cellular and tissue-scale equations.
- [ ] Alternative cell-based numerical formulations can be compared on shared experiments.
- [ ] Spatially varying environments can be specified without encoding target anatomy.
- [ ] Parameter uncertainty can be propagated across ensembles of simulations.
- [ ] Experimental perturbations can be represented and compared systematically.
- [ ] Large runs can checkpoint and resume without changing stochastic trajectories.
- [ ] Analysis outputs can be exported in documented, interoperable formats.
- [ ] Biological validation claims are tied to suitable empirical datasets and review.

## Quality Bars

- [x] Numerical properties are tested separately from biological plausibility.
- [x] Optimized algorithms are compared with simple reference implementations.
- [x] Random seeds and parameter sets are recorded for reproducibility.
- [ ] Performance claims include hardware-independent scaling and measured timings.
- [ ] Model changes never target a desired visual morphology without a stated hypothesis.

## Scope Boundaries

- No claim of complete or clinically predictive human embryogenesis.
- No anatomy-generation guarantee or hardcoded target morphology.
- No new biological mechanism without an explicit hypothesis and isolated experiment.
