# Embryo Sim

**An experimental framework for computational developmental biology and artificial morphogenesis.**

Embryo Sim explores how cellular behaviors and local interactions can give rise to multicellular organization. Its long-term direction is to build an extensible, testable simulation framework for processes such as cell proliferation, intercellular signaling, regulatory state changes, differentiation, and tissue mechanics—without prescribing a final anatomical structure.

The project is **not currently a biologically validated embryo simulation**. It does not reproduce human embryonic development, infer a complete developmental program from DNA, or predict real cell-fate outcomes. Current rules and parameters are simplified and illustrative. The name describes the research direction, not an achieved biological fidelity.

## Principles

- **Mechanisms before morphology.** Define local interactions and measure what patterns arise; do not encode a desired organ or tissue shape as the result.
- **Scientific model != numerical solver.** Keep biological assumptions, mathematical equations, integration methods, and experiment definitions distinct.
- **Verification != validation.** Numerical consistency does not alone establish biological accuracy.
- **Traceability.** Track ancestry, regulatory state, fate assignments, parameter sets, and notable events.
- **Reproducibility and convergence.** Record seeds, integration settings, boundary conditions, and versioned experiments; test resolution sensitivity.
- **Explicit uncertainty.** Label biological assumptions, simplifications, and mechanisms not yet implemented.

## Current prototype

The current 2D model initializes a founder cell and includes:

- Discrete cells with positions, growth and division, and parent-daughter lineage.
- A diffusing, degrading extracellular signaling field and a designated signal source.
- Simplified two-factor regulatory activity and rule-based fate transitions.
- Inherited terminal fate labels (Type A and Type C) and a transitional Type B label.
- Versioned finite-range center-based mechanics with an all-pairs reference implementation; **cell-type-dependent adhesion is not implemented**.
- Configurable minimal, standard, and trace recording, event and lineage records, mass accounting, ablations, and an SVG snapshot.

**Interpretation caveat:** Type A and Type C are abstract labels, not verified biological cell types. The current terminal labels can persist after the underlying gene activities change. Commitment is presently a rule for label inheritance, **not a demonstrated self-maintaining regulatory mechanism**. Type B can transition to Type A. These choices are provisional and explicitly under investigation.

Simulation time uses arbitrary model units. `dt = 0.08` advances model time by 0.08 per step; it does not represent minutes or hours. Values such as cell radius and diffusion coefficient likewise should not be interpreted as calibrated physical measurements without an explicit unit and calibration model.

## Quick start

Requires Python 3.12 and NumPy in a configured environment.

```bash
.venv/bin/python main.py --steps 240 --seed 7 --output output/development.svg
```

The command prints a JSON summary and writes an SVG snapshot. The visual encodings are: blue haze for signal concentration, cell colors for fate labels, and a pale ring for the signal source. `--steps`, `--seed`, and `--output` control run length, reproducibility seed, and output path. `-v`, `-vv`, and `-vvv` provide progressively detailed diagnostics on stderr.

Recording defaults to `standard`: periodic aggregate metrics, events, and cell snapshots. Choose `--recording-mode minimal` for aggregate metrics and events or `--recording-mode trace` for complete per-step cell and field diagnostics. `--trace-cell` and `-vvv` automatically use trace recording.

To trace a cell's ancestry and per-step state:

```bash
.venv/bin/python main.py --steps 240 --seed 7 --trace-cell 6
```

To run mechanism ablations, seed comparisons and a timestep comparison:

```bash
.venv/bin/python -m analysis.experiments
```

`biology/` contains the cellular rules; `core/` contains the simulation loop, environment, mechanics, and scheduler; `analysis/` contains experiments and diagnostics; `visualization/` handles outputs.

See [model specification](docs/model-specification.md), [numerical methods](docs/numerical-methods.md), [validation record](docs/validation.md), and [architecture](docs/architecture.md) for equations, assumptions, verification, and extension guidance.

## What existing experiments do and do not demonstrate

- **Demonstrated:** The implementation produces division, lineage inheritance, a diffusing signal field, abstract fate changes, and multicellular spatial arrangements.
- **Numerically checked:** A discrete morphogen mass ledger closes to floating-point precision in reported runs. This is a conservation/accounting check, not independent verification of every PDE discretization or model assumption.
- **Mechanism tests:** Disabling signaling removes Type A in the tested setting; disabling adhesion changes compaction and proliferation; disabling chemotaxis has little effect on the tested endpoint.
- **Not yet established:** Biologically realistic identity maintenance, differential cell sorting, parameter convergence, physiological timescales, validated embryonic structures, or organism-level development.

## Current numerical and model status

The chemical field uses a cell-centered grid, explicit reaction-diffusion update, and no-flux outer boundary. Source deposition conserves its integrated amount under grid refinement. Mechanics is `finite-range-v2`: overlap repulsion and finite-range attraction replace the historical unbounded attraction. This model revision changes trajectories and is not expected to reproduce previous outcomes. A uniform spatial hash accelerates mechanical and crowding neighbors; performance measurements and verification limits are recorded in `docs/validation.md`.

An isolated toggle-network experiment tests feedback and signal-pulse persistence. It is not coupled to default cells, whose regulatory activity remains a simple signal-response model. A numerical example of bistability is not evidence of biological commitment.

## Next scientific and engineering milestones

1. Quantify parameter sensitivity and identify which model outcomes are robust across plausible parameter ranges.
2. Expand timestep and grid convergence reports beyond isolated field modes to coupled runs and chosen observables.
3. Validate analysis measures on synthetic arrangements and add larger seed ensembles with uncertainty summaries.
4. Record model-version metadata in each experiment artifact and add reproducible checkpointing for long runs.
5. Consider new biological mechanisms only after defining an isolated hypothesis, measurable target, and appropriate empirical comparison.

## Scope and limitations

This is currently an exploratory mechanistic model, **not** an in-silico human embryo, organogenesis model, digital twin, clinical instrument, or tool for predicting real embryonic outcomes. Visual similarity to embryonic morphology is not evidence of biological validity. Every future mechanism should have a defined mathematical model, assumptions, units or declared dimensionless conventions, tests, and (when biological fidelity is claimed) appropriate empirical references.

A useful result may be either an emergent pattern **or a documented failure to produce one** under a specified model. The project's goal is to understand the consequences of the rules, not to guarantee an embryo-shaped output.
