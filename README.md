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
- Simple center-based mechanics with local forces; **cell-type-dependent adhesion is not yet implemented**.
- Optional mechanism ablations, cell histories, population summaries, morphogen mass accounting, and an SVG snapshot.

**Interpretation caveat:** Type A and Type C are abstract labels, not verified biological cell types. The current terminal labels can persist after the underlying gene activities change. Commitment is presently a rule for label inheritance, **not a demonstrated self-maintaining regulatory mechanism**. Type B can transition to Type A. These choices are provisional and explicitly under investigation.

Simulation time uses arbitrary model units. `dt = 0.08` advances model time by 0.08 per step; it does not represent minutes or hours. Values such as cell radius and diffusion coefficient likewise should not be interpreted as calibrated physical measurements without an explicit unit and calibration model.

## Quick start

Requires Python 3.12 and NumPy in a configured environment.

```bash
.venv/bin/python main.py --steps 240 --seed 7 --output output/development.svg
```

The command prints a JSON summary and writes an SVG snapshot. The visual encodings are: blue haze for signal concentration, cell colors for fate labels, and a pale ring for the signal source. `--steps`, `--seed`, and `--output` control run length, reproducibility seed, and output path. `-v`, `-vv`, and `-vvv` provide progressively detailed diagnostics on stderr.

To trace a cell's ancestry and per-step state:

```bash
.venv/bin/python main.py --steps 240 --seed 7 --trace-cell 6
```

To run mechanism ablations, seed comparisons and a timestep comparison:

```bash
.venv/bin/python -m analysis.experiments
```

`biology/` contains the cellular rules; `core/` contains the simulation loop, environment, mechanics, and scheduler; `analysis/` contains experiments and diagnostics; `visualization/` handles outputs.

## What existing experiments do and do not demonstrate

- **Demonstrated:** The implementation produces division, lineage inheritance, a diffusing signal field, abstract fate changes, and multicellular spatial arrangements.
- **Numerically checked:** A discrete morphogen mass ledger closes to floating-point precision in reported runs. This is a conservation/accounting check, not independent verification of every PDE discretization or model assumption.
- **Mechanism tests:** Disabling signaling removes Type A in the tested setting; disabling adhesion changes compaction and proliferation; disabling chemotaxis has little effect on the tested endpoint.
- **Not yet established:** Biologically realistic identity maintenance, differential cell sorting, parameter convergence, physiological timescales, validated embryonic structures, or organism-level development.

## Next scientific and engineering milestones

1. Specify precisely what regulatory state, fate label, developmental commitment, and phenotype each mean; introduce and test a mechanistic commitment model.
2. Test the numerical pieces independently, including diffusion conservation/convergence, regulatory fixed points, and mechanical force equivalence.
3. Replace all-pairs mechanical neighbor enumeration with a spatial index while retaining the reference solver and verifying equivalent forces for identical states.
4. Introduce configurable telemetry modes and benchmark mechanics, fields, biology, and logging separately.
5. Test timestep/grid convergence and robustness across repeated seeds; report quantitative morphology rather than relying on a snapshot.
6. Add cell-type-dependent mechanics or signaling feedback only with an explicit hypothesis and measurable validation target.

## Scope and limitations

This is currently an exploratory mechanistic model, **not** an in-silico human embryo, organogenesis model, digital twin, clinical instrument, or tool for predicting real embryonic outcomes. Visual similarity to embryonic morphology is not evidence of biological validity. Every future mechanism should have a defined mathematical model, assumptions, units or declared dimensionless conventions, tests, and (when biological fidelity is claimed) appropriate empirical references.

A useful result may be either an emergent pattern **or a documented failure to produce one** under a specified model. The project's goal is to understand the consequences of the rules, not to guarantee an embryo-shaped output.
