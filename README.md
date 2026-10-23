# Embryo Sim

A small, reproducible 2D agent-based model for exploring how local cell rules can produce spatial patterns. It starts from one virtual cell and combines cell division, a diffusing signal, simple gene expression, fate commitment, and soft cell mechanics.

This is an educational toy model. Its rules and parameters are illustrative and are not validated descriptions of human or other real embryonic development.

## Run

Python 3.12 and NumPy are required. From this directory:

```bash
.venv/bin/python main.py --steps 240 --seed 7 --output output/development.svg
```

The command prints population and fate counts as JSON and writes a standalone SVG snapshot. Open the SVG in a browser to view it. The blue haze shows morphogen concentration, cell colors show fate, and the pale ring marks the signal source.

Use `--steps`, `--seed`, and `--output` to change the run length, random seed, or snapshot path. Add `-v` for progress updates, `-vv` for per-step population and fate counts, or `-vvv` for per-step signal ranges and a final per-cell detail table. Verbose diagnostics go to stderr, leaving the JSON result on stdout. For example:

```bash
.venv/bin/python main.py --steps 240 -vv
```

The engine can also be used directly:

```python
from core.simulation import Simulation

simulation = Simulation(seed=7)
simulation.run(240)
print(simulation.summary())
```

The main rules live in `biology/`; the spatial environment, scheduler, mechanics, and simulation loop live in `core/`.
