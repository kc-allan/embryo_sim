# Verification and validation record

## Scope

These checks verify software and numerical properties of the declared models.
They do not validate the simulator against real embryological observations.

## Automated checks

Run with:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

The suite checks capacity and disabled-proliferation semantics; no-flux mass
conservation, first-order decay, source mass, signaling-ablation removal, and
positivity-bound rejection;
grid and time refinement of no-flux cosine modes; finite-range and coincident-cell
mechanics; reference/index force equivalence; spatial neighbor counts;
recording-mode trajectory equivalence; toggle fixed points and pulse
withdrawal and repeated signal changes. The grid refinement test decreases
analytical-mode error across 20, 40, and 80 cells per side; the timestep test
decreases Euler error against the semidiscrete solution.

## Reference behavior before this milestone

Before implementation changes, fixed-population timings were captured once on
the local Python 3.12 / NumPy environment using sparse random positions in a
1000-unit square. Approximate seconds per update:

| Cells | All-pairs mechanics | All-pairs scheduler |
|---:|---:|---:|
| 100 | 0.0031 | 0.0220 |
| 300 | 0.0143 | 0.1981 |
| 1,000 | 0.3625 | 2.2177 |
| 3,000 | 3.9930 | 20.3766 |

These are single-run pre-change timings and are machine-specific.

## Revised implementation performance

`analysis.benchmarks` measures fixed populations with three timed repetitions
after warm-up; reference pair forces are also checked against indexed forces
at 100 and 300 cells. On the same environment, indexed force calculation and
scheduler results were:

| Cells | Indexed mechanics forces | Indexed scheduler |
|---:|---:|---:|
| 100 | 0.00074 | 0.00074 |
| 300 | 0.00218 | 0.00237 |
| 1,000 | 0.00767 | 0.00825 |
| 3,000 | 0.03770 | 0.03229 |

Against the older single-run timings, the 3,000-cell position update is about
103x faster and scheduler about 631x faster. The scheduler semantics also
changed to age cells at capacity. These results support improved scalability,
not a universal complexity guarantee.

The benchmark also measures neighbor queries, field updates, one gene update,
synthetic event creation, JSON serialization, and SVG rendering. At 3,000
cells these measured approximately 0.0173, 0.00087, 0.000012 per gene update,
0.00053, 0.01045, and 0.00906 seconds respectively. Serialization and
rendering depend on payload size and output environment.

## Developmental experiments

`analysis.experiments` reports fate counts and proportions; signal mean and
quartiles sampled at each cell's final position, grouped by current fate;
signal at the update where fate transitions were recorded, grouped by
destination fate; cell-type centroids; the mean distance of cells from the
colony centroid; and the mean pairwise distance between observed fate
centroids. The centroid distance is a descriptive spatial summary, not a
correlation, sorting coefficient, or causal estimate. No numeric ordering is
assigned to fate labels.

For the revised mechanics and field baseline, seed 7, 240 steps, `dt=0.08`
produced 81 cells: 10 Stem, 11 Type A, no Type B, and 60 Type C. There were 80
divisions and 24 fate transitions. Across 19.2 model-time units integrated
secretion was 30.47430095, degradation was -18.55378738, diffusion delta was
numerically zero, and maximum ledger residual was `3.55e-15`. The difference
from the historical 59-cell run is expected because unbounded-attraction
mechanics was replaced; this is a scientific model revision, not an
equivalent optimization.

Mechanism ablations and robustness runs are generated with
`.venv/bin/python -m analysis.experiments`. In the revised seed-7 baseline,
signaling-off removed Type A assignments. Proliferation-off kept the founder
only. Chemotaxis-off had little endpoint effect. Adhesion-off reduced mean
fate-centroid separation from 6.64 to 8.26 model-length units in this run and
changed the population; it does not demonstrate
cell-type-dependent sorting, because that mechanism is absent. Seeds 1 and 17,
and a half timestep, are included in the generated report. Stochastic
trajectories are not expected to match across different timesteps.

## Regulatory experiment

The separate toggle model defaults yield numerical fixed points approximately
`(0.118, 2.000)`, `(1,1)`, and `(2.000,0.118)`. Jacobian eigenvalues classify
the outer states as stable and the middle state as a saddle. A unit signal
pulse of 10 model-time units from `(0.05,2.0)` produces an A-high state; after
signal withdrawal and 30 more units it remains near `(2.000,0.118)`. This
demonstrates hysteresis-like persistence for these equations and parameters
only.

## Not yet verified

The full model has not been calibrated or validated biologically. Automated
stochastic ensembles, arrangement sensitivity, event-rate convergence, broad
2D spatial convergence, and a complete timestep convergence study remain
future work. Source deposition and field integration are verified for mass
and an analytical diffusion mode, not every coupled dynamic.
