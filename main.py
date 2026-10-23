import argparse
import json
import sys

from core.simulation import Simulation
from biology.cell import FATE_NAMES
from visualization.renderer import render_svg


def _format_progress(simulation: Simulation) -> str:
    summary = simulation.summary()
    fates = ", ".join(
        f"{name.lower().replace(' ', '_')}={summary[name]}"
        for name in ("Stem", "Type A", "Type B", "Type C")
    )
    return f"step {summary['tick']}: cells={summary['cells']} ({fates})"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a small 2D cell development simulation.")
    parser.add_argument("--steps", type=int, default=240, help="number of simulation steps")
    parser.add_argument("--seed", type=int, default=7, help="random seed")
    parser.add_argument("--output", default="output/development.svg", help="SVG snapshot path")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase run detail (repeat up to three times)")
    args = parser.parse_args()
    verbosity = min(args.verbose, 3)

    simulation = Simulation(seed=args.seed)
    steps = max(args.steps, 0)
    if verbosity == 0:
        simulation.run(steps)
    else:
        if verbosity >= 3:
            config = vars(simulation.config)
            print(f"configuration: {json.dumps(config, sort_keys=True)}", file=sys.stderr)

        progress_interval = max(1, (steps + 9) // 10)
        for step in range(1, steps + 1):
            simulation.step()
            if verbosity >= 2:
                message = _format_progress(simulation)
                if verbosity >= 3:
                    field = simulation.world.morphogen
                    message += f" signal(min={field.min():.4f}, max={field.max():.4f}, mean={field.mean():.4f})"
                print(message, file=sys.stderr)
            elif step % progress_interval == 0 or step == steps:
                percent = round(step * 100 / steps)
                print(f"{percent:3d}% {_format_progress(simulation)}", file=sys.stderr)

    output = render_svg(simulation.world, simulation.tick, args.output)
    print(json.dumps({**simulation.summary(), "snapshot": str(output)}, indent=2))

    if verbosity >= 3:
        print("final cells (id parent fate x y age signal gene_a gene_b):", file=sys.stderr)
        for cell in simulation.world.cells:
            signal = simulation.world.sample_morphogen(cell.position)
            print(
                f"{cell.id:4d} {str(cell.parent_id):>6s} {FATE_NAMES[cell.fate]:>6s} "
                f"{cell.position[0]:7.2f} {cell.position[1]:7.2f} {cell.age:6.2f} "
                f"{signal:6.3f} {cell.gene_expression[0]:6.3f} {cell.gene_expression[1]:6.3f}",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
