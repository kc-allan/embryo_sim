import argparse
import json
import sys

from biology.cell import FATE_NAMES
from core.simulation import Simulation
from visualization.renderer import render_svg


def _format_progress(simulation: Simulation) -> str:
    summary = simulation.summary()
    fates = ", ".join(
        f"{name.lower().replace(' ', '_')}={summary[name]}"
        for name in ("Stem", "Type A", "Type B", "Type C")
    )
    return f"step {summary['tick']}: cells={summary['cells']} ({fates})"


def _print_cell_trace(simulation: Simulation, cell_id: int) -> None:
    cells = {cell.id: cell for cell in simulation.world.cells}
    cell = cells.get(cell_id)
    if cell is None:
        print(f"cell {cell_id} was not created in this run", file=sys.stderr)
        return

    lineage = []
    ancestor = cell
    while ancestor is not None:
        lineage.append(ancestor.id)
        ancestor = cells.get(ancestor.parent_id)
    lineage.reverse()
    lineage_set = set(lineage)
    print(f"cell {cell_id} lineage: {' -> '.join(map(str, lineage))}", file=sys.stderr)

    lineage_events = [
        event for event in simulation.events
        if (
            event["type"] == "fate_transition"
            and event.get("cell_id") in lineage_set
        ) or (
            event["type"] == "division"
            and (
                event.get("daughter_id") in lineage_set
                or event.get("parent_id") == cell_id
            )
        )
    ]
    for event in lineage_events:
        print(json.dumps(event, sort_keys=True), file=sys.stderr)

    print("tick time age fate x y signal gene_a gene_b cycle", file=sys.stderr)
    for state in cell.history:
        x, y = state["position"]
        gene_a, gene_b = state["gene_expression"]
        print(
            f"{state['tick']:4d} {state['time']:5.2f} {state['age']:5.2f} "
            f"{state['fate']:6s} {x:7.3f} {y:7.3f} {state['signal']:7.4f} "
            f"{gene_a:7.4f} {gene_b:7.4f} {state['cycle_progress']:6.3f}",
            file=sys.stderr,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a small 2D cell development simulation.")
    parser.add_argument("--steps", type=int, default=240, help="number of simulation steps")
    parser.add_argument("--seed", type=int, default=7, help="random seed")
    parser.add_argument("--output", default="output/development.svg", help="SVG snapshot path")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase run detail (repeat up to three times)")
    parser.add_argument("--trace-cell", type=int, help="print the complete lineage and state history for a cell ID")
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
                    mass = simulation.world.mass_history[-1]
                    message += (
                        f" signal(min={field.min():.4f}, max={field.max():.4f}, "
                        f"mean={field.mean():.4f}, mass={mass['final_mass']:.4f})"
                        f" mass_delta(diffusion={mass['diffusion_delta']:+.5f}, "
                        f"degradation={mass['degradation_delta']:+.5f}, "
                        f"secretion={mass['secretion_added']:+.5f}, "
                        f"uptake={mass['uptake_removed']:.5f}, "
                        f"clipping={mass['clipping_removed']:.5f}, "
                        f"boundary={mass['boundary_removed']:.5f}, "
                        f"residual={mass['residual']:+.2e})"
                    )
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

    if args.trace_cell is not None:
        _print_cell_trace(simulation, args.trace_cell)


if __name__ == "__main__":
    main()
