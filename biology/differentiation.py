from biology.cell import Cell


def update_fate(cell: Cell, morphogen: float) -> None:
    """Assign a model fate, allowing Type B to refine into Type A."""
    if cell.fate not in (0, 2):
        return
    if morphogen >= 0.64 and cell.gene_expression[0] >= 0.45:
        cell.fate = 1
    elif cell.age >= 5.0 and morphogen >= 0.29:
        cell.fate = 2
    elif cell.age >= 5.0:
        cell.fate = 3
