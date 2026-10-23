from dataclasses import dataclass

import numpy as np

from biology.cell import Cell, FATE_NAMES
from biology.differentiation import update_fate
from biology.gene_regulation import update_expression
from core.mechanics import update_positions
from core.scheduler import advance_cell_cycles
from core.world import World


@dataclass
class SimulationConfig:
    world_size: float = 80.0
    grid_size: int = 128
    dt: float = 0.08
    diffusion: float = 0.8
    decay: float = 0.12
    secretion: float = 4.0
    division_period: float = 1.7
    cell_radius: float = 1.7
    mobility: float = 0.85
    adhesion: float = 0.035
    max_cells: int = 300
    signaling_enabled: bool = True
    proliferation_enabled: bool = True
    chemotaxis_enabled: bool = True


class Simulation:
    """Reproducible 2D cell model for exploring local rules and patterning."""

    def __init__(self, seed: int = 7, config: SimulationConfig | None = None):
        self.seed = seed
        self.config = config or SimulationConfig()
        self.rng = np.random.default_rng(seed)
        self.world = World(size=self.config.world_size, grid_size=self.config.grid_size)
        self.tick = 0
        self.next_cell_id = 1
        self.events: list[dict] = []
        self.world.cells.append(Cell(
            id=0,
            parent_id=None,
            position=np.zeros(2, dtype=float),
            radius=self.config.cell_radius,
            age=0.0,
            cycle_progress=0.0,
            gene_expression=np.zeros(2, dtype=float),
            organizer=True,
        ))

    def step(self) -> World:
        config = self.config
        dt = config.dt
        self.world.update_morphogen(
            dt,
            config.diffusion,
            config.decay,
            config.secretion,
            tick=self.tick + 1,
            enabled=config.signaling_enabled,
        )
        for cell in self.world.cells:
            if not cell.alive:
                continue
            signal = self.world.sample_morphogen(cell.position)
            previous_fate = cell.fate
            previous_expression = cell.gene_expression.copy()
            cell.gene_expression = update_expression(cell.gene_expression, signal, dt)
            update_fate(cell, signal)
            if cell.fate != previous_fate:
                self.events.append({
                    "type": "fate_transition",
                    "tick": self.tick + 1,
                    "time": (self.tick + 1) * dt,
                    "cell_id": cell.id,
                    "from": FATE_NAMES[previous_fate],
                    "to": FATE_NAMES[cell.fate],
                    "age": cell.age,
                    "position": cell.position.tolist(),
                    "signal": signal,
                    "gene_expression_before": previous_expression.tolist(),
                    "gene_expression_after": cell.gene_expression.tolist(),
                    "parent_id": cell.parent_id,
                })

        ready = advance_cell_cycles(
            self.world,
            dt,
            config.division_period,
            config.max_cells,
            enabled=config.proliferation_enabled,
        )
        for parent in ready:
            if len(self.world.cells) >= config.max_cells:
                break
            parent.cycle_progress = max(parent.cycle_progress - 1.0, 0.0)
            angle = self.rng.uniform(0.0, 2.0 * np.pi)
            direction = np.array([np.cos(angle), np.sin(angle)])
            offset = direction * parent.radius * 1.12
            parent.position = parent.position - offset * 0.5
            self.world.cells.append(Cell(
                id=self.next_cell_id,
                parent_id=parent.id,
                position=parent.position + offset,
                radius=parent.radius,
                age=0.0,
                cycle_progress=0.0,
                gene_expression=np.clip(parent.gene_expression + self.rng.normal(0.0, 0.025, 2), 0.0, 1.0),
                fate=parent.fate,
                organizer=False,
            ))
            daughter = self.world.cells[-1]
            self.events.append({
                "type": "division",
                "tick": self.tick + 1,
                "time": (self.tick + 1) * dt,
                "parent_id": parent.id,
                "daughter_id": daughter.id,
                "inherited_fate": FATE_NAMES[daughter.fate],
                "parent_age": parent.age,
                "parent_position": parent.position.tolist(),
                "daughter_position": daughter.position.tolist(),
                "parent_gene_expression": parent.gene_expression.tolist(),
                "daughter_gene_expression": daughter.gene_expression.tolist(),
            })
            self.next_cell_id += 1

        update_positions(
            self.world,
            dt,
            config.mobility,
            config.adhesion,
            chemotaxis=config.chemotaxis_enabled,
        )
        self.tick += 1
        for cell in self.world.cells:
            if not cell.alive:
                continue
            cell.history.append({
                "tick": self.tick,
                "time": self.tick * dt,
                "age": cell.age,
                "position": cell.position.tolist(),
                "fate": FATE_NAMES[cell.fate],
                "signal": self.world.sample_morphogen(cell.position),
                "gene_expression": cell.gene_expression.tolist(),
                "cycle_progress": cell.cycle_progress,
            })
        return self.world

    def run(self, steps: int) -> World:
        for _ in range(max(0, steps)):
            self.step()
        return self.world

    def summary(self) -> dict[str, int | float]:
        living = [cell for cell in self.world.cells if cell.alive]
        counts = {name: 0 for name in ("Stem", "Type A", "Type B", "Type C")}
        for cell in living:
            counts[FATE_NAMES[cell.fate]] += 1
        return {
            "tick": self.tick,
            "time": round(self.tick * self.config.dt, 2),
            "cells": len(living),
            "mean_morphogen": round(float(np.mean(self.world.morphogen)), 4),
            "morphogen_mass": round(float(self.world.morphogen.sum() * self.world.spacing**2), 6),
            **counts,
        }
