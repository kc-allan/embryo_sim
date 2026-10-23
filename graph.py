from core.simulation import Simulation
import matplotlib.pyplot as plt

simulation = Simulation(seed=42)
# simulation.config.max_cells = 3000
steps = 500
sample_interval = 5

x = []

for step in range(steps):
	simulation.step()
	if step % sample_interval == 0:
		print(f"Step {step}: {len(simulation.world.cells)} cells")
		summary = simulation.summary()
		x.append(summary["cells"])

y = [(step + 1) * simulation.config.dt for step in range(0, steps, sample_interval)]

plt.plot(x, y)
plt.show()