from pathlib import Path
from xml.sax.saxutils import escape

from biology.cell import FATE_COLORS, FATE_NAMES
from core.world import World


def render_svg(world: World, tick: int, output: str | Path) -> Path:
    """Write a standalone SVG with the morphogen field and cell states."""
    path = Path(output)
    width, height = 900, 760
    margin, plot_size = 60, 640
    scale = plot_size / world.size
    cell_size = world.size / world.grid_size
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#101821"/>',
        '<text x="60" y="34" fill="#f0f2ed" font-family="sans-serif" font-size="20" font-weight="600">Cell development</text>',
        f'<text x="840" y="34" text-anchor="end" fill="#aeb9c5" font-family="sans-serif" font-size="13">tick {tick} · {len(world.cells)} cells</text>',
        f'<rect x="{margin}" y="{margin}" width="{plot_size}" height="{plot_size}" fill="#141f2a" stroke="#344353"/>',
    ]
    field = world.morphogen
    stride = max(1, world.grid_size // 56)
    rect_size = cell_size * scale * stride + 0.5
    for iy in range(0, world.grid_size, stride):
        for ix in range(0, world.grid_size, stride):
            intensity = float(field[iy, ix])
            if intensity < 0.025:
                continue
            x = margin + (ix + 0.5) * cell_size * scale
            y = margin + (world.grid_size - iy - 0.5) * cell_size * scale
            opacity = min(0.44, intensity * 0.45)
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{rect_size:.1f}" height="{rect_size:.1f}" fill="#53b7d2" opacity="{opacity:.3f}"/>')

    for cell in world.cells:
        if not cell.alive:
            continue
        x = margin + (cell.position[0] + world.size / 2) * scale
        y = margin + (world.size / 2 - cell.position[1]) * scale
        radius = cell.radius * scale
        color = FATE_COLORS[cell.fate]
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{color}" stroke="#101821" stroke-width="1.5"/>')
        if cell.organizer:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius + 3:.2f}" fill="none" stroke="#f3e5a2" stroke-width="1.5"/>')

    legend_x = 740
    parts.append('<text x="740" y="94" fill="#f0f2ed" font-family="sans-serif" font-size="14" font-weight="600">Cell fates</text>')
    for index, name in enumerate(FATE_NAMES):
        y = 122 + index * 28
        parts.append(f'<circle cx="{legend_x + 7}" cy="{y - 4}" r="6" fill="{FATE_COLORS[index]}"/>')
        parts.append(f'<text x="{legend_x + 22}" y="{y}" fill="#c5cfd8" font-family="sans-serif" font-size="13">{escape(name)}</text>')
    parts.append('<circle cx="747" cy="252" r="9" fill="none" stroke="#f3e5a2" stroke-width="1.5"/>')
    parts.append('<text x="762" y="256" fill="#c5cfd8" font-family="sans-serif" font-size="13">Signal source</text>')
    parts.append('<rect x="740" y="286" width="14" height="14" fill="#53b7d2" opacity="0.4"/>')
    parts.append('<text x="762" y="298" fill="#c5cfd8" font-family="sans-serif" font-size="13">Morphogen field</text>')
    parts.append('</svg>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(parts), encoding='utf-8')
    return path
