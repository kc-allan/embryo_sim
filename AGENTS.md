# Agent Operating Rules

## Purpose

Keep work reproducible, scientifically explicit, and aligned with the project roadmap.

- Read `README.md`, `VISION.md`, `MISSION.md`, and `PLAN.md` before development.
- Separate biological hypotheses, equations, numerical methods, and infrastructure changes.
- Add focused `unittest` coverage for behavior and numerical claims; run it with `.venv/bin/python -m unittest discover -s tests`.
- Run relevant experiment scripts and record configuration, seed, and numerical settings.
- Keep the model specification, numerical methods, validation, architecture, mission, and plan synchronized with actual code.
- Do not claim biological validation from code tests, conservation checks, or visual output.
- Do not commit unless the user explicitly asks.
