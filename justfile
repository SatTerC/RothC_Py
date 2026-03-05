docs:
  uv run --with matplotlib docs/plots.py
  zensical build

test:
  pytest

time:
  pytest -k timing -s

lint:
  ruff format
  ruff check

