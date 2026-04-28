# Default recipe to display available commands.
_:
  @just --list

# Format and lint the package using ruff.
lint:
  ruff format
  ruff check

# Build the documentation using Zensical
docs:
  uv run --with matplotlib docs/plots.py
  zensical build

# Run the test suite.
test:
  pytest

# Run a timing benchmark.
time:
  pytest -k timing -s
