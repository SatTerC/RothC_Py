docs:
  zensical build

test:
  pytest

time:
  pytest -k timing -s

lint:
  ruff format
  ruff check

