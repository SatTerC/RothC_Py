# RothC_Py: Rothamsted Carbon Model in Python

RothC_Py is an implementation of the **Rothamsted Carbon Model (RothC)** in Python. It is a fork of [github.com/Rothamsted-Models/RothC_Py](https://github.com/Rothamsted-Models/RothC_Py) `v1.0.0`, which it itself a translation from a Fortran model that has been around in some form or another since the 1970s.

## Comparison with [Rothamsted-Models/RothC_Py](https://github.com/Rothamsted-Models/RothC_Py)

From the perspective of a Python user, the main advantages offered by this version with respect to the original are:

- Properly packaged and hence installable via `pip` or `uv`
- Faster (~20x)
- Dependency-free
- More concise, clear and documented code

This version does not introduce anything new or fancy, either from a scientific or software perspective.

There are a couple of tests in `tests/test_rothc.py` to check for consistency with the original implementation, which can be run using `pytest` from the root of the repository.

If you care about the details of what's changed it's mostly contained in [Pull Request #1](https://github.com/SatTerC/RothC_Py/pull/1).

The following files have been preserved from the original repo:

- FUNDING.md
- LICENSE
- README_v1.md 
- RothC_description.pdf


## Documentation

The documentation (WIP) is hosted at [satterc.github.io/rothc_py](https://satterc.github.io/RothC_Py/index.html).

## Scientific Background

The RothC model splits soil organic carbon (SOC) into five distinct compartments, four active and one inert. The model accounts for soil type (clay content), temperature, moisture, and plant cover to calculate the decay rates of these pools.

### Carbon Pools:

* **DPM**: Decomposable Plant Material
* **RPM**: Resistant Plant Material
* **BIO**: Microbial Biomass
* **HUM**: Humified Organic Matter
* **IOM**: Inert Organic Matter

The decay follows first-order kinetics, where each active pool $i$ evolves according to:


$$C_i(t) = C_{i,0} e^{-k_i \rho t}$$


where $k$ is the decomposition constant and $\rho$ represents the combined rate-modifying factors.


## Developer Instructions

This project uses **[uv](https://docs.astral.sh/uv/)** for dependency management and packaging.

### Prerequisites

* Python 3.13
* `uv` installed (see [docs](https://docs.astral.sh/uv/getting-started/installation/))

### Setup for Development

1. **Clone the repository:**

```bash
git clone https://github.com/SatTerC/RothC_py.git
cd RothC_Py
```


2. **Create a virtual environment and install dependencies:**

```bash
uv sync
```


3. **Activate the environment:**

```bash
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```


### Pre-commit actions

1. **Run the tests:**

```bash
uv run pytest
```

2. **After making changes, run `ruff`:**

```bash
ruff format
ruff check
```

### Building the docs

Build the docs with 

```bash
zensical build
```

Next, open `site/index.html` in your browser.

See [zensical.org](https://zensical.org/) for more details.


### Useful short-cuts

The awesome [`just`](https://github.com/casey/just) is a development dependency that will be installed when you run `uv sync`.

You can run the following commands anywhere in the repository:

```bash
just test  # run the tests (pytest)
just time  # time RothC on the example data
just lint  # run the formatter and linter (ruff)
just docs  # build the docs (zensical)
```

## Contributing

At this point I'm not sure whether this fork will ever rejoin Rothamsted-Model/RothC_Py, so I'm not sure whether it makes sense to contribute to this repository or that one.
In case there is interest in contributing maybe just reach out to the author (see the noreply address in pyproject.toml).


## References

* **Bolinder MA, et al. (2007).** An approach for estimating net primary productivity and annual carbon inputs to soil for common agricultural crops in Canada. *Agriculture, Ecosystems & Environment*, 118: 29-42.
* **Farina R, et al. (2013).** Modification of the RothC model for simulations of soil organic C dynamics in dryland regions. *Geoderma*, 200: 18-30.
* **Giongo V, et al. (2020).** Optimizing multifunctional agroecosystems in irrigated dryland agriculture to restore soil carbon - Experiments and modelling. *Science of the Total Environment*, 725.
* **Jenkinson DS. (1990).** The Turnover of Organic-Carbon and Nitrogen in Soil. *Philosophical Transactions of the Royal Society of London, Series B: Biological Sciences*, 329: 361-368.
* **Jenkinson DS, et al. (1987).** Modelling the turnover of organic matter in long-term experiments at Rothamsted. *INTECOL Bulletin*, 15: 1-8.
* **Jenkinson DS, Rayner JH. (1977).** Turnover of soil organic matter in some of the Rothamsted classical experiments. *Soil Science*, 123: 298-305.
