---
title: RothC_Py
icon: lucide/worm
---

# RothC_Py: Rothamsted Carbon Model in Python

`RothC_Py` is an implementation of the **Rothamsted Carbon Model (RothC)** in Python.
It is a fork of [github.com/Rothamsted-Models/RothC_Py](https://github.com/Rothamsted-Models/RothC_Py) `v1.0.0`, which it itself a translation from a Fortran model that has been around in some form or another since the 1970s.

From the perspective of a Python user, the main advantages offered by this version with respect to the original are:

- Properly packaged and hence installable via `pip` or `uv`
- Faster (~20x)
- Dependency-free
- More concise, clear and documented code

This version does not introduce anything new or fancy, either from a scientific or software perspective.


## Installation

Install the package via `pip` or `uv`.
Currently it is only available from GitHub.

=== "pip"

    ``` sh
    pip install git+https://github.com/satterc/rothc_py
    ```

=== "uv"

    ``` sh
    uv add git+https://github.com/satterc/rothc_py
    ```

This will install a package called `rothc_py` into your environment.

## Basic usage

`rothc_py` provides a class, `RothC`, which can be used to initialise, spin-up and run the model, given user-provided parameters and data.

```python
from rothc_py import RothC

# To do
params = ...
spinup_data = ...
forward_data = ...

model = RothC(**params)

spun_up_state, n_cycles = model.spin_up(spinup_data)
final_state, results = model.forward(spun_up_state, forward_data)
```

Alternatively, one can call the `RothC` instance, which internally just calls `spin_up` followed by `forward`.

```python
model = RothC(**params)

final_state, results = model(spinup_data, forward_data)
```

State objects (`spun_up_state`, `final_state`) are instances of `rothc_py.containers.CarbonState`, which is just a `dataclasses.dataclass` containing the instantaneous state of the model.

The second returned value, `results`, is a dict of lists containing the monthly values.
This can easily be used to create a `pandas.DataFrame` and hence written to e.g. a CSV file.

```python
import pandas as pd

results_df = pd.DataFrame(results)

results_df.to_csv("results.csv")
```

For more details on the contents of this package, see the API reference.
