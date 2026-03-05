import csv
from datetime import datetime
from math import exp
from pathlib import Path
from typing import Generator

csv_path = Path(__file__).parent / "data" / "modern.csv"

DECAY_LAMBDA = 0.00513
"""Characteristic timescale for the exponential decay of carbon post 1964 in months^(-1)."""


def _read_modern_csv() -> tuple[list[datetime], list[float]]:
    dates = []
    moderns = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dates.append(datetime.strptime(row["date"], "%Y-%m-%d"))
            moderns.append(float(row["modern"]))
    return dates, moderns


def _modern_values(start_year: int, start_month: int) -> Generator[float, None, None]:
    obs_dates, obs_values = _read_modern_csv()

    start_date = datetime(start_year, start_month, 1)
    first_obs = obs_dates[0]
    last_obs = obs_dates[-1]

    if start_date < first_obs:
        months_before_first_obs = (12 * first_obs.year + first_obs.month) - (
            12 * start_date.year + start_date.month
        )
        historic_values = [100.0 for _ in range(months_before_first_obs)] + obs_values
    elif start_date <= last_obs:
        months_before_last_obs = (12 * last_obs.year + last_obs.month) - (
            12 * start_date.year + start_date.month
        )
        historic_values = reversed(
            [
                obs
                for _, obs in zip(
                    range(months_before_last_obs + 1), reversed(obs_values)
                )
            ]
        )
    else:
        historic_values = []

    for value in historic_values:
        yield value

    excess_at_last_obs = obs_values[-1] - 100

    if start_date <= last_obs:
        months_since_last_obs = 1  # start 1 month after last obs
    else:
        months_since_last_obs = (12 * start_date.year + start_date.month) - (
            12 * last_obs.year + last_obs.month
        )

    while True:
        yield 100 + excess_at_last_obs * exp(-DECAY_LAMBDA * months_since_last_obs)
        months_since_last_obs += 1


def modern_c(start_date: datetime, n_months: int) -> list[float]:
    """
    Returns the % modern radiocarbon values for a given time period.

    Uses the actual atmospheric radiocarbon data extracted from example_inputs.dat,
    which reproduces the curve of % modern C as in Fig 5 of the RothC description.

    Parameters
    ----------
    start_date : datetime
        The start date for the time series
    n_months : int
        Number of months to return

    Returns
    -------
    list[float]
        List of % modern radiocarbon values for each month

    Notes
    -----
    Only the `year` and `month` attributes of `start_date` are used. Any `day`
    or `time` information will be ignored.
    """
    assert n_months >= 1
    return [
        value
        for _, value in zip(
            range(n_months), _modern_values(start_date.year, start_date.month)
        )
    ]
