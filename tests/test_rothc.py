import numpy as np
import pandas as pd
import pytest
import time
from pathlib import Path

from rothc_py import RothC
from rothc_py.main import CarbonState


FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def input_file():
    return FIXTURE_DIR / "example_inputs.dat"


@pytest.fixture
def rothc_params(input_file):
    df_head = pd.read_csv(
        input_file,
        skiprows=3,
        header=0,
        nrows=1,
        index_col=None,
        sep=r"\s+",
    )
    return {
        "clay": df_head.loc[0, "clay"],
        "depth": df_head.loc[0, "depth"],
        "iom": float(df_head.loc[0, "iom"]),
    }


@pytest.fixture
def rothc_data(input_file) -> tuple[dict, dict]:
    df = pd.read_csv(input_file, skiprows=6, header=0, index_col=None, sep=r"\s+")
    df.columns = [
        "t_year",
        "t_month",
        "t_mod",
        "t_tmp",
        "t_rain",
        "t_evap",
        "t_C_Inp",
        "t_FYM_Inp",
        "t_PC",
        "t_DPM_RPM",
    ]
    data = {col: df[col].tolist() for col in df.columns}

    # Split: first 12 months for spin-up, rest for forward
    spinup_data = {k: v[:12] for k, v in data.items()}
    forward_data = {k: v[12:] for k, v in data.items()}

    return spinup_data, forward_data


@pytest.fixture
def expected_results():
    return pd.read_csv(FIXTURE_DIR / "month_results.csv")


@pytest.fixture
def expected_spun_up_state():
    return CarbonState(
        dpm=0.14546618698414293,
        rpm=5.678120858752452,
        bio=0.7405937979752077,
        hum=27.642769420831222,
        iom=3.0041,
        soc=37.211050264543026,
        dpm_rc_age=0.348556574230139,
        rpm_rc_age=7.776466544039556,
        bio_rc_age=22.431439734397333,
        hum_rc_age=137.3962355231855,
        iom_age=50000.0,
        total_rc_age=787.4192951123349,
        swc=0.0,
    )


def test_output_matches_expected(rothc_params, rothc_data, expected_results):
    spinup_data, forward_data = rothc_data
    model = RothC(**rothc_params)
    _, actual_results = model(spinup_data, forward_data)

    actual_results = pd.DataFrame(actual_results)

    np.testing.assert_allclose(
        actual_results.values, expected_results.values, rtol=1e-10
    )


def test_spin_up_matches_expected(rothc_params, rothc_data, expected_spun_up_state):
    spinup_data, _ = rothc_data
    model = RothC(**rothc_params)
    state, n_cycles = model.spin_up(spinup_data)

    # NOTE: original counted n_cycles in a very strange way: one iter == one cycle, but cycles
    # were initialised at -1 for no apparent reason. Now, we count one pass through the spin-up
    # data as one cycle, which is 12 months. Hence, n_cycles * 12 - 1 must equal the value for
    # n_cycles from the original code.
    assert n_cycles * 12 - 1 == 19871

    assert state == expected_spun_up_state


def test_timing(rothc_params, rothc_data):
    spinup_data, forward_data = rothc_data
    start = time.perf_counter()
    RothC(**rothc_params)(forward_data, spinup_data)
    end = time.perf_counter()

    elapsed = end - start
    print(f"\nRothC execution time: {elapsed:.4f} seconds")
