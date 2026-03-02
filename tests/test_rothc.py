import numpy as np
import pandas as pd
import pytest
import time
from pathlib import Path

from rothc_py import RothC


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
def expected_years():
    return pd.read_csv(FIXTURE_DIR / "year_results.csv")


@pytest.fixture
def expected_months():
    return pd.read_csv(FIXTURE_DIR / "month_results.csv")


def test_output_matches_expected(
    rothc_params, rothc_data, expected_years, expected_months
):
    spinup_data, forward_data = rothc_data
    year_results, month_results = RothC(**rothc_params)(forward_data, spinup_data)

    actual_years = pd.DataFrame(year_results)
    actual_months = pd.DataFrame(month_results)

    np.testing.assert_allclose(actual_years.values, expected_years.values, rtol=1e-10)
    np.testing.assert_allclose(actual_months.values, expected_months.values, rtol=1e-10)


def test_final_year_values(rothc_params, rothc_data):
    spinup_data, forward_data = rothc_data
    year_results, _ = RothC(**rothc_params)(forward_data, spinup_data)
    actual_years = pd.DataFrame(year_results)
    final_row = actual_years.iloc[-1]

    assert final_row["Year"] == 2007
    assert final_row["Month"] == 12
    assert final_row["DPM_t_C_ha"] == pytest.approx(0.1858764195450398, abs=1e-10)
    assert final_row["RPM_t_C_ha"] == pytest.approx(6.370323733617322, abs=1e-10)
    assert final_row["BIO_t_C_ha"] == pytest.approx(0.8285615862243552, abs=1e-10)
    assert final_row["HUM_t_C_ha"] == pytest.approx(27.802662060467732, abs=1e-10)
    assert final_row["IOM_t_C_ha"] == pytest.approx(3.0041, abs=1e-10)
    assert final_row["SOC_t_C_ha"] == pytest.approx(38.19152379985445, abs=1e-10)
    assert final_row["deltaC"] == pytest.approx(-1.6364720117949538, abs=1e-10)


def test_spin_up_output(rothc_params, rothc_data):
    spinup_data, _ = rothc_data
    model = RothC(**rothc_params)
    state, n_cycles = model.spin_up(spinup_data)

    # NOTE: original counted n_cycles in a very strange way: one iter == one cycle, but cycles
    # were initialised at -1 for no apparent reason. Now, we count one pass through the spin-up
    # data as one cycle, which is 12 months. Hence, n_cycles * 12 - 1 must equal the value for
    # n_cycles from the original code.
    assert n_cycles * 12 - 1 == 19871

    assert state.dpm == pytest.approx(0.14546618698414293, abs=1e-10)
    assert state.rpm == pytest.approx(5.678120858752452, abs=1e-10)
    assert state.bio == pytest.approx(0.7405937979752077, abs=1e-10)
    assert state.hum == pytest.approx(27.642769420831222, abs=1e-10)
    assert state.iom == pytest.approx(3.0041, abs=1e-10)
    assert state.soc == pytest.approx(37.211050264543026, abs=1e-10)
    assert state.dpm_rc_age == pytest.approx(0.348556574230139, abs=1e-10)
    assert state.rpm_rc_age == pytest.approx(7.776466544039556, abs=1e-10)
    assert state.bio_rc_age == pytest.approx(22.431439734397333, abs=1e-10)
    assert state.hum_rc_age == pytest.approx(137.3962355231855, abs=1e-10)
    assert state.iom_age == pytest.approx(50000.0, abs=1e-10)
    assert state.total_rc_age == pytest.approx(787.4192951123349, abs=1e-10)
    assert state.swc == pytest.approx(0.0, abs=1e-10)


def test_timing(rothc_params, rothc_data):
    spinup_data, forward_data = rothc_data
    start = time.perf_counter()
    RothC(**rothc_params)(forward_data, spinup_data)
    end = time.perf_counter()

    elapsed = end - start
    print(f"\nRothC execution time: {elapsed:.4f} seconds")
