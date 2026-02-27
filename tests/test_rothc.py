import numpy as np
import pandas as pd
import pytest
from pathlib import Path
import tempfile

from rothc_py.main import main


FIXTURE_DIR = Path(__file__).parent / "fixtures"


def test_output_matches_expected():
    """Test that main() produces expected output matching fixtures."""
    input_path = FIXTURE_DIR / "example_inputs.dat"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        main(input_path, output_dir)

        expected_years = pd.read_csv(FIXTURE_DIR / "year_results.csv")
        actual_years = pd.read_csv(output_dir / "year_results.csv")

        expected_months = pd.read_csv(FIXTURE_DIR / "month_results.csv")
        actual_months = pd.read_csv(output_dir / "month_results.csv")

        np.testing.assert_allclose(
            actual_years.values, expected_years.values, rtol=1e-10
        )
        np.testing.assert_allclose(
            actual_months.values, expected_months.values, rtol=1e-10
        )


def test_final_year_values():
    """Test specific final year values match expected results."""
    input_path = FIXTURE_DIR / "example_inputs.dat"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        main(input_path, output_dir)

        actual_years = pd.read_csv(output_dir / "year_results.csv")
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
