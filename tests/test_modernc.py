from datetime import datetime
from rothc_py.modernc import percent_modern_c, DECAY_LAMBDA


class TestModernC:
    """Tests for the percent_modern_c function."""

    def test_pre_1939_returns_all_100(self):
        """Starting before 1939 should return 100% for all months."""
        result = percent_modern_c(datetime(1860, 1, 1), 12)
        assert all(v == 100.0 for v in result)

    def test_pre_1939_then_1939_data(self):
        """Starting in 1938 should have 100s then actual 1939 data."""
        result = percent_modern_c(datetime(1938, 1, 1), 24)
        # First 12 months should be 100
        assert all(v == 100.0 for v in result[:12])
        # Next 12 months should be 97.5 (1939 data)
        assert all(v == 97.5 for v in result[12:])

    def test_within_data_range_1990(self):
        """Starting in 1990 should return correct observed values."""
        result = percent_modern_c(datetime(1990, 1, 1), 12)
        # 1990 values are all 115.0
        assert all(v == 115.0 for v in result)

    def test_within_data_range_2000(self):
        """Starting in 2000 should return correct observed values."""
        result = percent_modern_c(datetime(2000, 1, 1), 12)
        # 2000 values are all 109.0
        assert all(v == 109.0 for v in result)

    def test_bomb_peak_1963_1964(self):
        """The bomb peak should show high values around 1963-1964."""
        result = percent_modern_c(datetime(1963, 1, 1), 24)
        # 1963 should be ~194, 1964 should be ~196
        assert all(v == 194.0 for v in result[:12])
        assert all(v == 196.0 for v in result[12:])

    def test_extrapolation_starts_from_2007_12(self):
        """Extrapolation should start from the 2007-12 baseline (106.8)."""
        result = percent_modern_c(datetime(2007, 12, 1), 3)
        # First value should be the observation (106.8), then decay
        assert result[0] == 106.8
        # Second value should be slightly decayed
        assert 106.7 < result[1] < 106.8
        assert result[1] < result[0]

    def test_extrapolation_decay(self):
        """Extrapolated values should decay exponentially toward 100."""
        result = percent_modern_c(datetime(2008, 1, 1), 36)
        # Values should progressively decrease
        assert result[0] > result[11] > result[23] > result[35]
        # All should be above 100
        assert all(v > 100.0 for v in result)
        # At 36 months, should be around 105.6
        assert 105.0 < result[35] < 106.0

    def test_extrapolation_approaches_100(self):
        """Over long time periods, values should approach 100."""
        result = percent_modern_c(datetime(2020, 1, 1), 120)  # 10 years
        # Should be declining toward 100
        assert result[0] > result[-1] > 100.0
        # After 10 years (120 months), should be around 101-102
        assert 101.0 < result[-1] < 103.0


class TestDecayLambda:
    """Tests for the DECAY_LAMBDA constant."""

    def test_decay_lambda_is_positive(self):
        """DECAY_LAMBDA should be positive for decay behavior."""
        assert DECAY_LAMBDA > 0

    def test_decay_lambda_reasonable_magnitude(self):
        """DECAY_LAMBDA should be on order 0.005 per month."""
        assert 0.001 < DECAY_LAMBDA < 0.01

    def test_exponential_decay_formula(self):
        """Verify the exponential decay formula gives expected results."""
        from math import exp

        excess = 6.8  # 106.8 - 100
        # After 1 month
        expected_1 = 100 + excess * exp(-DECAY_LAMBDA * 1)
        # After 12 months (1 year)
        expected_12 = 100 + excess * exp(-DECAY_LAMBDA * 12)

        result = percent_modern_c(datetime(2007, 12, 1), 13)
        assert abs(result[1] - expected_1) < 0.001
        assert abs(result[12] - expected_12) < 0.001


class TestEdgeCases:
    """Edge case tests."""

    def test_single_month(self):
        """Should work with n_months=1."""
        result = percent_modern_c(datetime(2000, 1, 1), 1)
        assert len(result) == 1
        assert result[0] == 109.0

    def test_exactly_at_first_obs(self):
        """Starting exactly at first observation (1939-01) should work."""
        result = percent_modern_c(datetime(1939, 1, 1), 12)
        assert result[0] == 97.5

    def test_exactly_at_last_obs(self):
        """Starting exactly at last observation (2007-12) should work."""
        result = percent_modern_c(datetime(2007, 12, 1), 12)
        assert result[0] == 106.8
