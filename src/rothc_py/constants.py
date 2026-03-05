"""Model constants for the RothC soil carbon model.

This module contains all physical and mathematical constants used by the RothC
model, including radiocarbon parameters, decomposition rates, carbon flow
fractions, and rate modifying factors.
"""

# =============================================================================
# Radiocarbon Constants
# =============================================================================

RADIO_HALFLIFE = 5568.0
"""Conventional half-life of radiocarbon-14 (years)."""

RADIO_MEAN_LIFETIME = 8035.0
"""Mean lifetime of radiocarbon-14, used in δ14C calculations (years)."""

IOM_INITIAL_AGE = 50000.0
"""Default initial age of inert organic matter, implying negligible 14C (years)."""

MODERN_C_DECAY_LAMBDA = 0.00513
r"""Empirical timescale for the removal of atmospheric ^14^C in months^-1^.

    This is based on fitting an single exponential decary,

    $ x(t) = x(0) e^{-\lambda t} $

    to two (time, excess modern C) points: (1964, 96%) and (2007, 6.8%), which
    yields λ = 0.616 years^-1^ or λ = 0.00513 months^-1^.
"""


# =============================================================================
# Decomposition Rate Constants
# =============================================================================

ZERO_THRESHOLD = 1e-8
"""Minimum value threshold for numerical stability in radiocarbon age calculations."""

DPM_RATE = 10.0
"""Decomposition rate constant for Decomposable Plant Material (yr⁻¹)."""

RPM_RATE = 0.3
"""Decomposition rate constant for Resistant Plant Material (yr⁻¹)."""

BIO_RATE = 0.66
"""Decomposition rate constant for Microbial Biomass (yr⁻¹)."""

HUM_RATE = 0.02
"""Decomposition rate constant for Humified Organic Matter (yr⁻¹)."""


# =============================================================================
# Carbon Flow Fractions
# =============================================================================

CLAW_A = 1.67
"""Scaling factor in the clay-dependent CO2/(BIO+HUM) partitioning equation."""

CLAW_B = 1.85
"""Coefficient in the clay-dependent CO2/(BIO+HUM) partitioning equation."""

CLAW_C = 1.60
"""Coefficient in the clay-dependent CO2/(BIO+HUM) partitioning equation."""

CLAW_D = 0.0786
"""Exponential decay coefficient in the clay-dependent partitioning equation."""

FRAC_TO_BIO = 0.46
"""Fraction of (BIO+HUM) that becomes Microbial Biomass."""

FRAC_TO_HUM = 0.54
"""Fraction of (BIO+HUM) that becomes Humified Organic Matter."""

FYM_FRAC_DPM = 0.49
"""Fraction of farmyard manure carbon that enters DPM pool."""

FYM_FRAC_RPM = 0.49
"""Fraction of farmyard manure carbon that enters RPM pool."""

FYM_FRAC_HUM = 0.02
"""Fraction of farmyard manure carbon that enters HUM pool."""


# =============================================================================
# Moisture Rate Modifier Constants
# =============================================================================

RMF_MOIST_MAX = 1.0
"""Maximum rate modifying factor for moisture (no water stress)."""

RMF_MOIST_MIN = 0.2
"""Minimum rate modifying factor for moisture (maximum water stress)."""

SMD_COEFF_A = 20.0
"""Coefficient A in the maximum soil moisture deficit equation."""

SMD_COEFF_B = 1.3
"""Coefficient B in the maximum soil moisture deficit equation."""

SMD_COEFF_C = 0.01
"""Coefficient C in the maximum soil moisture deficit equation."""

SMD_DEPTH_DIVISOR = 23.0
"""Reference soil depth for scaling maximum soil moisture deficit (cm)."""

SMD_1BAR_FRAC = 0.444
"""Fraction of maximum TSMD at which rate modifier equals 1.0."""

SMD_BARE_FRAC = 0.556
"""Inverse of bare soil divisor (1/1.8) for reduced evaporation from bare soil."""

EVAP_FACTOR = 0.75
"""Factor to convert open-pan evaporation to evapotranspiration."""


# =============================================================================
# Temperature Rate Modifier Constants
# =============================================================================

TEMP_MIN = -5.0
"""Minimum temperature (°C) below which decomposition rate is zero."""

JENKINSON_A = 47.91
"""Coefficient A in the Jenkinson temperature rate modifier equation."""

JENKINSON_B = 106.06
"""Coefficient B in the Jenkinson temperature rate modifier equation."""

JENKINSON_C = 18.27
"""Coefficient C in the Jenkinson temperature rate modifier equation."""


# =============================================================================
# Plant Cover Rate Modifier Constants
# =============================================================================

RMF_PC_BARE = 1.0
"""Rate modifying factor for bare soil (no plant cover)."""

RMF_PC_COVERED = 0.6
"""Rate modifying factor for covered soil (with plant cover)."""


# =============================================================================
# Simulation Constants
# =============================================================================

MONTHS_PER_YEAR = 12
"""Number of months per year."""

EQUILIBRIUM_THRESHOLD = 1e-6
"""Threshold for spin-up convergence: maximum annual TOC change (t C/ha)."""
