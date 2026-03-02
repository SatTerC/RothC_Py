"""
Rothamsted Carbon Model (RothC) implemented in Python.

What follows is the original description from https://github.com/Rothamsted-Models/RothC_Py/RothC_Py.py,
at commit bd90ce3cf616d5316042b73a3b1f09c5b6e3b361 (Mar 12, 2025).

######################################################################################################################
#
#  RothC python version
#
#  This python version was translated from the Fortran code by Alice Milne, Jonah Prout and Kevin Coleman 29/02/2024
#
#  The Rothamsted Carbon Model: RothC
#  Developed by David Jenkinson and Kevin Coleman
#
#  INPUTS:
#
#  clay:  clay content of the soil (units: %)
#  depth: depth of topsoil (units: cm)
#  IOM: inert organic matter (t C /ha)
#  nsteps: number of timesteps
#
#  year:    year
#  month:   month (1-12)
#  modern:   %modern
#  TMP:      Air temperature (C)
#  Rain:     Rainfall (mm)
#  Evap:     open pan evaporation (mm)
#  C_inp:    carbon input to the soil each month (units: t C /ha)
#  FYM:      Farmyard manure input to the soil each month (units: t C /ha)
#  PC:       Plant cover (0 = no cover, 1 = covered by a crop)
#  DPM/RPM:  Ratio of DPM to RPM for carbon additions to the soil (units: none)
#
#  OUTPUTS:
#
#  All pools are carbon and not organic matter
#
#  DPM:   Decomposable Plant Material (units: t C /ha)
#  RPM:   Resistant Plant Material    (units: t C /ha)
#  Bio:   Microbial Biomass           (units: t C /ha)
#  Hum:   Humified Organic Matter     (units: t C /ha)
#  IOM:   Inert Organic Matter        (units: t C /ha)
#  SOC:   Soil Organic Matter / Total organic Matter (units: t C / ha)
#
#  DPM_Rage:   radiocarbon age of DPM
#  RPM_Rage:   radiocarbon age of RPM
#  Bio_Rage:   radiocarbon age of Bio
#  HUM_Rage:   radiocarbon age of Hum
#  Total_Rage: radiocarbon age of SOC (/ TOC)
#
#  SWC:       soil moisture deficit (mm per soil depth)
#  RM_TMP:    rate modifying fator for temperature (0.0 - ~5.0)
#  RM_Moist:  rate modifying fator for moisture (0.0 - 1.0)
#  RM_PC:     rate modifying fator for plant retainment (0.6 or 1.0)

######################################################################################################################
"""

from dataclasses import dataclass
from math import exp, log, isclose
from typing import Self, TypedDict
import logging


# =============================================================================
# Radiocarbon Constants
# =============================================================================

RADIO_HALFLIFE = 5568.0
"""Conventional half-life of radiocarbon-14 (years)."""

RADIO_MEAN_LIFETIME = 8035.0
"""Mean lifetime of radiocarbon-14, used in δ14C calculations (years)."""

IOM_INITIAL_AGE = 50000.0
"""Default initial age of inert organic matter, implying negligible 14C (years)."""


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


# =============================================================================
# Input Types
# =============================================================================


class RothCInputData(TypedDict):
    """Input data dictionary for the RothC model.

    All values are lists of monthly data, one entry per month.
    The lists should all have the same length.

    Attributes:
        t_tmp: Monthly mean air temperature (°C).
        t_rain: Monthly rainfall (mm).
        t_evap: Monthly open pan evaporation (mm).
        t_PC: Plant cover (0 = bare, 1 = covered).
        t_DPM_RPM: Ratio of DPM to RPM in plant inputs.
        t_C_Inp: Monthly plant carbon input (t C/ha).
        t_FYM_Inp: Monthly farmyard manure input (t C/ha).
        t_mod: Radiocarbon content as percentage modern carbon (%).
        t_year: Calendar year for each month.
        t_month: Month number (1-12) for each month.
    """

    t_tmp: list[float]
    t_rain: list[float]
    t_evap: list[float]
    t_PC: list[int]
    t_DPM_RPM: list[float]
    t_C_Inp: list[float]
    t_FYM_Inp: list[float]
    t_mod: list[float]
    t_year: list[int]
    t_month: list[int]


# =============================================================================
# Data Classes
# =============================================================================


@dataclass(eq=False)
class CarbonState:
    """Soil carbon pool state for the RothC model.

    Represents the state of all carbon pools and their radiocarbon ages
    at a given timestep.

    Attributes:
        dpm: Decomposable Plant Material (t C/ha).
        rpm: Resistant Plant Material (t C/ha).
        bio: Microbial Biomass (t C/ha).
        hum: Humified Organic Matter (t C/ha).
        iom: Inert Organic Matter (t C/ha).
        soc: Total Soil Organic Carbon (t C/ha).
        dpm_rc_age: Radiocarbon age of DPM pool (years).
        rpm_rc_age: Radiocarbon age of RPM pool (years).
        bio_rc_age: Radiocarbon age of BIO pool (years).
        hum_rc_age: Radiocarbon age of HUM pool (years).
        iom_age: Radiocarbon age of IOM pool (years).
        total_rc_age: Radiocarbon age of total SOC (years).
        swc: Soil water content/deficit (mm).
    """

    dpm: float
    rpm: float
    bio: float
    hum: float
    iom: float
    soc: float
    dpm_rc_age: float
    rpm_rc_age: float
    bio_rc_age: float
    hum_rc_age: float
    iom_age: float
    total_rc_age: float
    swc: float

    @classmethod
    def zero(cls) -> Self:
        """Create a CarbonState with all pools initialized to zero.

        The IOM age is set to the default value (50000 years).
        """
        return cls(
            dpm=0.0,
            rpm=0.0,
            bio=0.0,
            hum=0.0,
            iom=0.0,
            soc=0.0,
            dpm_rc_age=0.0,
            rpm_rc_age=0.0,
            bio_rc_age=0.0,
            hum_rc_age=0.0,
            iom_age=IOM_INITIAL_AGE,
            total_rc_age=0.0,
            swc=0.0,
        )

    def __eq__(
        self, other: object, rel_tol: float = 1e-10, abs_tol: float = 1e-10
    ) -> bool:
        if not isinstance(other, CarbonState):
            return False

        return all(
            [
                isclose(
                    getattr(self, attr),
                    getattr(other, attr),
                    rel_tol=rel_tol,
                    abs_tol=abs_tol,
                )
                for attr in [
                    "dpm",
                    "rpm",
                    "bio",
                    "iom",
                    "soc",
                    "dpm_rc_age",
                    "rpm_rc_age",
                    "bio_rc_age",
                    "hum_rc_age",
                    "iom_age",
                    "total_rc_age",
                    "swc",
                ]
            ]
        )


# =============================================================================
# Helper Functions
# =============================================================================


def temperature_rate_modifier(temp: float) -> float:
    """Calculate the rate modifying factor for temperature.

    Uses the Jenkinson equation to calculate the temperature rate modifier
    based on monthly mean air temperature.

    Args:
        temp: Monthly mean air temperature (°C).

    Returns:
        Rate modifying factor for temperature (typically 0.0 to ~5.0).
    """
    return (
        JENKINSON_A / (exp(JENKINSON_B / (temp + JENKINSON_C)) + 1.0)
        if temp > TEMP_MIN
        else 0.0
    )


def moisture_rate_modifier(
    rain: float,
    pevap: float,
    clay: float,
    depth: float,
    pc: bool,
    swc: float,
) -> tuple[float, float]:
    """Calculate the rate modifying factor for moisture.

    Calculates soil moisture deficit and derives a rate modifier based on
    the soil water balance, accounting for rainfall, evaporation, and
    plant cover.

    Args:
        rain: Monthly rainfall (mm).
        pevap: Open pan evaporation (mm).
        clay: Clay content of soil (%).
        depth: Depth of topsoil (cm).
        pc: Plant cover (False = no cover, True = covered).
        swc: Soil water content/deficit (mm).

    Returns:
        Tuple of (rate modifying factor for moisture, updated swc).
        rm_moist is typically between 0.2 and 1.0.
    """
    smd_max = -(SMD_COEFF_A + SMD_COEFF_B * clay - SMD_COEFF_C * (clay * clay))
    smd_max_adj = smd_max * depth / SMD_DEPTH_DIVISOR
    smd_1bar = SMD_1BAR_FRAC * smd_max_adj
    smd_bare = SMD_BARE_FRAC * smd_max_adj

    df = rain - EVAP_FACTOR * pevap

    min_swc_df = min(0.0, swc + df)

    if pc:
        swc_new = max(smd_max_adj, min_swc_df)
    else:
        min_smd_bare_swc = min(smd_bare, swc)
        swc_new = max(min_smd_bare_swc, min_swc_df)

    if swc_new > smd_1bar:
        rm_moist = RMF_MOIST_MAX
    else:
        rm_moist = RMF_MOIST_MIN + (RMF_MOIST_MAX - RMF_MOIST_MIN) * (
            smd_max_adj - swc_new
        ) / (smd_max_adj - smd_1bar)

    return rm_moist, swc_new


def plant_cover_rate_modifier(pc: bool) -> float:
    """Calculate the plant retainment modifying factor.

    Returns a reduced rate when the soil is covered by vegetation,
    representing reduced decomposition due to litter retention.

    Args:
        pc: Plant cover (False = no cover/bare soil, True = covered by crop).

    Returns:
        Rate modifying factor: 1.0 for bare soil, 0.6 for covered soil.
    """
    return RMF_PC_COVERED if pc else RMF_PC_BARE


def decompose_single_pool(
    pool: float, rate_k: float, rate_m: float, tstep: float
) -> tuple[float, float]:
    """Decompose a carbon pool using first-order decay kinetics.

    Args:
        pool: Carbon pool size (t C/ha).
        rate_k: Decomposition rate constant for this pool (per year).
        rate_m: Combined rate modifier (temperature × moisture × plant cover).
        tstep: Timestep (1/12 for monthly, 1/365 for daily).

    Returns:
        Tuple of (remaining_pool, decomposed_amount).
    """
    remaining = pool * exp(-rate_m * rate_k * tstep)
    decomposed = pool - remaining
    return remaining, decomposed


def partition_carbon_flows(decomposed: float, x: float) -> tuple[float, float, float]:
    """Partition decomposed carbon into CO2, BIO, and HUM fractions.

    The partitioning coefficient x depends on clay content and determines
    the proportion lost as CO2 vs incorporated into biomass/humus.

    Args:
        decomposed: Amount of carbon decomposed (t C/ha).
        x: Clay-dependent partitioning coefficient.

    Returns:
        Tuple of (co2, bio, hum) carbon amounts.
    """
    total = x + 1
    co2 = decomposed * (x / total)
    bio = decomposed * (FRAC_TO_BIO / total)
    hum = decomposed * (FRAC_TO_HUM / total)
    return co2, bio, hum


def calculate_radiocarbon_age(pool_new: float, ract_new: float, conr: float) -> float:
    """Calculate radiocarbon age from pool size and activity.

    Uses the radioactive decay equation inverted to solve for age:
    age = -ln(remaining/initial) / lambda
    where lambda = ln(2) / half_life

    Args:
        pool_new: Current pool size (t C/ha).
        ract_new: Radiocarbon activity (modern C equivalents).
        conr: Decay constant (1/mean_lifetime).

    Returns:
        Radiocarbon age in years.
    """
    return (
        log(pool_new / ract_new) / conr if pool_new > ZERO_THRESHOLD else ZERO_THRESHOLD
    )


def decompose_pools(
    state: CarbonState,
    modern_c: float,
    rate_m: float,
    clay: float,
    c_inp: float,
    fym_inp: float,
    dpm_rpm: float,
) -> CarbonState:
    """Calculate decomposition and radiocarbon age for soil carbon pools.

    Performs monthly carbon pool updates including: first-order decay
    kinetics, carbon flow between pools (DPM, RPM, BIO, HUM), CO2
    respiration, and radiocarbon age calculations.

    Args:
        state: Current carbon state (pools and ages).
        modern_c: Fraction of modern carbon (0.0 to 1.0).
        rate_m: Combined rate modifier.
        clay: Clay content of soil (%).
        c_inp: Plant carbon input (t C/ha).
        fym_inp: Farmyard manure carbon input (t C/ha).
        dpm_rpm: Ratio of DPM to RPM in plant inputs.

    Returns:
        Updated CarbonState.
    """
    dpm = state.dpm
    rpm = state.rpm
    bio = state.bio
    hum = state.hum
    iom = state.iom
    dpm_rc_age = state.dpm_rc_age
    rpm_rc_age = state.rpm_rc_age
    bio_rc_age = state.bio_rc_age
    hum_rc_age = state.hum_rc_age
    iom_age = state.iom_age

    conr = log(2.0) / RADIO_HALFLIFE

    tstep = 1.0 / MONTHS_PER_YEAR

    exc = exp(-conr * tstep)

    dpm1, dpm_d = decompose_single_pool(dpm, DPM_RATE, rate_m, tstep)
    rpm1, rpm_d = decompose_single_pool(rpm, RPM_RATE, rate_m, tstep)
    bio1, bio_d = decompose_single_pool(bio, BIO_RATE, rate_m, tstep)
    hum1, hum_d = decompose_single_pool(hum, HUM_RATE, rate_m, tstep)

    x = CLAW_A * (CLAW_B + CLAW_C * exp(-CLAW_D * clay))

    _, dpm_bio, dpm_hum = partition_carbon_flows(dpm_d, x)
    _, rpm_bio, rpm_hum = partition_carbon_flows(rpm_d, x)
    _, bio_bio, bio_hum = partition_carbon_flows(bio_d, x)
    _, hum_bio, hum_hum = partition_carbon_flows(hum_d, x)

    dpm_new = dpm1
    rpm_new = rpm1
    bio_new = bio1 + dpm_bio + rpm_bio + bio_bio + hum_bio
    hum_new = hum1 + dpm_hum + rpm_hum + bio_hum + hum_hum

    pi_c_dpm = dpm_rpm / (dpm_rpm + 1.0) * c_inp
    pi_c_rpm = 1.0 / (dpm_rpm + 1.0) * c_inp

    fym_c_dpm = FYM_FRAC_DPM * fym_inp
    fym_c_rpm = FYM_FRAC_RPM * fym_inp
    fym_c_hum = FYM_FRAC_HUM * fym_inp

    dpm_new = dpm_new + pi_c_dpm + fym_c_dpm
    rpm_new = rpm_new + pi_c_rpm + fym_c_rpm
    hum_new = hum_new + fym_c_hum

    dpm_ract = dpm1 * exp(-conr * dpm_rc_age)
    rpm_ract = rpm1 * exp(-conr * rpm_rc_age)

    bio_ract = bio1 * exp(-conr * bio_rc_age)
    dpm_bio_ract = dpm_bio * exp(-conr * dpm_rc_age)
    rpm_bio_ract = rpm_bio * exp(-conr * rpm_rc_age)
    bio_bio_ract = bio_bio * exp(-conr * bio_rc_age)
    hum_bio_ract = hum_bio * exp(-conr * hum_rc_age)

    hum_ract = hum1 * exp(-conr * hum_rc_age)
    dpm_hum_ract = dpm_hum * exp(-conr * dpm_rc_age)
    rpm_hum_ract = rpm_hum * exp(-conr * rpm_rc_age)
    bio_hum_ract = bio_hum * exp(-conr * bio_rc_age)
    hum_hum_ract = hum_hum * exp(-conr * hum_rc_age)

    iom_ract = iom * exp(-conr * iom_age)

    pi_dpm_ract = modern_c * pi_c_dpm
    pi_rpm_ract = modern_c * pi_c_rpm

    fym_dpm_ract = modern_c * fym_c_dpm
    fym_rpm_ract = modern_c * fym_c_rpm
    fym_hum_ract = modern_c * fym_c_hum

    dpm_ract_new = fym_dpm_ract + pi_dpm_ract + dpm_ract * exc
    rpm_ract_new = fym_rpm_ract + pi_rpm_ract + rpm_ract * exc

    bio_ract_new = (
        bio_ract + dpm_bio_ract + rpm_bio_ract + bio_bio_ract + hum_bio_ract
    ) * exc

    hum_ract_new = (
        fym_hum_ract
        + (hum_ract + dpm_hum_ract + rpm_hum_ract + bio_hum_ract + hum_hum_ract) * exc
    )

    soc_new = dpm_new + rpm_new + bio_new + hum_new + iom

    total_ract = dpm_ract_new + rpm_ract_new + bio_ract_new + hum_ract_new + iom_ract

    dpm_rc_age_new = calculate_radiocarbon_age(dpm_new, dpm_ract_new, conr)
    rpm_rc_age_new = calculate_radiocarbon_age(rpm_new, rpm_ract_new, conr)
    bio_rc_age_new = calculate_radiocarbon_age(bio_new, bio_ract_new, conr)
    hum_rc_age_new = calculate_radiocarbon_age(hum_new, hum_ract_new, conr)
    total_rc_age_new = calculate_radiocarbon_age(soc_new, total_ract, conr)

    return CarbonState(
        dpm=dpm_new,
        rpm=rpm_new,
        bio=bio_new,
        hum=hum_new,
        iom=iom,
        soc=soc_new,
        dpm_rc_age=dpm_rc_age_new,
        rpm_rc_age=rpm_rc_age_new,
        bio_rc_age=bio_rc_age_new,
        hum_rc_age=hum_rc_age_new,
        iom_age=iom_age,
        total_rc_age=total_rc_age_new,
        swc=state.swc,
    )


# =============================================================================
# RothC Class
# =============================================================================


class RothC:
    """Rothamsted Carbon Model.

    A class-based implementation of the RothC soil carbon model.
    """

    def __init__(self, clay: float, depth: float, iom: float):
        """Initialize the RothC model with soil parameters.

        Args:
            clay: Clay content of the soil (%).
            depth: Depth of topsoil (cm).
            iom: Inert organic matter (t C/ha).
        """
        self.clay = clay
        self.depth = depth
        self.iom = iom

    def run_timestep(
        self,
        state: CarbonState,
        temp: float,
        rain: float,
        pevap: float,
        pc: bool,
        dpm_rpm: float,
        c_inp: float,
        fym_inp: float,
        modern_c: float,
    ) -> CarbonState:
        """Run one timestep of the RothC model.

        Calculates rate modifying factors for temperature, moisture, and plant
        cover, then performs decomposition and radiocarbon age updates.

        Args:
            state: Current carbon state (pools and ages).
            temp: Monthly mean air temperature (°C).
            rain: Monthly rainfall (mm).
            pevap: Open pan evaporation (mm).
            pc: Plant cover (False = no cover, True = covered).
            dpm_rpm: Ratio of DPM to RPM in plant inputs.
            c_inp: Plant carbon input (t C/ha).
            fym_inp: Farmyard manure carbon input (t C/ha).
            modern_c: Fraction of modern carbon (0.0 to 1.0).

        Returns:
            Updated CarbonState.
        """
        rm_tmp = temperature_rate_modifier(temp)
        rm_moist, swc = moisture_rate_modifier(
            rain, pevap, self.clay, self.depth, pc, state.swc
        )
        rm_pc = plant_cover_rate_modifier(pc)

        rate_m = rm_tmp * rm_moist * rm_pc

        new_state = decompose_pools(
            state,
            modern_c,
            rate_m,
            self.clay,
            c_inp,
            fym_inp,
            dpm_rpm,
        )

        new_state.swc = swc

        return new_state

    def spin_up(self, data: RothCInputData) -> tuple[CarbonState, int]:
        """Spin up the RothC model to equilibrium.

        This method iteratively applies the same climate/input data until the
        annual change in total organic carbon falls below a threshold.
        Only the first 12 months of data are used, cycling through them repeatedly.

        Args:
            data: Dictionary containing monthly climate and input data.

        Returns:
            Tuple of (final carbon state at equilibrium, n_cycles).
        """
        state = CarbonState.zero()
        state.iom = self.iom

        months_per_cycle = len(data["t_tmp"])
        if months_per_cycle % 12 != 0:
            raise ValueError("Spin-up data should be a multiple of 12 months.")

        def data_iterator():
            return zip(
                data["t_tmp"],
                data["t_rain"],
                data["t_evap"],
                data["t_PC"],
                data["t_DPM_RPM"],
                data["t_C_Inp"],
                data["t_FYM_Inp"],
                data["t_mod"],
            )

        n_cycles = 0
        while True:
            toc_before_cycle = state.dpm + state.rpm + state.bio + state.hum

            # Cycle through entire series of spin-up driving data
            for (
                temp,
                rain,
                evap,
                pc,
                dpm_rpm,
                c_inp,
                fym_inp,
                modern_c,
            ) in data_iterator():
                state = self.run_timestep(
                    state,
                    temp,
                    rain,
                    evap,
                    pc,
                    dpm_rpm,
                    c_inp,
                    fym_inp,
                    modern_c / 100.0,
                )

            toc_after_cycle = state.dpm + state.rpm + state.bio + state.hum
            n_cycles += 1

            if abs(toc_after_cycle - toc_before_cycle) < EQUILIBRIUM_THRESHOLD:
                logging.info(
                    f"Spin-up converged after {n_cycles} cycles ({n_cycles * months_per_cycle} iterations)"
                )
                break

        return state, n_cycles

    def forward(
        self, state: CarbonState, data: RothCInputData
    ) -> tuple[CarbonState, dict[str, list]]:
        """Run the forward simulation from an initial state.

        Args:
            state: Initial carbon state (typically from spin_up).
            data: Dictionary containing monthly climate and input data.

        Returns:
            Tuple of (final carbon state, dict of monthly results where each key maps to a list of values).
        """

        def data_iterator():
            return zip(
                data["t_tmp"],
                data["t_rain"],
                data["t_evap"],
                data["t_PC"],
                data["t_DPM_RPM"],
                data["t_C_Inp"],
                data["t_FYM_Inp"],
                data["t_mod"],
            )

        (
            dpm_monthly,
            rpm_monthly,
            bio_monthly,
            hum_monthly,
            iom_monthly,
            soc_monthly,
            delta_monthly,
        ) = [], [], [], [], [], [], []

        for (
            temp,
            rain,
            evap,
            pc,
            dpm_rpm,
            c_inp,
            fym_inp,
            modern_c,
        ) in data_iterator():
            state = self.run_timestep(
                state,
                temp,
                rain,
                evap,
                pc,
                dpm_rpm,
                c_inp,
                fym_inp,
                modern_c / 100.0,
            )

            dpm_monthly.append(state.dpm)
            rpm_monthly.append(state.rpm)
            bio_monthly.append(state.bio)
            hum_monthly.append(state.hum)
            iom_monthly.append(state.iom)
            soc_monthly.append(state.soc)

            total_delta = (
                exp(-state.total_rc_age / RADIO_MEAN_LIFETIME) - 1.0
            ) * 1000.0
            delta_monthly.append(total_delta)

        results = dict(
            Year=data["t_year"],
            Month=data["t_month"],
            DPM_t_C_ha=dpm_monthly,
            RPM_t_C_ha=rpm_monthly,
            BIO_t_C_ha=bio_monthly,
            HUM_t_C_ha=hum_monthly,
            IOM_t_C_ha=iom_monthly,
            SOC_t_C_ha=soc_monthly,
            deltaC=delta_monthly,
        )

        return state, results

    def __call__(
        self, spinup_data: RothCInputData, forward_data: RothCInputData
    ) -> tuple[CarbonState, dict[str, list]]:
        """Run the full RothC simulation (spin-up + forward).

        Args:
            spinup_data: Dictionary containing monthly climate and input data for spin-up.
            forward_data: Dictionary containing monthly climate and input data for forward run.

        Returns:
            Tuple of (final carbon state, dict of monthly results where each key maps to a list of values).
        """
        state, _ = self.spin_up(spinup_data)
        return self.forward(state, forward_data)
