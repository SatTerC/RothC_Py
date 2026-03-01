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
from math import exp, log
from pathlib import Path
from typing import Self


# =============================================================================
# Radiocarbon Constants
# =============================================================================

# Radiocarbon half-life (years) - used in radiocarbon age calculations
RADIO_HALFLIFE = 5568.0

# Radiocarbon mean lifetime (years) - derived from half-life
# Mean lifetime = half-life / ln(2) ≈ 8033, using 8035.0 as in original
RADIO_MEAN_LIFETIME = 8035.0

# Initial radiocarbon age for IOM pool (years)
# Effectively infinite - IOM doesn't exchange with atmosphere
IOM_INITIAL_AGE = 50000.0


# =============================================================================
# Decomposition Rate Constants
# =============================================================================

# Threshold below which pool size is considered zero for age calculations
ZERO_THRESHOLD = 1e-8

# Decomposition rate constants for each pool (per year)
# DPM: Decomposable Plant Material - fast turnover
DPM_RATE = 10.0

# RPM: Resistant Plant Material - slow turnover
RPM_RATE = 0.3

# BIO: Microbial Biomass - intermediate turnover
BIO_RATE = 0.66

# HUM: Humified Organic Matter - very slow turnover
HUM_RATE = 0.02


# =============================================================================
# Carbon Flow Fractions
# =============================================================================

# Clay effect coefficients for CO2/BIO/HUM partitioning
# These determine the proportion of decomposed C lost as CO2 vs incorporated into BIO/HUM
CLAW_A = 1.67
CLAW_B = 1.85
CLAW_C = 1.60
CLAW_D = 0.0786

# Fraction of decomposed C flowing to microbial biomass (vs CO2)
FRAC_TO_BIO = 0.46

# Fraction of decomposed C flowing to humified organic matter (vs CO2)
FRAC_TO_HUM = 0.54

# Farmyard manure (FYM) composition fractions
FYM_FRAC_DPM = 0.49
FYM_FRAC_RPM = 0.49
FYM_FRAC_HUM = 0.02


# =============================================================================
# Moisture Rate Modifier Constants
# =============================================================================

# Maximum and minimum rate modifying factors for moisture
RMF_MOIST_MAX = 1.0
RMF_MOIST_MIN = 0.2

# Soil moisture deficit (SMD) formula coefficients
# SMDMax = -(20 + 1.3*clay - 0.01*clay^2)
SMD_COEFF_A = 20.0
SMD_COEFF_B = 1.3
SMD_COEFF_C = 0.01

# Depth adjustment divisor
SMD_DEPTH_DIVISOR = 23.0

# SMD at 1 bar pressure (fraction of SMDMax)
SMD_1BAR_FRAC = 0.444

# SMD for bare soil (fraction of SMDMax)
SMD_BARE_FRAC = 0.556

# Evaporation factor - proportion of PEVAP available for soil moisture
EVAP_FACTOR = 0.75


# =============================================================================
# Temperature Rate Modifier Constants
# =============================================================================

# Temperature below which decomposition rate is zero (°C)
TEMP_MIN = -5.0

# Jenkinson equation coefficients for temperature rate modifier
# RM_TMP = 47.91 / (exp(106.06 / (TEMP + 18.27)) + 1)
JENKINSON_A = 47.91
JENKINSON_B = 106.06
JENKINSON_C = 18.27


# =============================================================================
# Plant Cover Rate Modifier Constants
# =============================================================================

# Rate modifier for bare soil (no plant cover)
RMF_PC_BARE = 1.0

# Rate modifier for covered soil (plant cover reduces decomposition)
RMF_PC_COVERED = 0.6


# =============================================================================
# Simulation Constants
# =============================================================================

# Number of months per year
MONTHS_PER_YEAR = 12

# Equilibrium threshold - simulation stops when annual TOC change < this
EQUILIBRIUM_THRESHOLD = 1e-6


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class SoilParams:
    clay: float
    depth: float
    iom: float


@dataclass
class CarbonState:
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


def temperature_rate_modifier(temp: float, *, temp_min: float = TEMP_MIN) -> float:
    """Calculate the rate modifying factor for temperature.

    Uses the Jenkinson equation to calculate the temperature rate modifier
    based on monthly mean air temperature.

    Args:
        temp: Monthly mean air temperature (°C).
        temp_min: Temperature below which rate is zero (default: TEMP_MIN).

    Returns:
        Rate modifying factor for temperature (typically 0.0 to ~5.0).
    """
    if temp < temp_min:
        rm_tmp = 0.0
    else:
        rm_tmp = JENKINSON_A / (exp(JENKINSON_B / (temp + JENKINSON_C)) + 1.0)

    return rm_tmp


def moisture_rate_modifier(
    rain: float,
    pevap: float,
    clay: float,
    depth: float,
    pc: bool,
    swc: float,
    *,
    rmf_max: float = RMF_MOIST_MAX,
    rmf_min: float = RMF_MOIST_MIN,
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
        rmf_max: Maximum rate modifying factor (default: RMF_MOIST_MAX).
        rmf_min: Minimum rate modifying factor (default: RMF_MOIST_MIN).

    Returns:
        Tuple of (rate modifying factor for moisture, updated swc).
        rm_moist is typically between 0.2 and 1.0.
    """
    # calc soil water functions properties
    smd_max = -(SMD_COEFF_A + SMD_COEFF_B * clay - SMD_COEFF_C * (clay * clay))
    smd_max_adj = smd_max * depth / SMD_DEPTH_DIVISOR
    smd_1bar = SMD_1BAR_FRAC * smd_max_adj
    smd_bare = SMD_BARE_FRAC * smd_max_adj

    df = rain - EVAP_FACTOR * pevap

    min_swc_df = min(0.0, swc + df)
    min_smd_bare_swc = min(smd_bare, swc)

    if pc:
        swc_new = max(smd_max_adj, min_swc_df)
    else:
        swc_new = max(min_smd_bare_swc, min_swc_df)

    if swc_new > smd_1bar:
        rm_moist = 1.0
    else:
        rm_moist = rmf_min + (rmf_max - rmf_min) * (smd_max_adj - swc_new) / (
            smd_max_adj - smd_1bar
        )

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
    if not pc:
        rm_pc = RMF_PC_BARE
    else:
        rm_pc = RMF_PC_COVERED

    return rm_pc


def _decompose_single_pool(
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


def _partition_carbon_flows(decomposed: float, x: float) -> tuple[float, float, float]:
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


def _calculate_radiocarbon_age(pool_new: float, ract_new: float, conr: float) -> float:
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
    if pool_new <= ZERO_THRESHOLD:
        return ZERO_THRESHOLD
    return log(pool_new / ract_new) / conr


def decompose_pools(
    state: CarbonState,
    modern_c: float,
    rate_m: float,
    clay: float,
    c_inp: float,
    fym_inp: float,
    dpm_rpm: float,
    *,
    dpm_k: float = DPM_RATE,
    rpm_k: float = RPM_RATE,
    bio_k: float = BIO_RATE,
    hum_k: float = HUM_RATE,
    radio_halflife: float = RADIO_HALFLIFE,
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
        dpm_k: Decomposition rate constant for DPM (default: DPM_RATE).
        rpm_k: Decomposition rate constant for RPM (default: RPM_RATE).
        bio_k: Decomposition rate constant for BIO (default: BIO_RATE).
        hum_k: Decomposition rate constant for HUM (default: HUM_RATE).
        radio_halflife: Radiocarbon half-life in years (default: RADIO_HALFLIFE).

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

    conr = log(2.0) / radio_halflife

    tstep = 1.0 / MONTHS_PER_YEAR

    exc = exp(-conr * tstep)

    # decomposition
    dpm1, dpm_d = _decompose_single_pool(dpm, dpm_k, rate_m, tstep)
    rpm1, rpm_d = _decompose_single_pool(rpm, rpm_k, rate_m, tstep)
    bio1, bio_d = _decompose_single_pool(bio, bio_k, rate_m, tstep)
    hum1, hum_d = _decompose_single_pool(hum, hum_k, rate_m, tstep)

    x = CLAW_A * (CLAW_B + CLAW_C * exp(-CLAW_D * clay))

    # proportion C from each pool into CO2, BIO and HUM
    dpm_co2, dpm_bio, dpm_hum = _partition_carbon_flows(dpm_d, x)
    rpm_co2, rpm_bio, rpm_hum = _partition_carbon_flows(rpm_d, x)
    bio_co2, bio_bio, bio_hum = _partition_carbon_flows(bio_d, x)
    hum_co2, hum_bio, hum_hum = _partition_carbon_flows(hum_d, x)

    # update C pools
    dpm_new = dpm1
    rpm_new = rpm1
    bio_new = bio1 + dpm_bio + rpm_bio + bio_bio + hum_bio
    hum_new = hum1 + dpm_hum + rpm_hum + bio_hum + hum_hum

    # split plant C to DPM and RPM
    pi_c_dpm = dpm_rpm / (dpm_rpm + 1.0) * c_inp
    pi_c_rpm = 1.0 / (dpm_rpm + 1.0) * c_inp

    # split FYM C to DPM, RPM and HUM
    fym_c_dpm = FYM_FRAC_DPM * fym_inp
    fym_c_rpm = FYM_FRAC_RPM * fym_inp
    fym_c_hum = FYM_FRAC_HUM * fym_inp

    # add Plant C and FYM_C to DPM, RPM and HUM
    dpm_new = dpm_new + pi_c_dpm + fym_c_dpm
    rpm_new = rpm_new + pi_c_rpm + fym_c_rpm
    hum_new = hum_new + fym_c_hum

    # calc new ract of each pool
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

    # assign new C from plant and FYM the correct age
    pi_dpm_ract = modern_c * pi_c_dpm
    pi_rpm_ract = modern_c * pi_c_rpm

    fym_dpm_ract = modern_c * fym_c_dpm
    fym_rpm_ract = modern_c * fym_c_rpm
    fym_hum_ract = modern_c * fym_c_hum

    # update ract for each pool
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

    # calculate new radiocarbon age for each pool
    dpm_rc_age_new = _calculate_radiocarbon_age(dpm_new, dpm_ract_new, conr)
    rpm_rc_age_new = _calculate_radiocarbon_age(rpm_new, rpm_ract_new, conr)
    bio_rc_age_new = _calculate_radiocarbon_age(bio_new, bio_ract_new, conr)
    hum_rc_age_new = _calculate_radiocarbon_age(hum_new, hum_ract_new, conr)
    total_rc_age_new = _calculate_radiocarbon_age(soc_new, total_ract, conr)

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


def run_rothc_timestep(
    state: CarbonState,
    soil: SoilParams,
    temp: float,
    rain: float,
    pevap: float,
    pc: bool,
    dpm_rpm: float,
    c_inp: float,
    fym_inp: float,
    modern_c: float,
) -> CarbonState:
    """Run one timestep of the RothC carbon model.

    Calculates rate modifying factors for temperature, moisture, and plant
    cover, then performs decomposition and radiocarbon age updates.

    Args:
        state: Current carbon state (pools and ages).
        soil: Soil parameters (clay, depth, iom).
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
        rain, pevap, soil.clay, soil.depth, pc, state.swc
    )
    rm_pc = plant_cover_rate_modifier(pc)

    rate_m = rm_tmp * rm_moist * rm_pc

    new_state = decompose_pools(
        state,
        modern_c,
        rate_m,
        soil.clay,
        c_inp,
        fym_inp,
        dpm_rpm,
    )

    new_state.swc = swc

    return new_state


def spin_up_to_equilibrium(
    data: dict[str, list],
    state: CarbonState,
    soil: SoilParams,
) -> tuple[CarbonState, int, int]:
    """Spin up the RothC model to equilibrium using an acceleration technique.

    This function implements the RothC spin-up procedure, which accelerates the
    model to steady state by repeatedly cycling through a fixed period of climate
    and input data (typically one year) until the annual change in total organic
    carbon falls below a threshold.

    The acceleration technique avoids the need to simulate thousands of years of
    historical carbon inputs. Instead, the model runs much faster by iteratively
    applying the same climate/input year until equilibrium is reached.

    Args:
        data: Dictionary containing monthly climate and input data. Must have keys:
            t_tmp, t_rain, t_evap, t_PC, t_DPM_RPM, t_C_Inp, t_FYM_Inp, t_mod
        state: Initial carbon state.
        soil: Soil parameters (clay, depth, iom).

    Returns:
        Tuple of (final state, n_cycles, n_iterations).
    """
    toc_prev = 0.0
    n_cycles = -1
    k = -1
    total_iterations = 0

    while True:
        k = k + 1
        total_iterations = total_iterations + 1
        n_cycles = n_cycles + 1

        if k == MONTHS_PER_YEAR:
            k = 0

        temp = data["t_tmp"][k]
        rain = data["t_rain"][k]
        pevap = data["t_evap"][k]

        pc = bool(data["t_PC"][k])
        dpm_rpm = data["t_DPM_RPM"][k]

        c_inp = data["t_C_Inp"][k]
        fym_inp = data["t_FYM_Inp"][k]

        modern_c = data["t_mod"][k] / 100.0

        state = run_rothc_timestep(
            state,
            soil,
            temp,
            rain,
            pevap,
            pc,
            dpm_rpm,
            c_inp,
            fym_inp,
            modern_c,
        )

        if (k + 1) % MONTHS_PER_YEAR == 0:
            toc_curr = state.dpm + state.rpm + state.bio + state.hum
            if abs(toc_curr - toc_prev) < EQUILIBRIUM_THRESHOLD:
                break
            toc_prev = toc_curr

    return state, n_cycles, total_iterations


def run_simulation(
    data: dict[str, list],
    soil: SoilParams,
    nsteps: int,
) -> tuple[list[dict], list[dict]]:
    """Run the RothC model simulation.

    Args:
        data: Dictionary containing monthly climate and input data.
        soil: Soil parameters (clay, depth, iom).
        nsteps: Number of timesteps in the input data.

    Returns:
        Tuple of (year_results, month_results), each a list of dicts.
    """
    state = CarbonState.zero()
    state.iom = soil.iom

    print(0, state.dpm, state.rpm, state.bio, state.hum, state.iom, state.soc)

    (
        state,
        j,
        _,
    ) = spin_up_to_equilibrium(
        data,
        state,
        soil,
    )

    total_delta = (exp(-state.total_rc_age / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0
    print(
        j, state.dpm, state.rpm, state.bio, state.hum, state.iom, state.soc, total_delta
    )

    year_results = [
        {
            "Year": 1,
            "Month": j + 1,
            "DPM_t_C_ha": state.dpm,
            "RPM_t_C_ha": state.rpm,
            "BIO_t_C_ha": state.bio,
            "HUM_t_C_ha": state.hum,
            "IOM_t_C_ha": state.iom,
            "SOC_t_C_ha": state.soc,
            "deltaC": total_delta,
        }
    ]

    month_results = []

    for i in range(MONTHS_PER_YEAR, nsteps):
        temp = data["t_tmp"][i]
        rain = data["t_rain"][i]
        pevap = data["t_evap"][i]

        pc = bool(data["t_PC"][i])
        dpm_rpm = data["t_DPM_RPM"][i]

        c_inp = data["t_C_Inp"][i]
        fym_inp = data["t_FYM_Inp"][i]

        modern_c = data["t_mod"][i] / 100.0

        state = run_rothc_timestep(
            state,
            soil,
            temp,
            rain,
            pevap,
            pc,
            dpm_rpm,
            c_inp,
            fym_inp,
            modern_c,
        )

        total_delta = (exp(-state.total_rc_age / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0

        print(
            c_inp,
            fym_inp,
            temp,
            rain,
            pevap,
            state.swc,
            pc,
            state.dpm,
            state.rpm,
            state.bio,
            state.hum,
            state.iom,
            state.soc,
        )

        month_results.append(
            {
                "Year": data["t_year"][i],
                "Month": data["t_month"][i],
                "DPM_t_C_ha": state.dpm,
                "RPM_t_C_ha": state.rpm,
                "BIO_t_C_ha": state.bio,
                "HUM_t_C_ha": state.hum,
                "IOM_t_C_ha": state.iom,
                "SOC_t_C_ha": state.soc,
                "deltaC": total_delta,
            }
        )

        if data["t_month"][i] == MONTHS_PER_YEAR:
            year_results.append(
                {
                    "Year": data["t_year"][i],
                    "Month": data["t_month"][i],
                    "DPM_t_C_ha": state.dpm,
                    "RPM_t_C_ha": state.rpm,
                    "BIO_t_C_ha": state.bio,
                    "HUM_t_C_ha": state.hum,
                    "IOM_t_C_ha": state.iom,
                    "SOC_t_C_ha": state.soc,
                    "deltaC": total_delta,
                }
            )
            print(
                i,
                state.dpm,
                state.rpm,
                state.bio,
                state.hum,
                state.iom,
                state.soc,
                total_delta,
            )

    return year_results, month_results


def main(input_path: Path | str, output_dir: Path | str) -> None:
    """Run the RothC carbon model.

    Args:
        input_path: Path to the input data file.
        output_dir: Directory where output CSV files will be written.
    """
    import pandas as pd

    input_path = Path(input_path)
    output_dir = Path(output_dir)

    df_head = pd.read_csv(
        input_path,
        skiprows=3,
        header=0,
        nrows=1,
        index_col=None,
        sep=r"\s+",
    )
    clay = df_head.loc[0, "clay"]
    depth = df_head.loc[0, "depth"]
    iom = float(df_head.loc[0, "iom"])
    nsteps = df_head.loc[0, "nsteps"]
    soil = SoilParams(clay=clay, depth=depth, iom=iom)
    df = pd.read_csv(input_path, skiprows=6, header=0, index_col=None, sep=r"\s+")
    print(df)
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

    year_results, month_results = run_simulation(data, soil, nsteps)

    output_years = pd.DataFrame(year_results)
    output_months = pd.DataFrame(month_results)

    output_years.to_csv(output_dir / "year_results.csv", index=False)
    output_months.to_csv(output_dir / "month_results.csv", index=False)


if __name__ == "__main__":
    data_dir = Path(__file__).parent / "data"
    input_path = data_dir / "example_inputs.dat"
    output_dir = Path.cwd()
    main(input_path, output_dir)
