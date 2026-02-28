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

import math
from pathlib import Path

import pandas as pd


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


def rmf_tmp(temp: float) -> float:
    """Calculate the rate modifying factor for temperature.

    Uses the Jenkinson equation to calculate the temperature rate modifier
    based on monthly mean air temperature.

    Args:
        temp: Monthly mean air temperature (°C).

    Returns:
        Rate modifying factor for temperature (typically 0.0 to ~5.0).
    """
    if temp < TEMP_MIN:
        rm_tmp = 0.0
    else:
        rm_tmp = JENKINSON_A / (math.exp(JENKINSON_B / (temp + JENKINSON_C)) + 1.0)

    return rm_tmp


def rmf_moist(
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
    rmf_max = RMF_MOIST_MAX
    rmf_min = RMF_MOIST_MIN

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


def rmf_pc(pc: bool) -> float:
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


def decompose(
    time_step: float,
    dpm: float,
    rpm: float,
    bio: float,
    hum: float,
    iom: float,
    soc: float,
    dpm_rc_age: float,
    rpm_rc_age: float,
    bio_rc_age: float,
    hum_rc_age: float,
    iom_age: float,
    total_rc_age: float,
    modern_c: float,
    rate_m: float,
    clay: float,
    c_inp: float,
    fym_inp: float,
    dpm_rpm: float,
) -> tuple[
    float, float, float, float, float, float, float, float, float, float, float, float
]:
    """Calculate decomposition and radiocarbon age for soil carbon pools.

    Performs monthly carbon pool updates including: first-order decay
    kinetics, carbon flow between pools (DPM, RPM, BIO, HUM), CO2
    respiration, and radiocarbon age calculations.

    Args:
        time_step: Timestep factor (12 for monthly).
        dpm: Decomposable Plant Material pool (t C/ha).
        rpm: Resistant Plant Material pool (t C/ha).
        bio: Microbial Biomass pool (t C/ha).
        hum: Humified Organic Matter pool (t C/ha).
        iom: Inert Organic Matter pool (t C/ha).
        soc: Total Soil Organic Carbon (t C/ha).
        dpm_rc_age: Radiocarbon age of DPM pool (years).
        rpm_rc_age: Radiocarbon age of RPM pool (years).
        bio_rc_age: Radiocarbon age of BIO pool (years).
        hum_rc_age: Radiocarbon age of HUM pool (years).
        iom_age: Radiocarbon age of IOM pool (years).
        total_rc_age: Radiocarbon age of total SOC (years).
        modern_c: Fraction of modern carbon (0.0 to 1.0).
        clay: Clay content of soil (%).
        c_inp: Plant carbon input (t C/ha).
        fym_inp: Farmyard manure carbon input (t C/ha).
        dpm_rpm: Ratio of DPM to RPM in plant inputs.

    Returns:
        Tuple of (dpm, rpm, bio, hum, iom, soc, dpm_rc_age, rpm_rc_age, bio_rc_age, hum_rc_age, iom_age, total_rc_age).
    """
    zero = ZERO_THRESHOLD
    # rate constant are params so don't need to be passed
    dpm_k = DPM_RATE
    rpm_k = RPM_RATE
    bio_k = BIO_RATE
    hum_k = HUM_RATE

    conr = math.log(2.0) / RADIO_HALFLIFE

    tstep = 1.0 / time_step  # monthly 1/12, or daily 1/365

    exc = math.exp(-conr * tstep)

    # decomposition
    dpm1 = dpm * math.exp(-rate_m * dpm_k * tstep)
    rpm1 = rpm * math.exp(-rate_m * rpm_k * tstep)
    bio1 = bio * math.exp(-rate_m * bio_k * tstep)
    hum1 = hum * math.exp(-rate_m * hum_k * tstep)

    dpm_d = dpm - dpm1
    rpm_d = rpm - rpm1
    bio_d = bio - bio1
    hum_d = hum - hum1

    x = CLAW_A * (CLAW_B + CLAW_C * math.exp(-CLAW_D * clay))

    # proportion C from each pool into CO2, BIO and HUM
    dpm_co2 = dpm_d * (x / (x + 1))
    dpm_bio = dpm_d * (FRAC_TO_BIO / (x + 1))
    dpm_hum = dpm_d * (FRAC_TO_HUM / (x + 1))

    rpm_co2 = rpm_d * (x / (x + 1))
    rpm_bio = rpm_d * (FRAC_TO_BIO / (x + 1))
    rpm_hum = rpm_d * (FRAC_TO_HUM / (x + 1))

    bio_co2 = bio_d * (x / (x + 1))
    bio_bio = bio_d * (FRAC_TO_BIO / (x + 1))
    bio_hum = bio_d * (FRAC_TO_HUM / (x + 1))

    hum_co2 = hum_d * (x / (x + 1))
    hum_bio = hum_d * (FRAC_TO_BIO / (x + 1))
    hum_hum = hum_d * (FRAC_TO_HUM / (x + 1))

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
    dpm_ract = dpm1 * math.exp(-conr * dpm_rc_age)
    rpm_ract = rpm1 * math.exp(-conr * rpm_rc_age)

    bio_ract = bio1 * math.exp(-conr * bio_rc_age)
    dpm_bio_ract = dpm_bio * math.exp(-conr * dpm_rc_age)
    rpm_bio_ract = rpm_bio * math.exp(-conr * rpm_rc_age)
    bio_bio_ract = bio_bio * math.exp(-conr * bio_rc_age)
    hum_bio_ract = hum_bio * math.exp(-conr * hum_rc_age)

    hum_ract = hum1 * math.exp(-conr * hum_rc_age)
    dpm_hum_ract = dpm_hum * math.exp(-conr * dpm_rc_age)
    rpm_hum_ract = rpm_hum * math.exp(-conr * rpm_rc_age)
    bio_hum_ract = bio_hum * math.exp(-conr * bio_rc_age)
    hum_hum_ract = hum_hum * math.exp(-conr * hum_rc_age)

    iom_ract = iom * math.exp(-conr * iom_age)

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

    # calculate rage of each pool.
    if dpm_new <= zero:
        dpm_rc_age_new = zero
    else:
        dpm_rc_age_new = (math.log(dpm_new / dpm_ract_new)) / conr

    if rpm_new <= zero:
        rpm_rc_age_new = zero
    else:
        rpm_rc_age_new = (math.log(rpm_new / rpm_ract_new)) / conr

    if bio_new <= zero:
        bio_rc_age_new = zero
    else:
        bio_rc_age_new = (math.log(bio_new / bio_ract_new)) / conr

    if hum_new <= zero:
        hum_rc_age_new = zero
    else:
        hum_rc_age_new = (math.log(hum_new / hum_ract_new)) / conr

    if soc_new <= zero:
        total_rc_age_new = zero
    else:
        total_rc_age_new = (math.log(soc_new / total_ract)) / conr

    return (
        dpm_new,
        rpm_new,
        bio_new,
        hum_new,
        iom,
        soc_new,
        dpm_rc_age_new,
        rpm_rc_age_new,
        bio_rc_age_new,
        hum_rc_age_new,
        iom_age,
        total_rc_age_new,
    )


def run_rothc(
    time_step: float,
    dpm: float,
    rpm: float,
    bio: float,
    hum: float,
    iom: float,
    soc: float,
    dpm_rc_age: float,
    rpm_rc_age: float,
    bio_rc_age: float,
    hum_rc_age: float,
    iom_age: float,
    total_rc_age: float,
    modern_c: float,
    clay: float,
    depth: float,
    temp: float,
    rain: float,
    pevap: float,
    pc: bool,
    dpm_rpm: float,
    c_inp: float,
    fym_inp: float,
    swc: float,
) -> tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]:
    """Run one timestep of the RothC carbon model.

    Calculates rate modifying factors for temperature, moisture, and plant
    cover, then performs decomposition and radiocarbon age updates.

    Args:
        time_step: Timestep factor (12 for monthly).
        dpm: Decomposable Plant Material pool (t C/ha).
        rpm: Resistant Plant Material pool (t C/ha).
        bio: Microbial Biomass pool (t C/ha).
        hum: Humified Organic Matter pool (t C/ha).
        iom: Inert Organic Matter pool (t C/ha).
        soc: Total Soil Organic Carbon (t C/ha).
        dpm_rc_age: Radiocarbon age of DPM pool (years).
        rpm_rc_age: Radiocarbon age of RPM pool (years).
        bio_rc_age: Radiocarbon age of BIO pool (years).
        hum_rc_age: Radiocarbon age of HUM pool (years).
        iom_age: Radiocarbon age of IOM pool (years).
        total_rc_age: Radiocarbon age of total SOC (years).
        modern_c: Fraction of modern carbon (0.0 to 1.0).
        clay: Clay content of soil (%).
        depth: Depth of topsoil (cm).
        temp: Monthly mean air temperature (°C).
        rain: Monthly rainfall (mm).
        pevap: Open pan evaporation (mm).
        pc: Plant cover (False = no cover, True = covered).
        dpm_rpm: Ratio of DPM to RPM in plant inputs.
        c_inp: Plant carbon input (t C/ha).
        fym_inp: Farmyard manure carbon input (t C/ha).
        swc: Soil water content/deficit (mm).

    Returns:
        Tuple of (dpm, rpm, bio, hum, iom, soc, dpm_rc_age, rpm_rc_age, bio_rc_age, hum_rc_age, iom_age, total_rc_age, swc).
    """
    # Calculate RMFs
    rm_tmp = rmf_tmp(temp)
    rm_moist, swc = rmf_moist(rain, pevap, clay, depth, pc, swc)
    rm_pc = rmf_pc(pc)

    # Combine RMF's into one.
    rate_m = rm_tmp * rm_moist * rm_pc

    (
        dpm,
        rpm,
        bio,
        hum,
        iom,
        soc,
        dpm_rc_age,
        rpm_rc_age,
        bio_rc_age,
        hum_rc_age,
        iom_age,
        total_rc_age,
    ) = decompose(
        time_step,
        dpm,
        rpm,
        bio,
        hum,
        iom,
        soc,
        dpm_rc_age,
        rpm_rc_age,
        bio_rc_age,
        hum_rc_age,
        iom_age,
        total_rc_age,
        modern_c,
        rate_m,
        clay,
        c_inp,
        fym_inp,
        dpm_rpm,
    )

    return (
        dpm,
        rpm,
        bio,
        hum,
        iom,
        soc,
        dpm_rc_age,
        rpm_rc_age,
        bio_rc_age,
        hum_rc_age,
        iom_age,
        total_rc_age,
        swc,
    )


def main(input_path: Path | str, output_dir: Path | str) -> None:
    """Run the RothC carbon model.

    Args:
        input_path: Path to the input data file.
        output_dir: Directory where output CSV files will be written.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    ######################################################################################################
    # program RothC_Python

    # set initial pool values
    dpm = 0.0
    rpm = 0.0
    bio = 0.0
    hum = 0.0
    soc = 0.0

    dpm_rc_age = 0.0
    rpm_rc_age = 0.0
    bio_rc_age = 0.0
    hum_rc_age = 0.0
    iom_age = IOM_INITIAL_AGE

    # set initial soil water content (deficit)
    swc = 0.0
    toc1 = 0.0

    # read in RothC input data file
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

    # run RothC to equilibrium
    k = -1
    j = -1

    soc = dpm + rpm + bio + hum + iom

    print(j, dpm, rpm, bio, hum, iom, soc)

    time_step = MONTHS_PER_YEAR

    test = 100.0
    while test > EQUILIBRIUM_THRESHOLD:
        k = k + 1
        j = j + 1

        if k == time_step:
            k = 0

        temp = df.t_tmp[k]
        rain = df.t_rain[k]
        pevap = df.t_evap[k]

        pc = bool(df.t_PC[k])
        dpm_rpm = df.t_DPM_RPM[k]

        c_inp = df.t_C_Inp[k]
        fym_inp = df.t_FYM_Inp[k]

        modern_c = df.t_mod[k] / 100.0

        total_rc_age = 0.0

        (
            dpm,
            rpm,
            bio,
            hum,
            iom,
            soc,
            dpm_rc_age,
            rpm_rc_age,
            bio_rc_age,
            hum_rc_age,
            iom_age,
            total_rc_age,
            swc,
        ) = run_rothc(
            time_step,
            dpm,
            rpm,
            bio,
            hum,
            iom,
            soc,
            dpm_rc_age,
            rpm_rc_age,
            bio_rc_age,
            hum_rc_age,
            iom_age,
            total_rc_age,
            modern_c,
            clay,
            depth,
            temp,
            rain,
            pevap,
            pc,
            dpm_rpm,
            c_inp,
            fym_inp,
            swc,
        )

        # each a year calculates the difference between previous year and current year (counter =12 monthly model)
        if (k + 1) % time_step == 0:
            toc0 = toc1
            toc1 = dpm + rpm + bio + hum
            test = abs(toc1 - toc0)

    total_delta = (math.exp(-total_rc_age / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0

    print(j, dpm, rpm, bio, hum, iom, soc, total_delta)

    year_list = [[1, j + 1, dpm, rpm, bio, hum, iom, soc, total_delta]]

    month_list = []

    for i in range(time_step, nsteps):
        temp = df.t_tmp[i]
        rain = df.t_rain[i]
        pevap = df.t_evap[i]

        pc = bool(df.t_PC[i])
        dpm_rpm = df.t_DPM_RPM[i]

        c_inp = df.t_C_Inp[i]
        fym_inp = df.t_FYM_Inp[i]

        modern_c = df.t_mod[i] / 100.0

        (
            dpm,
            rpm,
            bio,
            hum,
            iom,
            soc,
            dpm_rc_age,
            rpm_rc_age,
            bio_rc_age,
            hum_rc_age,
            iom_age,
            total_rc_age,
            swc,
        ) = run_rothc(
            time_step,
            dpm,
            rpm,
            bio,
            hum,
            iom,
            soc,
            dpm_rc_age,
            rpm_rc_age,
            bio_rc_age,
            hum_rc_age,
            iom_age,
            total_rc_age,
            modern_c,
            clay,
            depth,
            temp,
            rain,
            pevap,
            pc,
            dpm_rpm,
            c_inp,
            fym_inp,
            swc,
        )

        total_delta = (math.exp(-total_rc_age / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0

        print(
            c_inp,
            fym_inp,
            temp,
            rain,
            pevap,
            swc,
            pc,
            dpm,
            rpm,
            bio,
            hum,
            iom,
            soc,
        )

        month_list.insert(
            i - time_step,
            [
                df.loc[i, "t_year"],
                df.loc[i, "t_month"],
                dpm,
                rpm,
                bio,
                hum,
                iom,
                soc,
                total_delta,
            ],
        )

        if df.t_month[i] == time_step:
            time_step_index = int(i / time_step)
            year_list.insert(
                time_step_index,
                [
                    df.loc[i, "t_year"],
                    df.loc[i, "t_month"],
                    dpm,
                    rpm,
                    bio,
                    hum,
                    iom,
                    soc,
                    total_delta,
                ],
            )
            print(i, dpm, rpm, bio, hum, iom, soc, total_delta)

    output_years = pd.DataFrame(
        year_list,
        columns=[
            "Year",
            "Month",
            "DPM_t_C_ha",
            "RPM_t_C_ha",
            "BIO_t_C_ha",
            "HUM_t_C_ha",
            "IOM_t_C_ha",
            "SOC_t_C_ha",
            "deltaC",
        ],
    )
    output_months = pd.DataFrame(
        month_list,
        columns=[
            "Year",
            "Month",
            "DPM_t_C_ha",
            "RPM_t_C_ha",
            "BIO_t_C_ha",
            "HUM_t_C_ha",
            "IOM_t_C_ha",
            "SOC_t_C_ha",
            "deltaC",
        ],
    )

    output_years.to_csv(output_dir / "year_results.csv", index=False)
    output_months.to_csv(output_dir / "month_results.csv", index=False)


if __name__ == "__main__":
    data_dir = Path(__file__).parent / "data"
    input_path = data_dir / "example_inputs.dat"
    output_dir = Path.cwd()
    main(input_path, output_dir)
