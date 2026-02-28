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


def RMF_Tmp(TEMP: float) -> float:
    """Calculate the rate modifying factor for temperature.

    Uses the Jenkinson equation to calculate the temperature rate modifier
    based on monthly mean air temperature.

    Args:
        TEMP: Monthly mean air temperature (°C).

    Returns:
        Rate modifying factor for temperature (typically 0.0 to ~5.0).
    """
    if TEMP < TEMP_MIN:
        rm_tmp = 0.0
    else:
        rm_tmp = JENKINSON_A / (math.exp(JENKINSON_B / (TEMP + JENKINSON_C)) + 1.0)

    return rm_tmp


def RMF_Moist(
    RAIN: float,
    PEVAP: float,
    clay: float,
    depth: float,
    PC: bool,
    SWC: float,
) -> tuple[float, float]:
    """Calculate the rate modifying factor for moisture.

    Calculates soil moisture deficit and derives a rate modifier based on
    the soil water balance, accounting for rainfall, evaporation, and
    plant cover.

    Args:
        RAIN: Monthly rainfall (mm).
        PEVAP: Open pan evaporation (mm).
        clay: Clay content of soil (%).
        depth: Depth of topsoil (cm).
        PC: Plant cover (False = no cover, True = covered).
        SWC: Soil water content/deficit (mm).

    Returns:
        Tuple of (rate modifying factor for moisture, updated SWC).
        RM_Moist is typically between 0.2 and 1.0.
    """
    RMFMax = RMF_MOIST_MAX
    RMFMin = RMF_MOIST_MIN

    # calc soil water functions properties
    SMDMax = -(SMD_COEFF_A + SMD_COEFF_B * clay - SMD_COEFF_C * (clay * clay))
    SMDMaxAdj = SMDMax * depth / SMD_DEPTH_DIVISOR
    SMD1bar = SMD_1BAR_FRAC * SMDMaxAdj
    SMDBare = SMD_BARE_FRAC * SMDMaxAdj

    DF = RAIN - EVAP_FACTOR * PEVAP

    minSWCDF = min(0.0, SWC + DF)
    minSMDBareSWC = min(SMDBare, SWC)

    if PC:
        SWC_new = max(SMDMaxAdj, minSWCDF)
    else:
        SWC_new = max(minSMDBareSWC, minSWCDF)

    if SWC_new > SMD1bar:
        RM_Moist = 1.0
    else:
        RM_Moist = RMFMin + (RMFMax - RMFMin) * (SMDMaxAdj - SWC_new) / (
            SMDMaxAdj - SMD1bar
        )

    return RM_Moist, SWC_new


def RMF_PC(PC: float) -> float:
    """Calculate the plant retainment modifying factor.

    Returns a reduced rate when the soil is covered by vegetation,
    representing reduced decomposition due to litter retention.

    Args:
        PC: Plant cover (False = no cover/bare soil, True = covered by crop).

    Returns:
        Rate modifying factor: 1.0 for bare soil, 0.6 for covered soil.
    """
    if not PC:
        rm_pc = RMF_PC_BARE
    else:
        rm_pc = RMF_PC_COVERED

    return rm_pc


def decomp(
    timeFact: float,
    DPM: float,
    RPM: float,
    BIO: float,
    HUM: float,
    IOM: float,
    SOC: float,
    DPM_Rage: float,
    RPM_Rage: float,
    BIO_Rage: float,
    HUM_Rage: float,
    IOM_Rage: float,
    Total_Rage: float,
    modernC: float,
    RateM: float,
    clay: float,
    C_Inp: float,
    FYM_Inp: float,
    DPM_RPM: float,
) -> tuple[
    float, float, float, float, float, float, float, float, float, float, float, float
]:
    """Calculate decomposition and radiocarbon age for soil carbon pools.

    Performs monthly carbon pool updates including: first-order decay
    kinetics, carbon flow between pools (DPM, RPM, BIO, HUM), CO2
    respiration, and radiocarbon age calculations.

    Args:
        timeFact: Timestep factor (12 for monthly).
        DPM: Decomposable Plant Material pool (t C/ha).
        RPM: Resistant Plant Material pool (t C/ha).
        BIO: Microbial Biomass pool (t C/ha).
        HUM: Humified Organic Matter pool (t C/ha).
        IOM: Inert Organic Matter pool (t C/ha).
        SOC: Total Soil Organic Carbon (t C/ha).
        DPM_Rage: Radiocarbon age of DPM pool (years).
        RPM_Rage: Radiocarbon age of RPM pool (years).
        BIO_Rage: Radiocarbon age of BIO pool (years).
        HUM_Rage: Radiocarbon age of HUM pool (years).
        IOM_Rage: Radiocarbon age of IOM pool (years).
        Total_Rage: Radiocarbon age of total SOC (years).
        modernC: Fraction of modern carbon (0.0 to 1.0).
        clay: Clay content of soil (%).
        C_Inp: Plant carbon input (t C/ha).
        FYM_Inp: Farmyard manure carbon input (t C/ha).
        DPM_RPM: Ratio of DPM to RPM in plant inputs.

    Returns:
        Tuple of (DPM, RPM, BIO, HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, IOM_Rage, Total_Rage).
    """
    zero = ZERO_THRESHOLD
    # rate constant are params so don't need to be passed
    DPM_k = DPM_RATE
    RPM_k = RPM_RATE
    BIO_k = BIO_RATE
    HUM_k = HUM_RATE

    conr = math.log(2.0) / RADIO_HALFLIFE

    tstep = 1.0 / timeFact  # monthly 1/12, or daily 1/365

    exc = math.exp(-conr * tstep)

    # decomposition
    DPM1 = DPM * math.exp(-RateM * DPM_k * tstep)
    RPM1 = RPM * math.exp(-RateM * RPM_k * tstep)
    BIO1 = BIO * math.exp(-RateM * BIO_k * tstep)
    HUM1 = HUM * math.exp(-RateM * HUM_k * tstep)

    DPM_d = DPM - DPM1
    RPM_d = RPM - RPM1
    BIO_d = BIO - BIO1
    HUM_d = HUM - HUM1

    x = CLAW_A * (CLAW_B + CLAW_C * math.exp(-CLAW_D * clay))

    # proportion C from each pool into CO2, BIO and HUM
    DPM_co2 = DPM_d * (x / (x + 1))
    DPM_BIO = DPM_d * (FRAC_TO_BIO / (x + 1))
    DPM_HUM = DPM_d * (FRAC_TO_HUM / (x + 1))

    RPM_co2 = RPM_d * (x / (x + 1))
    RPM_BIO = RPM_d * (FRAC_TO_BIO / (x + 1))
    RPM_HUM = RPM_d * (FRAC_TO_HUM / (x + 1))

    BIO_co2 = BIO_d * (x / (x + 1))
    BIO_BIO = BIO_d * (FRAC_TO_BIO / (x + 1))
    BIO_HUM = BIO_d * (FRAC_TO_HUM / (x + 1))

    HUM_co2 = HUM_d * (x / (x + 1))
    HUM_BIO = HUM_d * (FRAC_TO_BIO / (x + 1))
    HUM_HUM = HUM_d * (FRAC_TO_HUM / (x + 1))

    # update C pools
    DPM_new = DPM1
    RPM_new = RPM1
    BIO_new = BIO1 + DPM_BIO + RPM_BIO + BIO_BIO + HUM_BIO
    HUM_new = HUM1 + DPM_HUM + RPM_HUM + BIO_HUM + HUM_HUM

    # split plant C to DPM and RPM
    PI_C_DPM = DPM_RPM / (DPM_RPM + 1.0) * C_Inp
    PI_C_RPM = 1.0 / (DPM_RPM + 1.0) * C_Inp

    # split FYM C to DPM, RPM and HUM
    FYM_C_DPM = FYM_FRAC_DPM * FYM_Inp
    FYM_C_RPM = FYM_FRAC_RPM * FYM_Inp
    FYM_C_HUM = FYM_FRAC_HUM * FYM_Inp

    # add Plant C and FYM_C to DPM, RPM and HUM
    DPM_new = DPM_new + PI_C_DPM + FYM_C_DPM
    RPM_new = RPM_new + PI_C_RPM + FYM_C_RPM
    HUM_new = HUM_new + FYM_C_HUM

    # calc new ract of each pool
    DPM_Ract = DPM1 * math.exp(-conr * DPM_Rage)
    RPM_Ract = RPM1 * math.exp(-conr * RPM_Rage)

    BIO_Ract = BIO1 * math.exp(-conr * BIO_Rage)
    DPM_BIO_Ract = DPM_BIO * math.exp(-conr * DPM_Rage)
    RPM_BIO_Ract = RPM_BIO * math.exp(-conr * RPM_Rage)
    BIO_BIO_Ract = BIO_BIO * math.exp(-conr * BIO_Rage)
    HUM_BIO_Ract = HUM_BIO * math.exp(-conr * HUM_Rage)

    HUM_Ract = HUM1 * math.exp(-conr * HUM_Rage)
    DPM_HUM_Ract = DPM_HUM * math.exp(-conr * DPM_Rage)
    RPM_HUM_Ract = RPM_HUM * math.exp(-conr * RPM_Rage)
    BIO_HUM_Ract = BIO_HUM * math.exp(-conr * BIO_Rage)
    HUM_HUM_Ract = HUM_HUM * math.exp(-conr * HUM_Rage)

    IOM_Ract = IOM * math.exp(-conr * IOM_Rage)

    # assign new C from plant and FYM the correct age
    PI_DPM_Ract = modernC * PI_C_DPM
    PI_RPM_Ract = modernC * PI_C_RPM

    FYM_DPM_Ract = modernC * FYM_C_DPM
    FYM_RPM_Ract = modernC * FYM_C_RPM
    FYM_HUM_Ract = modernC * FYM_C_HUM

    # update ract for each pool
    DPM_Ract_new = FYM_DPM_Ract + PI_DPM_Ract + DPM_Ract * exc
    RPM_Ract_new = FYM_RPM_Ract + PI_RPM_Ract + RPM_Ract * exc

    BIO_Ract_new = (
        BIO_Ract + DPM_BIO_Ract + RPM_BIO_Ract + BIO_BIO_Ract + HUM_BIO_Ract
    ) * exc

    HUM_Ract_new = (
        FYM_HUM_Ract
        + (HUM_Ract + DPM_HUM_Ract + RPM_HUM_Ract + BIO_HUM_Ract + HUM_HUM_Ract) * exc
    )

    SOC_new = DPM_new + RPM_new + BIO_new + HUM_new + IOM

    Total_Ract = DPM_Ract_new + RPM_Ract_new + BIO_Ract_new + HUM_Ract_new + IOM_Ract

    # calculate rage of each pool.
    if DPM_new <= zero:
        DPM_Rage_new = zero
    else:
        DPM_Rage_new = (math.log(DPM_new / DPM_Ract_new)) / conr

    if RPM_new <= zero:
        RPM_Rage_new = zero
    else:
        RPM_Rage_new = (math.log(RPM_new / RPM_Ract_new)) / conr

    if BIO_new <= zero:
        BIO_Rage_new = zero
    else:
        BIO_Rage_new = (math.log(BIO_new / BIO_Ract_new)) / conr

    if HUM_new <= zero:
        HUM_Rage_new = zero
    else:
        HUM_Rage_new = (math.log(HUM_new / HUM_Ract_new)) / conr

    if SOC_new <= zero:
        Total_Rage_new = zero
    else:
        Total_Rage_new = (math.log(SOC_new / Total_Ract)) / conr

    return (
        DPM_new,
        RPM_new,
        BIO_new,
        HUM_new,
        IOM,
        SOC_new,
        DPM_Rage_new,
        RPM_Rage_new,
        BIO_Rage_new,
        HUM_Rage_new,
        IOM_Rage,
        Total_Rage_new,
    )


def RothC(
    timeFact: float,
    DPM: float,
    RPM: float,
    BIO: float,
    HUM: float,
    IOM: float,
    SOC: float,
    DPM_Rage: float,
    RPM_Rage: float,
    BIO_Rage: float,
    HUM_Rage: float,
    IOM_Rage: float,
    Total_Rage: float,
    modernC: float,
    clay: float,
    depth: float,
    TEMP: float,
    RAIN: float,
    PEVAP: float,
    PC: bool,
    DPM_RPM: float,
    C_Inp: float,
    FYM_Inp: float,
    SWC: float,
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
        timeFact: Timestep factor (12 for monthly).
        DPM: Decomposable Plant Material pool (t C/ha).
        RPM: Resistant Plant Material pool (t C/ha).
        BIO: Microbial Biomass pool (t C/ha).
        HUM: Humified Organic Matter pool (t C/ha).
        IOM: Inert Organic Matter pool (t C/ha).
        SOC: Total Soil Organic Carbon (t C/ha).
        DPM_Rage: Radiocarbon age of DPM pool (years).
        RPM_Rage: Radiocarbon age of RPM pool (years).
        BIO_Rage: Radiocarbon age of BIO pool (years).
        HUM_Rage: Radiocarbon age of HUM pool (years).
        IOM_Rage: Radiocarbon age of IOM pool (years).
        Total_Rage: Radiocarbon age of total SOC (years).
        modernC: Fraction of modern carbon (0.0 to 1.0).
        clay: Clay content of soil (%).
        depth: Depth of topsoil (cm).
        TEMP: Monthly mean air temperature (°C).
        RAIN: Monthly rainfall (mm).
        PEVAP: Open pan evaporation (mm).
        PC: Plant cover (False = no cover, True = covered).
        DPM_RPM: Ratio of DPM to RPM in plant inputs.
        C_Inp: Plant carbon input (t C/ha).
        FYM_Inp: Farmyard manure carbon input (t C/ha).
        SWC: Soil water content/deficit (mm).

    Returns:
        Tuple of (DPM, RPM, BIO, HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, IOM_Rage, Total_Rage, SWC).
    """
    # Calculate RMFs
    RM_TMP = RMF_Tmp(TEMP)
    RM_Moist, SWC = RMF_Moist(RAIN, PEVAP, clay, depth, PC, SWC)
    RM_PC = RMF_PC(PC)

    # Combine RMF's into one.
    RateM = RM_TMP * RM_Moist * RM_PC

    (
        DPM,
        RPM,
        BIO,
        HUM,
        IOM,
        SOC,
        DPM_Rage,
        RPM_Rage,
        BIO_Rage,
        HUM_Rage,
        IOM_Rage,
        Total_Rage,
    ) = decomp(
        timeFact,
        DPM,
        RPM,
        BIO,
        HUM,
        IOM,
        SOC,
        DPM_Rage,
        RPM_Rage,
        BIO_Rage,
        HUM_Rage,
        IOM_Rage,
        Total_Rage,
        modernC,
        RateM,
        clay,
        C_Inp,
        FYM_Inp,
        DPM_RPM,
    )

    return (
        DPM,
        RPM,
        BIO,
        HUM,
        IOM,
        SOC,
        DPM_Rage,
        RPM_Rage,
        BIO_Rage,
        HUM_Rage,
        IOM_Rage,
        Total_Rage,
        SWC,
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
    DPM = 0.0
    RPM = 0.0
    BIO = 0.0
    HUM = 0.0
    SOC = 0.0

    DPM_Rage = 0.0
    RPM_Rage = 0.0
    BIO_Rage = 0.0
    HUM_Rage = 0.0
    IOM_Rage = IOM_INITIAL_AGE

    # set initial soil water content (deficit)
    SWC = 0.0
    TOC1 = 0.0

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
    IOM = float(df_head.loc[0, "iom"])
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

    SOC = DPM + RPM + BIO + HUM + IOM

    print(j, DPM, RPM, BIO, HUM, IOM, SOC)

    timeFact = MONTHS_PER_YEAR

    test = 100.0
    while test > EQUILIBRIUM_THRESHOLD:
        k = k + 1
        j = j + 1

        if k == timeFact:
            k = 0

        TEMP = df.t_tmp[k]
        RAIN = df.t_rain[k]
        PEVAP = df.t_evap[k]

        PC = bool(df.t_PC[k])
        DPM_RPM = df.t_DPM_RPM[k]

        C_Inp = df.t_C_Inp[k]
        FYM_Inp = df.t_FYM_Inp[k]

        modernC = df.t_mod[k] / 100.0

        Total_Rage = 0.0

        (
            DPM,
            RPM,
            BIO,
            HUM,
            IOM,
            SOC,
            DPM_Rage,
            RPM_Rage,
            BIO_Rage,
            HUM_Rage,
            IOM_Rage,
            Total_Rage,
            SWC,
        ) = RothC(
            timeFact,
            DPM,
            RPM,
            BIO,
            HUM,
            IOM,
            SOC,
            DPM_Rage,
            RPM_Rage,
            BIO_Rage,
            HUM_Rage,
            IOM_Rage,
            Total_Rage,
            modernC,
            clay,
            depth,
            TEMP,
            RAIN,
            PEVAP,
            PC,
            DPM_RPM,
            C_Inp,
            FYM_Inp,
            SWC,
        )

        # each a year calculates the difference between previous year and current year (counter =12 monthly model)
        if (k + 1) % timeFact == 0:
            TOC0 = TOC1
            TOC1 = DPM + RPM + BIO + HUM
            test = abs(TOC1 - TOC0)

    Total_Delta = (math.exp(-Total_Rage / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0

    print(j, DPM, RPM, BIO, HUM, IOM, SOC, Total_Delta)

    year_list = [[1, j + 1, DPM, RPM, BIO, HUM, IOM, SOC, Total_Delta]]

    month_list = []

    for i in range(timeFact, nsteps):
        TEMP = df.t_tmp[i]
        RAIN = df.t_rain[i]
        PEVAP = df.t_evap[i]

        PC = bool(df.t_PC[i])
        DPM_RPM = df.t_DPM_RPM[i]

        C_Inp = df.t_C_Inp[i]
        FYM_Inp = df.t_FYM_Inp[i]

        modernC = df.t_mod[i] / 100.0

        (
            DPM,
            RPM,
            BIO,
            HUM,
            IOM,
            SOC,
            DPM_Rage,
            RPM_Rage,
            BIO_Rage,
            HUM_Rage,
            IOM_Rage,
            Total_Rage,
            SWC,
        ) = RothC(
            timeFact,
            DPM,
            RPM,
            BIO,
            HUM,
            IOM,
            SOC,
            DPM_Rage,
            RPM_Rage,
            BIO_Rage,
            HUM_Rage,
            IOM_Rage,
            Total_Rage,
            modernC,
            clay,
            depth,
            TEMP,
            RAIN,
            PEVAP,
            PC,
            DPM_RPM,
            C_Inp,
            FYM_Inp,
            SWC,
        )

        Total_Delta = (math.exp(-Total_Rage / RADIO_MEAN_LIFETIME) - 1.0) * 1000.0

        print(
            C_Inp,
            FYM_Inp,
            TEMP,
            RAIN,
            PEVAP,
            SWC,
            PC,
            DPM,
            RPM,
            BIO,
            HUM,
            IOM,
            SOC,
        )

        month_list.insert(
            i - timeFact,
            [
                df.loc[i, "t_year"],
                df.loc[i, "t_month"],
                DPM,
                RPM,
                BIO,
                HUM,
                IOM,
                SOC,
                Total_Delta,
            ],
        )

        if df.t_month[i] == timeFact:
            timeFact_index = int(i / timeFact)
            year_list.insert(
                timeFact_index,
                [
                    df.loc[i, "t_year"],
                    df.loc[i, "t_month"],
                    DPM,
                    RPM,
                    BIO,
                    HUM,
                    IOM,
                    SOC,
                    Total_Delta,
                ],
            )
            print(i, DPM, RPM, BIO, HUM, IOM, SOC, Total_Delta)

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
