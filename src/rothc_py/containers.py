"""Data containers for the RothC soil carbon model.

This module contains the core data structures used by the RothC model,
including input data types and the carbon state representation.
"""

from dataclasses import dataclass
from math import isclose
from typing import Self, TypedDict

from rothc_py.constants import IOM_INITIAL_AGE


# =============================================================================
# Input Types
# =============================================================================


class InputData(TypedDict):
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
