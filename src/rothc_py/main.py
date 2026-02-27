######################################################################################################
# program RothC_Python
from pathlib import Path

import numpy as np
import pandas as pd

data_dir = Path(__file__).parent / "data"
input_file = data_dir / "example_inputs.dat"

print(Path.cwd())

# set initial pool values
DPM = [0.0]
RPM = [0.0]
BIO = [0.0]
HUM = [0.0]
SOC = [0.0]

DPM_Rage = [0.0]
RPM_Rage = [0.0]
BIO_Rage = [0.0]
HUM_Rage = [0.0]
IOM_Rage = [50000.0]

# set initial soil water content (deficit)
SWC = [0.0]
TOC1 = 0.0

# read in RothC input data file
df_head = pd.read_csv(
    input_file,
    skiprows=3,
    header=0,
    nrows=1,
    index_col=None,
    delim_whitespace=True,
)
clay = df_head.loc[0, "clay"]
depth = df_head.loc[0, "depth"]
IOM = [df_head.loc[0, "iom"]]
nsteps = df_head.loc[0, "nsteps"]
df = pd.read_csv(
    input_file, skiprows=6, header=0, index_col=None, delim_whitespace=True
)
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
    "t_pc",
    "t_dpm_rpm",
]


# run rothc to equilibrium
k = -1
j = -1

soc[0] = dpm[0] + rpm[0] + bio[0] + hum[0] + iom[0]

print(j, dpm[0], rpm[0], bio[0], hum[0], iom[0], soc[0])

timefact = 12

test = 100.0
while test > 1e-6:
    k = k + 1
    j = j + 1

    if k == timefact:
        k = 0

    temp = df.t_tmp[k]
    rain = df.t_rain[k]
    pevap = df.t_evap[k]

    pc = df.t_pc[k]
    dpm_rpm = df.t_dpm_rpm[k]

    c_inp = df.t_c_inp[k]
    fym_inp = df.t_fym_inp[k]

    modernc = df.t_mod[k] / 100.0

    total_rage = [0.0]

    rothc(
        timefact,
        dpm,
        rpm,
        bio,
        hum,
        iom,
        soc,
        dpm_rage,
        rpm_rage,
        bio_rage,
        hum_rage,
        total_rage,
        modernc,
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
    if np.mod(k + 1, timefact) == 0:
        toc0 = toc1
        toc1 = dpm[0] + rpm[0] + bio[0] + hum[0]
        test = abs(toc1 - toc0)

total_delta = (np.exp(-total_rage[0] / 8035.0) - 1.0) * 1000.0

print(j, dpm[0], rpm[0], bio[0], hum[0], iom[0], soc[0], total_delta)

year_list = [[1, j + 1, dpm[0], rpm[0], bio[0], hum[0], iom[0], soc[0], total_delta[0]]]

month_list = []

for i in range(timefact, nsteps):
    temp = df.t_tmp[i]
    rain = df.t_rain[i]
    pevap = df.t_evap[i]

    pc = df.t_pc[i]
    dpm_rpm = df.t_dpm_rpm[i]

    c_inp = df.t_c_inp[i]
    fym_inp = df.t_fym_inp[i]

    modernc = df.t_mod[i] / 100.0

    rothc(
        timefact,
        dpm,
        rpm,
        bio,
        hum,
        iom,
        soc,
        dpm_rage,
        rpm_rage,
        bio_rage,
        hum_rage,
        total_rage,
        modernc,
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

    total_delta = (np.exp(-total_rage[0] / 8035.0) - 1.0) * 1000.0

    print(
        c_inp,
        fym_inp,
        temp,
        rain,
        pevap,
        swc[0],
        pc,
        dpm[0],
        rpm[0],
        bio[0],
        hum[0],
        iom[0],
        soc[0],
    )

    month_list.insert(
        i - timefact,
        [
            df.loc[i, "t_year"],
            df.loc[i, "t_month"],
            dpm[0],
            rpm[0],
            bio[0],
            hum[0],
            iom[0],
            soc[0],
            total_delta[0],
        ],
    )

    if df.t_month[i] == timefact:
        timefact_index = int(i / timefact)
        year_list.insert(
            timefact_index,
            [
                df.loc[i, "t_year"],
                df.loc[i, "t_month"],
                dpm[0],
                rpm[0],
                bio[0],
                hum[0],
                iom[0],
                soc[0],
                total_delta[0],
            ],
        )
        print(i, dpm, rpm, bio, hum, iom, soc, total_delta)

output_years = pd.dataframe(
    year_list,
    columns=[
        "year",
        "month",
        "dpm_t_c_ha",
        "rpm_t_c_ha",
        "bio_t_c_ha",
        "hum_t_c_ha",
        "iom_t_c_ha",
        "soc_t_c_ha",
        "deltac",
    ],
)
output_months = pd.dataframe(
    month_list,
    columns=[
        "year",
        "month",
        "dpm_t_c_ha",
        "rpm_t_c_ha",
        "bio_t_c_ha",
        "hum_t_c_ha",
        "iom_t_c_ha",
        "soc_t_c_ha",
        "deltac",
    ],
)

output_years.to_csv("year_results.csv", index=false)
output_months.to_csv("month_results.csv", index=false)
