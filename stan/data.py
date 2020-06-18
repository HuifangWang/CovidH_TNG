import numpy as np
import pylab as pl
import pandas as pd
import datetime
import stanio

# start with March
t0 = datetime.datetime(year=2020, month=3, day=10)


def fix_t0(df, date_key):
    df = df.set_index(df[date_key])
    df = df.set_index((df.index - t0).days)
    return df


def load_rki_r():
    dk = "date_of_disease_onset"
    "Load Robert Koch Institute reproduction rate values."
    df = pd.read_csv(
        "../data/processed/germany_rki_nowcasting.csv", parse_dates=[dk]
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="^R_.*point")
    return df


def load_mob():
    dk = "date"
    df = pd.read_csv(
        "../data/processed/germany_google_mobility_report.csv",
        parse_dates=[dk],
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="retail")
    return df


def load_phi():
    dk = "date"
    df = pd.read_csv(
        "../data/processed/germany_cosmo_phi.csv", parse_dates=[dk]
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="phi")
    return df


def main():
    df_r = load_rki_r().loc[0:60]
    df_m = load_mob().loc[0:60]
    df_p = load_phi().loc[1:]  # not daily, but 8 weeks
    pl.plot(df_r.index, np.r_[df_r])
    pl.plot(df_m.index, np.r_[df_m] / 100)
    pl.plot(df_p.index, np.r_[df_p] / np.r_[df_p].ptp())
    pl.legend(["RKI R 7 day", "mobility/100%", "CosmoPhi"])
    pl.savefig("beta-data.png")
    stanio.rdump(
        "beta.R",
        {
            "nt": len(df_r.index),
            "r": np.r_[df_r][:, 0],
            "m": np.r_[df_m][:, 0] / 100,
            "iP": np.r_[df_p.index] + 1,
            "P": np.r_[df_p][:, 0],
        },
    )


if __name__ == "__main__":
    main()
