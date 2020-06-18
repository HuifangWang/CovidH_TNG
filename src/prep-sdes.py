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


def load_r0t():
    dk = "date"
    "Load Denmark Epiforecast reproduction rate values."
    df = pd.read_csv(
        "../data/processed/denmark_epiforecast_rt.csv", parse_dates=[dk]
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="^median$")
    return df


def load_mob():
    dk = "date"
    df = pd.read_csv(
        "../data/processed/denmark_google_mobility_report.csv",
        parse_dates=[dk],
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="retail")
    return df


def load_phi():
    dk = "date"
    df = pd.read_csv(
        "../data/processed/denmark_cosmo_phi.csv", parse_dates=[dk]
    )
    df = fix_t0(df, dk)
    df = df.filter(regex="phi")
    return df


def main():
    df_r = load_r0t().loc[0:60]
    df_m = load_mob().loc[0:60]
    df_p = load_phi().iloc[1:9]  # not daily, but 8 weeks
    pl.subplot(311)
    pl.plot(df_r.index, np.r_[df_r])
    pl.title("RKI R 7 day")
    pl.grid(1)
    pl.subplot(312)
    pl.plot(df_m.index, np.r_[df_m] / 100)
    pl.title("mobility/100%")
    pl.grid(1)
    pl.subplot(313)
    pl.plot(df_p.index, np.r_[df_p])
    pl.title("CosmoPhi")
    pl.grid(1)
    pl.tight_layout()
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
