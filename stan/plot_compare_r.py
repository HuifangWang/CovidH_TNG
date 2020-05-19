import numpy as np
import pylab as pl
import pandas as pd
import os
import datetime

de_sir_csv = '../data/processed/germany_imperial_college_sir.csv'
de_sir = pd.read_csv(de_sir_csv, parse_dates=['date'])
t0 = datetime.datetime(year=2020, month=2, day=15)
de_rci = de_sir.filter(regex='date|time.*repro.*$')
de_rci.set_index(de_rci.date, inplace=True)
de_rci.set_index((de_rci.index - t0).days, inplace=True)
print(de_rci.head())

de_pci = de_sir.filter(regex='date|pred.*infec.*cumulative$')
_, kmu, kqlo, kqhi = de_pci.columns
de_imu = np.r_[de_pci[kmu]]
de_isd = (np.r_[de_pci[kqhi]] - np.r_[de_pci[kqlo]]) / 4 # 95% is z -1.9 to 1.9

_, kmu, kqlo, kqhi = de_rci.columns
de_idx = np.r_[de_rci.index]
de_rmu = np.r_[de_rci[kmu]]
de_rsd = (np.r_[de_rci[kqhi]] - np.r_[de_rci[kqlo]]) / 4 # 95% is z
# TODO convert to mean & std for use in fitting


de_rki_csv = '../data/processed/germany_rki_nowcasting.csv'
de_klm_csv = '../data/external/kalman-filtered.csv'

de_rki = pd.read_csv(de_rki_csv, parse_dates=['date_of_disease_onset'])
de_rki.set_index(de_rki.date_of_disease_onset, inplace=True)
de_rki.set_index((de_rki.index - t0).days, inplace=True)
de_rki = de_rki.filter(regex='^repro.*$')
de_rki_idx = np.r_[de_rki.index][4:]
kmu, kqlo, kqhi = de_rki.columns
de_rki_rmu = np.r_[de_rki[kmu]][4:]
de_rki_rsd = (np.r_[de_rki[kqhi]] - np.r_[de_rki[kqlo]])[4:] / 4 # 95% is z

de_klm = pd.read_csv(de_klm_csv, parse_dates=['Date'])
de_klm = de_klm[de_klm["Country/Region"]=="Germany"]
de_klm = de_klm[de_klm.days_infectious==10]
de_klm.set_index(de_klm.Date, inplace=True)
de_klm.set_index((de_klm.index - t0).days, inplace=True)
print(de_klm.columns)

de_klm = de_klm.filter(regex='^R|^ci_95.*')
de_klm_idx = np.r_[de_klm.index]
kmu, kqlo, kqhi = de_klm.columns
de_klm_rmu = np.r_[de_klm[kmu]]
de_klm_rsd = (np.r_[de_klm[kqhi]] - np.r_[de_klm[kqlo]]) / 4 # 95% is z

pl.figure()
pl.plot(de_idx, de_rmu, 'b')
pl.plot(de_idx, de_rmu - de_rsd, 'b--')
pl.plot(de_idx, de_rmu + de_rsd, 'b--')
pl.plot(de_rki_idx, de_rki_rmu, 'g')
pl.plot(de_rki_idx, de_rki_rmu - de_rki_rsd, 'g--')
pl.plot(de_rki_idx, de_rki_rmu + de_rki_rsd, 'g--')
pl.plot(de_klm_idx, de_klm_rmu, 'r')
pl.plot(de_klm_idx, de_klm_rmu - de_klm_rsd, 'r--')
pl.plot(de_klm_idx, de_klm_rmu + de_klm_rsd, 'r--')
pl.show()
