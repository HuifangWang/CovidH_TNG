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

de_mob_csv = '../data/processed/germany_google_mobility_report.csv'
de_mob = pd.read_csv(de_mob_csv, parse_dates=['date'])
de_mob.set_index(de_mob.date, inplace=True)
de_mob.set_index((de_mob.index - t0).days, inplace=True)
de_mob = de_mob.filter(regex='^retail')
de_mobi = de_mob.index
de_mobp = np.r_[de_mob][:, 0]/100

de_phi_csv = '../data/processed/germany_cosmo_phi.csv'
de_phi = pd.read_csv(de_phi_csv, parse_dates=['date'])
de_phi.set_index(de_phi.date, inplace=True)
de_phi.set_index((de_phi.index - t0).days, inplace=True)
de_iphi = de_phi.index
de_vphi = np.r_[de_phi][:, 1]

pl.plot(de_idx, de_rmu, 'ko')
pl.plot(de_idx, de_rsd, 'kx')
pl.plot(de_idx, np.log(de_imu), 'go')
pl.plot(de_idx, np.log(de_isd), 'gx')
pl.plot(de_mobi, de_mobp, 'bo')
pl.plot(de_iphi, de_vphi, 'ro')
pl.savefig('german-data.png')
# pl.show()


import stanio
cs = '/Volumes/crypt/cmdstan-2.22.1'

stanio.compile_model(cs, 'model')

data = {
    'nii': 10,
    'nr': len(de_idx),
    'rid': de_idx + 1,
    'rmu': de_rmu,
    'rsd': de_rsd,
    'imu': de_imu,
    'isd': de_isd,
    'nm': len(de_mobi),
    'mid': np.r_[de_mobi] + 1,
    'mpr': de_mobp,
    'np': len(de_iphi),
    'pid': np.r_[de_iphi] + 1,
    'phi': de_vphi,
    'I0': 1e-6,
    }

init = {
    'gamma': 0.1,
    'Gamma': np.r_[0.5,3.0],
    'c': np.r_[0.5,5.0],
    'alpha': 0.06,
    'J': 4.0,
}
stanio.rdump('model.R', data)
stanio.rdump('init.R', init)
os.system('./model sample algorithm=fixed_param num_samples=1 init=init.R data file=model.R output file=fit.csv')
csv = stanio.parse_csv('fit.csv')
y = np.transpose(csv['yt_'][0].reshape((-1, 5, 2)), (1,0,2))
names = r'$S(t)$ $I(t)$ $\beta(t)$ $\phi(t)$ $u(t)$'.split(' ')

pl.figure(figsize=(12, 8))
for i, yi in enumerate(y):
    pl.subplot(5,2,2*i+1)
    #pl.semilogy(yi) if i==1 else pl.plot(yi)
    pl.plot(yi)
    pl.ylabel(names[i])
    pl.grid(1)
    if i==0:
        pl.title(str(init))
pl.subplot(5,2,2)
pl.plot(np.tanh(2*y[4]))
pl.grid(1)
pl.ylabel(r'$\tanh 2u(t)$')
pl.subplot(5,2,4)
pl.semilogy(y[1]*86e6)
pl.semilogy(data['imu'], 'b')
pl.semilogy(data['imu']-data['isd'], 'b--')
pl.semilogy(data['imu']+data['isd'], 'b--')
pl.grid(1)
pl.ylabel(r'$\log I(t)$')
pl.subplot(5,2,6)
pl.plot(csv['rh'][0], 'k')
#pl.plot(y[2].mean(axis=1), 'k')
pl.plot(data['rmu'], 'b')
pl.plot(data['rmu']+data['rsd'], 'b--')
pl.plot(data['rmu']-data['rsd'], 'b--')
pl.axhline(1.0, color='r', linestyle='--')
pl.ylabel(r'$\beta(t)/\alpha$')
pl.xlabel('Days after 02/15')
pl.grid(1)
nt_yl = pl.ylim()
pl.subplot(5,2,8)
pl.plot(data['mid'], data['mpr'], 'x')
#betah,phih,uh=y[2:].mean(axis=-1)
#pl.plot(((betah-betah[0])*5+uh)/5, 'k')
pl.grid(1)
pl.xlim([-10, 220])
pl.ylabel('Google mobility reduction')
pl.subplot(5,2,10)
pl.plot(data['pid'], data['phi'], 'o')
#pl.plot(y[2]*5+y[4], 'k')
pl.ylabel('COSMO Phi')
pl.xlim([-10, 220])
#pl.ylim(nt_yl)
pl.grid(1)
pl.tight_layout()
pl.savefig('german-data-sim-bsf.png')
pl.show()

# TODO y[1,29].sum() should still hit 40k before J goes up
