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

# do an optimization first
os.system('./model optimize algorithm=newton data file=model.R output file=opt.csv refresh=1')
csv = stanio.parse_csv('opt.csv')
init2 = {}
for key in init.keys():
    init2[key] = csv[key][0]
for key in 'mpr_a mpr_b phi_a phi_b'.split():
    init2[key] = csv[key][0]
print(init2)
stanio.rdump('init2.R', init2)

os.system('./model sample init=init2.R data file=model.R output file=fit.csv refresh=1')
os.system('../../cmdstan-2.22.1/bin/diagnose fit.csv')
csv = stanio.parse_csv('fit.csv')

for k, v in csv.items():
    print(k, v.shape)

pars = np.c_[csv['gamma'],csv['Gamma'], csv['c'], csv['alpha'], csv['mpr_a'], csv['mpr_b'], csv['phi_a'], csv['phi_b']].T
pl.figure(figsize=(10,10))
for i, ip in enumerate(pars):
    for j, jp in enumerate(pars):
        ax = pl.subplot(len(pars),len(pars),i*len(pars)+j+1)
        for k,(kip,kjp) in enumerate(zip(ip.reshape((10,-1)), jp.reshape((10,-1)))):
            k_ = k / 10.0
            c = (k_,1-k_,0.5,0.4)
            if i==j:
                pl.hist(kip, fc=c)
            else:
                pl.plot(jp, ip, color=c, symbol=',')
        ax.axis('off')
pl.tight_layout()
pl.savefig('german-fit-pairs.png')
pl.show()

print(pars[:,:10].mean(axis=1))
