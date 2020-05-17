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
    'nii': 4,
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
    'eps': 10.0,
    }
init = {
    'gamma': 0.1,
    'Gamma_': np.r_[1.0],
    'c': np.r_[0.5,5.0],
    'alpha': 0.06,
    'J': 4.0,
}
stanio.rdump('model.R', data)
stanio.rdump('init.R', init)

stanio.compile_model(cs,'model')
os.remove('fit.csv')
assert 0==os.system('./model sample algorithm=hmc engine=nuts max_depth=15 save_warmup=1 init=init.R data file=model.R output file=fit.csv')
os.system('../../cmdstan-2.22.1/bin/diagnose fit.csv')
csv = stanio.parse_csv('fit.csv')

for k, v in csv.items():
    print(k, v.shape)

from scipy import stats
pars = np.c_[csv['lp__'], csv['gamma'],csv['Gamma'], csv['c'], csv['alpha'], csv['J']][-1000:].T#, csv['mpr_a'], csv['mpr_b'], csv['phi_a'], csv['phi_b']].T
names = (r'$\log{p(\theta|y)}$ $\gamma$ $\Gamma_1$ $\Gamma_2$ $c_1$ $c_2$ $\alpha$ $J_b$').split(' ')
ranges = [(_.min(), _.max()) for _ in pars]
eps = data['eps']
priors = [(None, None), (0.1, 0.01),
          # (0.3, 0.01), (3.0, 0.1), # asymmetric case
          (1.0, 0.1), (1.0, 0.1), # symmetric case Gamma
          (0.5,0.05), (5.0,0.5), (0.06,0.006), (4.0,0.1)]
pl.figure(figsize=(7,7))
for i, ip in enumerate(pars):
    for j, jp in enumerate(pars):
        ax = pl.subplot(len(pars),len(pars),i*len(pars)+j+1)
        for k,(kip,kjp) in enumerate(zip(ip.reshape((10,-1)), jp.reshape((10,-1)))):
            k_ = k / 10.0
            c = k_, 1-k_, 0.5
            if i==j:
                pl.hist(kip, fc=c+(0.3,), density=True)
                if i>0:
                    _ = np.r_[pl.xlim()[0]:pl.xlim()[1]:100j]
                    pl.plot(_, stats.norm.pdf(_, priors[i][0], priors[i][1]*eps), 'k', alpha=0.5)
            else:
                pl.plot(kjp, kip, color=c+(0.2,), linestyle='none', marker='.')
        #ax.axis('off')
        if i!=j:
            pl.ylim(ranges[i])
            pl.grid(1)
        pl.xlim(ranges[j])
        if i==(len(pars)-1):
            pl.xlabel(names[j])
        pl.xticks(pl.xticks()[0], [])
        if j==0:
            pl.ylabel(names[i])
        else:
            pl.yticks(pl.yticks()[0], [])
pl.tight_layout()
pl.savefig('german-fit-pairs.png')
pl.show()

# now ppc
csv['pp_ih'] = np.log(stats.norm.rvs(csv['ih'],np.r_[data['isd'],np.ones((210-len(data['isd'])))*data['isd'][-1]][None,:]))
pl.figure()
for i, k in enumerate('ih rmu mpr phi'.split(' ')):
    pl.subplot(4, 1, i + 1)
    ppk = f'pp_{k}'
    q1,q2,q3=np.percentile(csv[ppk],[5,50,95],axis=0)
    pl.plot(q1, 'k--')
    pl.plot(q2, 'k')
    pl.plot(q3, 'k--')
    if i < 3:
        pl.xticks(pl.xticks()[0],[])
    pl.grid(1)
    if i == 0:
        pl.plot(data['rid'],np.log(data['imu']),'r')
        pl.plot(data['rid'],np.log(data['imu']-data['isd']),'r--')
        pl.plot(data['rid'],np.log(data['imu']+data['isd']),'r--')
    elif i==1:
        pl.plot(data['rid'],data['rmu'],'r')
        pl.plot(data['rid'],data['rmu']-data['rsd'],'r--')
        pl.plot(data['rid'],data['rmu']+data['rsd'],'r--')
    elif i==2:
        pl.plot(data['mid'],data['mpr'],'r')
    else:
        pl.plot(data['pid'],data['phi'],'r')
    pl.ylabel('log(I(t)) R0(t) Mobility CosmoPhi'.split(' ')[i])
pl.xlabel('Days after Feb 14th 2020')
pl.tight_layout()
pl.savefig('germany-ppc.png')
pl.show()
