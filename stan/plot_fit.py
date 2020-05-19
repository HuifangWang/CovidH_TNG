import numpy as np
import pylab as pl
import pandas as pd
import os
import datetime


## load all the data


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
pl.title("Data for Germany")
pl.xlabel("Days after Feb 15th")
pl.legend("$\mu_{R0(t)}$ $\sigma^2_{R0(t)}$ $\log{\mu_{I(t)}}$ $\log{\sigma^2_{I(t)}}$ $Mobility$ $Cosmo\phi$".split(' '))
pl.grid(1)
pl.savefig('german-data.png', dpi=300)
pl.show()


import stanio
cs = '/Volumes/crypt/cmdstan-2.22.1'

stanio.compile_model(cs, 'model')

data = {
    'nii': 4,
    'nr': len(de_idx),
    'rid': de_idx + 1,
    'rmu': de_rmu,
    'rsd': de_rsd,
    # only loosely constrain I with ICL estimates
    'ni': len(de_idx),
    'iid': de_idx + 1,
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
    # model comparision params
    'mc_use_pse': 1,
    'mc_use_groups': 1,
    'mc_Gamma_cov': 3.0,
    'mc_ind_Gamma': 0,
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


# model comparison cases
defaults = {
    'mc_use_pse': 1,
    'mc_use_groups': 1,
    'mc_Gamma_cov': 3.0,
    'mc_ind_Gamma': 0,
}
# values deviate from defaults
cases = {
    'no_pse': {'mc_use_pse': 0, 'mc_use_groups': 0},
    'no_group': {'mc_use_groups': 0},
    'eq_Gamma': {'mc_Gamma_cov': 1.0},
    'cov_Gamma': {'mc_Gamma_cov': 3.0},
    'ind_Gamma': {'mc_ind_Gamma': 1}
}
# concrete cases hold all values
ccases = {key: dict(**defaults) for key in cases.keys()}
for key, val in cases.items():
    print(key, val)
    ccases[key].update(val)

# run all variants of the model
fits = {}
for key, val in ccases.items():
    data.update(val)
    stanio.rdump(f'data-{key}.R', data)
    assert 0==os.system(f'./model random seed=1337 sample save_warmup=1 algorithm=hmc engine=nuts max_depth=15 init=init.R data file=data-{key}.R output file=fit.csv'), f'{key} failed!'
    os.system(f'cp fit.csv fit-{key}.csv')
    os.system(f'../../cmdstan-2.22.1/bin/diagnose fit-{key}.csv')
    fits[key] = stanio.parse_csv(f'fit-{key}.csv')


### model comparison Bayes factors, just printing to command line and included in paper
from scipy import stats
def log_lik_IR(fit):
    rh = fit['rh'][-1000:, data['rid']]
    ih = fit['ih'][-1000:, data['iid']]
    llr = np.log(stats.norm.pdf(data['rmu'], rh, data['rsd']))
    lli = np.log(stats.norm.pdf(data['imu'], rh, data['isd']))
    return np.c_[llr, lli]
def loo_ir(fit):
    from psis import psisloo
    loo, loos, ks = psisloo(log_lik_IR(fit))
    return loo
loo_no_pse = loo_ir(fits['no_pse'])
for key, val in fits.items():
    if key == 'no_pse':
        continue
    l=loo_ir(val)
    print(key, np.exp(l - loo_no_pse))

# pairs plots (not in paper)
csv = fits['cov_Gamma']
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


# Figure X2: posteriors & posterior predictive
fig = pl.figure()
csv_ = {k:v[-1000:] for k,v in fits['cov_Gamma'].items()}
# posterior distributions 3x2+4x1
def php(x,mu,sd):
    "plot a posterior histogram and prior"
    pl.hist(x, 40, alpha=0.5, fc='red', density=True, orientation="horizontal")
    lo = x.min()# min(x.min(), mu-sd)
    hi = x.max() #max(x.max(), mu+sd)
    _ = np.r_[lo:hi:100j]
    pl.plot(stats.norm.pdf(_, mu, sd), _, 'k', alpha=0.5, linewidth=5)
    pl.yticks(np.percentile(x,[2.5,50,97.5]))
    # try to have reasonable number of significant digits:
    from matplotlib.ticker import FormatStrFormatter
    pl.gca().yaxis.set_major_formatter(FormatStrFormatter(f'%.{max(2,int(2-np.log10(x.mean())))}f'))
    # x axis is irrelevant here
    pl.xticks([])
    pl.grid(1)
fig.set_constrained_layout(True)
gs = fig.add_gridspec(4,4)
ax = fig.add_subplot(gs[0,0]); php(csv_['gamma'],0.1,0.1); pl.title(r'$\gamma$')
ax = fig.add_subplot(gs[0,1]); php(csv_['alpha'],0.1,0.1); pl.title(r'$\alpha$')
ax = fig.add_subplot(gs[1,0]); php(csv_['Gamma'][:,0],1/3,1/3); pl.title(r'$\Gamma_1$')
ax = fig.add_subplot(gs[1,1]); php(csv_['Gamma'][:,1],1*3,1*3); pl.title(r'$\Gamma_2$')
ax = fig.add_subplot(gs[2,0]); php(csv_['c'][:,0],0.5,0.5); pl.title(r'$c_1$')
ax = fig.add_subplot(gs[2,1]); php(csv_['c'][:,1],5,5); pl.title(r'$c_2$')
ax = fig.add_subplot(gs[3,0]); php(csv_['mpr_a'],1,1); pl.title(r'$m_{mob}$')
ax = fig.add_subplot(gs[3,1]); php(csv_['phi_a'],1,1); pl.title(r'$m_{\phi}$')
# PPC time series
csv_['pp_ih'] = np.log(stats.norm.rvs(csv['ih'],np.r_[data['isd'],np.ones((210-len(data['isd'])))*data['isd'][-1]][None,:]))
for i, k in enumerate('ih rmu mpr phi'.split(' ')):
    #pl.subplot(4, 2, 2*i + 2)
    ax = fig.add_subplot(gs[i,2:])
    ppk = f'pp_{k}'
    #q=np.percentile(csv[ppk],[2.5,10,25,75,90,97.5],axis=0)
    if 0:
        pl.plot(q1, 'k--')
        pl.plot(q2, 'k')
        pl.plot(q3, 'k--')
    else:
        for a,p1,p2 in [(0.1,2.5,97.5), (0.3,25,75)]:
            q1, q2 = np.percentile(csv_[ppk], [p1,p2], axis=0)
            pl.fill_between(np.r_[:len(q1)], q1, q2, color='k', alpha=a)
        # pl.plot(q2, 'k')
    if i < 3:
        pl.xticks(pl.xticks()[0],[])
    pl.grid(1)
    if i == 0:
        pl.plot(data['rid'],np.log(data['imu']),'r', alpha=0.5, linewidth=2)
        #pl.plot(data['rid'],np.log(data['imu']-data['isd']),'r--')
        #pl.plot(data['rid'],np.log(data['imu']+data['isd']),'r--')
        #pl.title('Predictive distributions for Germany')
    elif i==1:
        pl.plot(data['rid'],data['rmu'],'r', linewidth=2, alpha=0.5)
        #pl.plot(data['rid'],data['rmu']-data['rsd'],'r--')
        #pl.plot(data['rid'],data['rmu']+data['rsd'],'r--')
        pl.yticks([1,4])
    elif i==2:
        pl.plot(data['mid'],data['mpr'],'r', alpha=0.5, linewidth=2)
    else:
        pl.plot(data['pid'],data['phi'],'r', alpha=0.5, linewidth=2)
    pl.ylabel('log(I(t)) R0(t) Mobility Cosmo$\phi$'.split(' ')[i])
pl.xlabel('Days after Feb 14th 2020')
pl.tight_layout()
pl.savefig('germany-inversion.png', dpi=300)
pl.show()



# TODO figure showing forecast CI for different interventions?



"""
## another nice figure to show the different among models in terms of the data
# normalize each I(t) and R(t) estimate to the data as a z score
# plot CI bands
# expect to see tightest bands for
pl.figure()
leg = []
for i,k in enumerate('no_pse cov_Gamma'.split()):
    v = fits[k]
    l,m,h=np.percentile(v['rh'][-1000:],[0.5,50,99.5],axis=0)
    pl.fill_between(np.r_[:l.size], l,h, color='bk'[i], alpha=0.2)
    # pl.plot(m, 'bk'[i], alpha=0.5)
    leg.append(k)
pl.plot(data['rid'], data['rmu'], 'r', alpha=0.4)
pl.legend('no_pse cov_Gamma R0'.split())
pl.show()
"""

