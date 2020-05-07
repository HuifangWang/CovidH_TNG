# simulation based calibration

import os
import numpy as np
import matplotlib.pyplot as pl

import stanio

def zs(draws, true, prior_sd):
    "Betancourt's posterior z-score & shrinkage."
    z = np.abs((draws.mean() - true)/draws.std())
    s = 1 - draws.std()/prior_sd
    return z, s

cs = "/home/user/Downloads/cmdstan-2.23.0"

def run_sim(noise=0.1, nt=210, ur=2.0, gamma=0.1, Gamma=np.r_[0.1, 0.5]):
    os.system(f"make -C {cs} $PWD/sim")
    data = {
        'noise': noise,
        'nt': nt,
        'ur': ur,
        'gamma': gamma,
        'Gamma': Gamma,
    }
    stanio.rdump('sim.R', data)
    os.system("./sim sample num_warmup=0 num_samples=1 algorithm=fixed_param "
              "data file=sim.R output file=sim.csv")
    sim = {k:v[0] for k,v in stanio.parse_csv('sim.csv').items()}
    pse = sim['pse']
    soI = sim['soI']
    return data, pse, soI


def fit_sim(params):
    sim_data, sim_pse, sim_soI = run_sim()
    fit_data = {
        'noise': sim_data['noise'],
        'nt': sim_data['nt'],
        'no': 90,
        'soI': sim_soI,
        'use_pse': 1,
        'use_soI': 1,
        'pse': sim_pse.T,
        'pseh_sd': 0.3,
        'ur': 2.0,
        'ur_sd': 1.0,
        'gamma_mu': 0.1,
        'gamma_sd': 1.0,
        'Gamma_mu': np.r_[0.1, 0.5],
        'Gamma_sd': np.r_[1.0, 1.0],
    }
    fit_data.update(params)
    stanio.compile_model('fit', cs)
    stanio.rdump('fit.R', fit_data)
    os.system('./fit sample data file=fit.R output file=fit.csv')
    fit = stanio.parse_csv('fit.csv')
    return fit

fit_wo_soi = [fit_sim({'use_soI': 0, 'use_pse': 0}) for _ in range(100)]
fit_wo_pse = [fit_sim({'use_pse': 0}) for _ in range(100)]
fit_w_pse = [fit_sim({'use_pse': 1}) for _ in range(100)]
fit_only_pse = [fit_sim({'use_pse': 1, 'use_soI': 0}) for _ in range(100)]

# calibration with z, s, rhat, psis
fitzs = {}
for key in 'fit_wo_soi fit_wo_pse fit_w_pse fit_only_pse'.split():
    fits = eval(key)
    fitzs[key] = []
    for fit in fits:
        fitzs[key].append([
            zs(fit['gammah'], 0.1, 1.0),
            zs(fit['Gammah'][:,0], [0.1, 0.5][0], 1.0),
            zs(fit['Gammah'][:,1], [0.1, 0.5][1], 1.0),
        ])
    fitzs[key] = np.array(fitzs[key])


# TODO add diagnose checks & PSIS LOO checks & PPC checks


pl.figure()
fit_names = 'fit_wo_soi fit_wo_pse fit_w_pse fit_only_pse'.split()
for i,l in enumerate(r'$\hat{\gamma}$ $\hat{\Gamma}_1$ $\hat{\Gamma}_2$'.split()):
    pl.subplot(2, 2, i+2)
    for j, k in enumerate(fit_names):
        c = 'kbgr'[j]
        z,s = fitzs[k][:, i].T
        pl.plot(s, z, '.'+c)
    pl.ylim([0, 2])
    pl.xlim([0, 1])
    pl.xlabel('Shrinkage')
    pl.ylabel('Posterior z-score')
    pl.grid(1)
    pl.title(l)
pl.legend(['fit_no_data' if _=='fit_wo_soi' else _ for _ in fit_names])
pl.tight_layout()
pl.savefig('sbc-data-sensitivity.png')
pl.show()
