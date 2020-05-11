import os
import numpy as np
import matplotlib.pyplot as pl
from importlib import reload
import stanio; reload(stanio)

pz = lambda x, xs: np.abs((xs.mean(axis=0) - x[0])/(1e-16+xs.std(axis=0)))

def process_csv(csv):
    s_g = 1 - csv['gammah_z'].std(axis=0)
    s_G = 1 - csv['Gammah_z'].std(axis=0)
    z_g = pz(csv['gamma_'], csv['gammah'])
    z_G = pz(csv['Gamma_'], csv['Gammah'])
    # z-score for forecast
    z_pse = pz(csv['pse_'], csv['psep'])
    z_soI = pz(csv['soI_'], csv['soIp'])
    return s_g, s_G, z_g, z_G, z_pse[:,(90//7):].mean(), z_soI[90:].mean()

cs = '/home/user/Downloads/cmdstan-2.23.0'
stanio.compile_model(cs, 'sbc')

def do_one(i):
    csv = stanio.run('sbc', sampler_args='random seed=%d'%(1337+i,))
    zs = process_csv(csv)
    print('############   ', i, 'done!')
    return zs

# ran on cluster, results in zsa.npz

zsa = np.load('zsa.npz')
for k, v in zsa.items():
    print(k, v.shape)

globals().update(zsa)

pl.figure(figsize=(10,6))
for i, (s_, z_) in enumerate(zip(s_g.T, z_g.T)):
    pl.subplot(4,5,i+1)
    pl.plot(s_, z_, 'k.')
    pl.xlabel('Shrinkage')
    pl.ylabel('Posterior z-score')
    pl.xlim([-0.2, 1.2])
    pl.ylim([0, 2])
    pl.title('$g_%d$' % (i,))
for i, (s_, z_) in enumerate(zip(s_G.T, z_G.T)):
    pl.subplot(4,5,i+6)
    pl.plot(s_, z_, 'k.')
    pl.xlabel('Shrinkage')
    pl.ylabel('Posterior z-score')
    pl.xlim([-0.2, 1.2])
    pl.ylim([0, 2])
    pl.title('$\Gamma_%d$' % (i,))
for i, (s_, z_) in enumerate(zip(s_ic.T, z_ic.T)):
    pl.subplot(4,5,i+11)
    pl.plot(s_, z_, 'k.')
    pl.xlabel('Shrinkage')
    pl.ylabel('Posterior z-score')
    pl.xlim([-0.2, 1.2])
    pl.ylim([0, 2])
    pl.title(['$I(t=0)$', r'$\beta(t=0)$', r'$\phi(t=0)$', '$u(t=0)$'][i])
#pl.subplot(3,3,8); pl.hist(z_pse); pl.title('Posterior predictive z-score PSE')
pl.subplot(4,5,16); pl.hist(z_soI); pl.title('Posterior predictive z-score sum(I)'); 
pl.tight_layout()
pl.savefig('sbc-fixed-sensitivity.png')
pl.show()
