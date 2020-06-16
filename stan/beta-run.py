import numpy as np
import pylab as pl
import stanio
import os
import scipy.stats


# cs = '~/Downloads/cmdstan-2.22.1'
# os.system('./beta sample data file=beta.R output file=beta.csv')
# os.system(f'{cs}/bin/stansummary beta.csv | grep -v ^z | grep -v ^b | grep -v ^p')
# stanio.diagnose_csvs(cs, 'beta.csv')

r = stanio.rload('beta.R')
csv = stanio.parse_csv('beta.csv')

alpha=1/14

pl.figure(figsize=(8,4))

pl.subplot(322)
b = csv['b'].mean(axis=0)
pl.plot(b, 'r--')
b0 = b[0]
p = [b[0]]
u = [0.0*b[0]]
for i, bi in enumerate(b):
	dp = bi - p[-1] + u[-1]*0.9
	du = b0 - bi - u[-1]
	p.append(p[-1] + dp/3)
	u.append(u[-1] + du/21)
u, p = [np.array(_) for _ in (u,p)]
pl.plot(u)
pl.plot(p)
pl.legend(('beta', 'phi', 'u'))

pl.subplot(321)
a,b=np.percentile(csv['b']/alpha,[5,95],axis=0)
#pl.fill_between(np.r_[:61],a,b, alpha=0.4)
pl.plot(r['r'], 'k')
pl.plot(a, 'k--')
pl.plot(b, 'k--')
pl.plot(r['iP'], r['P'])
pl.plot(np.r_[:62], (p.T-b0)*35)
pl.grid(1)


pl.subplot(312)
a,b=np.percentile(csv['z'],[5,95],axis=0)
pl.fill_between(np.r_[:60],a,b,alpha=0.4)
# pl.fill_between(np.r_[:60],a[1],b[1],alpha=0.4)
pl.title('z~N(0,1)')
pl.grid(1)

pl.subplot(349)
pl.hist(csv['a'], 30)
pl.title('a~N(0,1)')
pl.grid(1)

pl.subplot(3,4,10)
pl.hist(csv['s'], 30)
pl.title('s~lN(1,1)')
pl.grid(1)

pl.subplot(3,4,11)
pl.hist(csv['l'], 30)
pl.title('l~lN(1,1)')
pl.grid(1)

# pl.subplot(3,4,12)
# pl.hist(csv['c'], 30)
# pl.title('c~N(0,1)')
# pl.grid(1)

pl.tight_layout()
pl.show()
# pl.savefig('beta-run.png')