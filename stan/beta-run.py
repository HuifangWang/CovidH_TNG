import numpy as np
import pylab as pl
import stanio
import os
import scipy.stats


# the Stan run can sometimes get poorly initialized and
# results are crap, just rerun it.
cs = "~/crypt/cmdstan-2.22.1"
os.system("./beta sample data file=beta.R output file=beta.csv")
stanio.diagnose_csvs(cs, "beta.csv")

r = stanio.rload("beta.R")
csv = stanio.parse_csv("beta.csv")

alpha = 1 / 14


def ts_ci(v, ci=95, axis=0, it=None):
    p = (100 - ci) / 2
    a, b = np.percentile(v, [p, 100 - p], axis=axis)
    it = np.r_[: a.size] if it is None else it
    pl.fill_between(it, a, b, alpha=0.5)


pl.figure(figsize=(5, 7))

pl.subplot(311)
pl.plot(r["r"])
ts_ci(csv["gq_r"])
pl.plot(r["iP"], r["P"])
ts_ci(csv["gq_P"], it=r["iP"])
pl.legend(("$R_0(t)$", r"$\phi(t)$", r"$\hat{R_0}(t)$", r"$\hat{\phi}(t)$"))
pl.grid(1)

pl.subplot(323)
pl.hist(csv["a"], 30)
pl.title("a~N(0,1)")
pl.grid(1)

pl.subplot(324)
pl.hist(csv["s"], 30)
pl.title("s~lN(0,1)")
pl.grid(1)

pl.subplot(325)
pl.hist(csv["l"], 30)
pl.title("l~lN(1,1)")
pl.grid(1)

pl.subplot(326)
x = np.r_[:0.5:0.01]
y = scipy.stats.lognorm.pdf(x, 1, loc=0)
pl.plot(x, y, "k")
pl.hist(csv["g"], 50, density=True)
pl.legend((r"$p(\gamma)=\log N(0.1,1)$", r"$\hat{\gamma}$",))
# pl.title('g~lN(0.1,1)')
pl.xlabel(r"$\gamma$")
pl.title(r"$p(\gamma\|R0(t),\phi(t))$")
pl.yticks(pl.yticks()[0], [])
pl.grid(1)


pl.tight_layout()
pl.show()
# pl.savefig('beta-run.png')
