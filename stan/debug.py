stanio.compile_model(cs, 'model')
init = {
    'gamma': 0.1,
    'Gamma_': np.r_[1.0, 1.0],
    'c': np.r_[0.5,5.0],
    'alpha': 0.1,
    'J': 0.7,
    'ur': -1
}
stanio.rdump('init.R', init)
try:
    os.unlink('fit.csv')
except:
    pass
assert 0==os.system(f'./model sample adapt delta=0.95 init=init.R data file=model.R output refresh=1 file=fit.csv'), f'failed!'
stanio.diagnose_csvs(cs, 'fit.csv')

csv_ = stanio.parse_csv(f'fit.csv')
csv = {k:v[-1000:] for k,v in csv_.items()}


# imu = csv['pp_imu']
# prc = csv['pp_cases']
# pl.figure()
# pl.semilogy(de_idx, de_irc)
# pl.semilogy(de_idx, de_imu)
# pl.semilogy(np.percentile(imu,2.5,axis=0), 'k--')
# pl.semilogy(np.percentile(imu,97.5,axis=0), 'k--')
# pl.semilogy(np.percentile(imu,2.5,axis=0), 'k--')
# pl.semilogy(np.percentile(imu,97.5,axis=0), 'k--')
# pl.legend('icl_reports icl_predicted pp2.5% pp97.5%'.split(' '))
# pl.title('Daily new infections')
# pl.savefig('debug.png',dpi=100)
# pl.close('all')

# pl.figure()
# pl.plot(csv['gamma'], csv['alpha'], 'k.')
# pl.show()
    

# import datetime as dt
# ts = [
#     dt.datetime(year=2020, month=2, day=15),
#     dt.datetime(year=2020, month=3, day=6),
#     dt.datetime(year=2020, month=3, day=12),
#     dt.datetime(year=2020, month=3, day=14),
#     dt.datetime(year=2020, month=3, day=22),
#     dt.datetime(year=2020, month=5, day=6)
#     ]
# t0 = ts[0]
# for t in ts:
#     print(t, 'is day', (t - t0).days)
