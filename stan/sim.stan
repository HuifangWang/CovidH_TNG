// generate synthetic data

#include "functions.stan"

data {
  real noise;
  int nt;
  real ur;
  real gamma;
  vector[2] Gamma;
}

generated quantities {
  matrix[2,nt-1] dwI = mat_rng(2, nt - 1);
  vector[5] arI = mat_rng(5,1)[,1] / 10;
  matrix[10,nt] yt = run(gamma, Gamma, noise, dwI, arI);
  matrix[3,nt/7] pse = sub_pse(yt);
  vector[nt] soI = obs_I(yt, ur);
}
