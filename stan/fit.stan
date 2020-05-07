#include "functions.stan"

data {
  real noise;
  int nt; // total days
  int no;
  vector[nt] soI;
  int use_pse;
  int use_soI;
  matrix[3,nt/7] pse;
  real pseh_sd;
  real ur;
  real ur_sd;
  real gamma_mu;
  real gamma_sd;
  vector[2] Gamma_mu;
  vector[2] Gamma_sd;
}

transformed data {
  real noiseh = noise; // no estimate noise for now
  int arp = 5;
}

parameters {
  real<lower=0> gammah;
  vector<lower=0>[2] Gammah;
  matrix[2, no - 1] dwIh;
  vector[arp] arIh;
}

transformed parameters {
  matrix[10, no] yth = run(gammah, Gammah, noiseh, dwIh, arIh);
  matrix[3,no/7] pseh = sub_pse(yth);
  row_vector[no] soIh = exp(yth[3,]) + exp(yth[4,]);
  vector[no] urh = soIh' ./ soI[1:no];
}

model {
  gammah ~ normal(gamma_mu, gamma_sd) T[0,];
  Gammah[1] ~ normal(Gamma_mu[1], Gamma_sd[1]) T[0,];
  Gammah[2] ~ normal(Gamma_mu[2], Gamma_sd[2]) T[0,];
  to_vector(dwIh) ~ std_normal();
  to_vector(arIh) ~ normal(0, 0.1);
  if (use_soI)
    target += normal_lpdf(to_vector(urh) | ur, ur_sd);
  if (use_pse)
    to_vector(pse[,1:cols(pseh)]) ~ normal(to_vector(pseh), pseh_sd);
}

generated quantities {
  /* matrix[10, nt] ytp = run(gammah, Gammah, noiseh, append_col(dwIh, mat_rng(2,nt-no)), arIh); */
  /* matrix[2, nt] Ip = ytp[3:4,] + lur; */
  /* matrix[3, nt/7] psep = sub_pse(ytp); */
  /* vector[no] log_lik; */
  /* for (t in 1:no) */
  /*   log_lik[t] = normal_lpdf(lurh[t] | lur, lur_sd); */
}
