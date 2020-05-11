// SBC according to Betancourt's workflow:
// 1. generate synthetic datasets
// 2. generate one fit per dataset
// 3. collect posterior z-score and shrinkage for each fit
// 4. plot z-score vs shrinkage, ideally lower right corner

#include "functions.stan"

transformed data {
  real noise = 0.1;
  int nt = 210;
  int nt_7 = nt / 7;
  int no = 90;
  int no_7 = no / 7;
  real ur = 2.0;
  real ur_sd = 0.3;
  real pseh_sd = 0.3;
  vector[4] ic;
  vector[5] gamma;
  vector[2] Gamma;
  matrix[2,nt-1] dwI;
  vector[5] arI;
  matrix[10,nt] yt;
  matrix[3,nt_7] pse;
  vector[nt] soI;
  {
    int have_nan = 1;
    int draws = 0;
    while ((have_nan==1) && (draws < 100)) {
      draws += 1;
      gamma = gamma_xfm(normal_vec_rng(5));
      Gamma = Gamma_xfm(normal_vec_rng(2));
      ic = ic_xfm(normal_vec_rng(4));
      dwI = mat_rng(2, nt - 1);
      arI = normal_vec_rng(5)/10;
      yt = run(gamma, Gamma, noise, dwI, arI, ic);
      pse = sub_pse(yt);
      soI = obs_I(yt, ur);
      if (count_nan_mat(yt) == 0)
	have_nan = 0;
    }
    if (draws == 100)
      reject("unable to generate stable simulation with 100 draws, giving up.");
    print("stable simulation generated with ", draws, "draws.");
  }
}

parameters {
  vector[5] gammah_z;
  vector[2] Gammah_z;
  vector[4] ic_z;
  matrix[2, no - 1] dwIh;
  vector[5] arIh;
}

transformed parameters {
  vector[5] gammah = gamma_xfm(gammah_z);
  vector[2] Gammah = Gamma_xfm(Gammah_z);
  vector[4] ich = ic_xfm(ic_z);
  matrix[10,no] yth = run(gammah, Gammah, noise, dwIh, arIh, ich);
  matrix[3,no_7] pseh = sub_pse(yth);
  row_vector[no] soIh = exp(yth[3,]) + exp(yth[4,]);
  vector[no] urh = soIh' ./ soI[1:no];
  reject_nan_mat([urh']);
}

model {
  gammah_z ~ std_normal();
  Gammah_z ~ std_normal();
  ic_z ~ std_normal();
  to_vector(dwIh) ~ std_normal();
  to_vector(arIh) ~ normal(0, 0.1);
  to_vector(pse[,1:cols(pseh)]) ~ normal(to_vector(pseh), pseh_sd);
  target += normal_lpdf(to_vector(urh) | ur, ur_sd);
}

generated quantities {
  // recover true values
  vector[5] gamma_ = gamma;
  vector[2] Gamma_ = Gamma;
  vector[4] ic_ = ic;
  matrix[3,nt_7] pse_ = pse;
  vector[nt] soI_ = soI;
  // log_lik for waic/psis
  vector[3*(no_7)+no] log_lik;
  {
    for (j in 1:no_7)
      for (i in 1:3)
	log_lik[(j-1)*3+i] = normal_lpdf(pse[i,j]|pseh[i,j],pseh_sd);
    for (i in 1:no)
      log_lik[3*no_7-1+i] = normal_lpdf(urh[i]|ur,ur_sd);
  }
  // ppc
  matrix[3,nt_7] psep;
  vector[nt] soIp;
  {
    matrix[10,nt] ytp;
    int have_nan = 1;
    while (have_nan==1) {
      ytp = run(gammah, Gammah, noise,
	        append_col(dwIh, mat_rng(2, nt - no)),
	        arIh, ich);
      psep = sub_pse(ytp);
      for (i in 1:3)
	for (t in 1:nt_7)
	  psep[i,t] += normal_rng(0,pseh_sd);
      soIp = obs_I(yt, 1);
      for (t in 1:nt)
	soIp[t] /= normal_rng(ur, ur_sd);
      if (count_nan_mat(yt) == 0)
	have_nan = 0;
    }
  }
}
