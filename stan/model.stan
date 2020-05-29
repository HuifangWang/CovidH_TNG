functions {
  row_vector rvdiff(row_vector x) {
    return x[2:cols(x)] - x[1:cols(x)-1];
  }
  
  matrix run(int nt, real gamma, vector c, vector Gamma, vector ic, real alpha, real J, int nii, real[] Jb, int reject_) {
    real dt = 1.0;
    real sqrt_dt = sqrt(dt);
    real Lambda = 1.0 / 80.0 / 365.0 / 2.0;
    real mu = 1.0 / 80.0/ 365.0;
    // real alpha = 0.2;
    real phi0 = 1.0;
    real sigma = 30.0;
    real u0 = 0.0;
    real beta0 = ic[2];
    real h = 0.2;
    real tau = 1.0;
    // initial conditions
    real I0 = ic[1] / 2;
    real S0 = (1.0/2) - I0;
    matrix [2,5] y;
    for (g in 1:2) y[g,] = append_col([S0,I0], ic[2:4]');
    matrix[2, 5] k[4]; // rk4
    matrix[2, 5] y_; // rk4
    matrix[10, nt] yt;
    yt[,1] = to_vector(y);
    for (t in 2:nt) {
      if (Jb[t] < Jb[t-1])
	y[,5] = atanh(tanh(y[,5]) * square(Jb[t]/Jb[t-1]));
      real Jb_ = J * Jb[t];
      real dt_ = dt / nii;
      for (ii in 1:nii) { // dt=0.25 to handle bigger J jumps
	for (i in 1:4) {
	  /* intermediate state for RK4 i'th step */
	  if (i == 1)
	    y_ = y;
	  else if (i == 2 || i == 3)
	    y_ = y + dt_/2*k[i-1];
	  else // i == 4
	    y_ = y + dt_*k[i-1];
	  /* ODE RHS each group */
	  real sum_I = sum(y_[,2]);
	  real sum_S = sum(y_[,1]);
	  for (g in 1:2) {
	    real S	= y_[g,1];
	    real I	= y_[g,2];
	    real beta	= y_[g,3];
	    real phi	= y_[g,4];
	    real u	= y_[g,5];
        
	    real dS	= Lambda - beta*sum_I*S - mu*S;
	    // d log I = beta*(sum_S-alpha) - mu
	    real dI	= beta*(sum_S-alpha)*I - mu*I;
	    real dbeta	= -(beta - beta0 + Jb_*beta) + gamma*(phi-phi0)*(beta - beta0);
	    real dphi	= -(phi - phi0 + h*Jb_)/tau - c[g]*tanh(2*u)*(Jb_>0);
	    real du       = -Gamma[g]*(beta - beta0)/sigma; // u0 zero

        if ((beta + dt_*dbeta)<0) {
          reject("beta<0: dbeta=", dbeta, " = ", "-(beta - beta0 + Jb_*beta)=", -(beta - beta0 + Jb_*beta), "  - gamma*(phi-phi0)*(beta - beta0))=",- gamma*(phi-phi0)*(beta - beta0), "\n\tgamma = ", gamma, ", phi - phi0 = ", phi, " - ", phi0);
        }

	    k[i][g,] = [dS, dI, dbeta, dphi, du];
	  }
	}
	{
	  matrix[2,5] yu = (dt_/6)*(k[1] + 2*k[2] + 2*k[3] + k[4]); // rk4 update
      /*
	  if ((y[1,3]+yu[1,3])<0 || (y[2,3]+yu[2,3])<0) {
	    print(" !! beta about to go negative...");
	    print("t = ", t, "; ii = ", ii);
	    print("beta - beta0 = ", y_[,3]' - beta0);
	    print("Jb_ * beta", Jb_ * y_[,3]');
	    print("phi - phi0", y_[,4]' - phi0);
	    print("beta = ", y_[,3]', " += dbeta = ", k[4][,3]');
	    reject(" !! FATALITYYYY !!");
	  }
      */
	  y += yu;
	}
      }
      for (g in 1:2)
	for (j in 1:3)
	  if (y[g,j] < 0) {
	    if (reject_)
	      reject("t=",t," y[",g,",",j,"] < 0, but should not be!");
	    y[g,j] = 0.0;
	  }
      yt[,t] = to_vector(y); // save
    }
    return yt;
  }

  matrix yt2pse(matrix yt) {
    int nt = cols(yt);
    matrix[3,nt] pse;
    for (i in 1:3) {
      int j = (i - 1)*2 + 5;
      pse[i] = (yt[j,] + yt[j+1,]) * 0.5;
    }
    // pse[2] should be phi_dot, smoothed weekly
    for (t_ in 0:(nt-1)) {
      int t = nt - t_;
      pse[2,t] = mean(pse[2,max({1,t-6}):t]);
    }
    pse[2] = append_col([0], rvdiff(pse[2]));
    return pse;
  }
}

data {
  int nr;
  int rid[nr];
  real rmu[nr]; // r(t) mean
  real rsd[nr]; // r(t) ci/4
  int ni;
  int iid[ni];
  real imu[ni]; // i(t) mean
  real isd[ni]; // i(t) ci/4
  int nc;
  int cid[nc];
  int cases[nc];
  int nm;
  int mid[nm];
  real mpr[nm]; // mobility proporation reduction
  int np;
  int pid[np];
  real phi[np]; // phi values
  real I0;
  int nii; // steps per day to take
  real eps; // loosen priors
  // model comparison flags
  int mc_use_pse;
  int mc_use_groups;
  real mc_Gamma_cov;
  int mc_ind_Gamma;
}

transformed data {
  int nt = 210;
  // from Spase for Germany
  real Jv[7] =	{0.0, 0.05, 0.1, 0.15,  3,   2,   1};
  int Ji[7] =	{  1,   21,  27,   29, 37,  82, 82+60};
  real Jb[nt] = rep_array(-1.0, nt);
  Jb[Ji] = Jv;
  for (t in 1:nt)
    if (Jb[t] < 0)
      Jb[t] = Jb[t - 1];
  for (t in 80:86)
    print("Jb[", t, "] = ", Jb[t]);
}

parameters {
  //real<lower=0> gamma_;
  //vector<lower=0>[2] Gamma_;
  vector<lower=0>[2] c_;
  real<lower=0> alpha;
  //real<lower=0> J;
  real mpr_a;
  real mpr_b;
  real phi_a;
  real phi_b;
  // real<lower=0> beta0;
  real ur;
}


transformed parameters {
  vector[4] ic = [I0, 0.4, 1, 0]'; // I0, beta, phi u
  real gamma_ = 0.1;
  real gamma = mc_use_pse ? gamma_ : 0.0;
  vector[2] Gamma;
  vector[2] c;
  real J = 1.5;
  vector[2] Gamma_ = [0.3,3.0]';
  if (mc_use_groups) {
    c = c_;
    Gamma = mc_ind_Gamma ? Gamma_ : [Gamma_[1]/mc_Gamma_cov, Gamma_[1]*mc_Gamma_cov]';
  } else {
    Gamma = [Gamma_[1], Gamma_[1]]';
    c = [c_[1]*5, c_[1]*0.2]';
  }
  matrix[10, nt] yt = run(nt, gamma*mc_use_pse, c, Gamma, ic, alpha, J, nii, Jb, 1);
  row_vector[nt] rh = (yt[5,]+yt[6,]) / alpha * 0.5;
  row_vector[nt] It = (1.0-(yt[1,]+yt[2,]))*86e6;
  row_vector/*<lower=0>*/[nt-1] ih = It[2:nt]-It[1:nt-1];
  matrix[3,nt] pse = yt2pse(yt);
  real eur = exp(ur);
}

model {
  // SIR
  // J ~ normal(1.0,0.0001*eps) T[0,];
  alpha ~ normal(0.1,0.01*eps) T[0,];
  // beta0 ~ normal(0.4,0.04*eps) T[0,];
  // gamma ~ normal(0.1, 0.01*eps) T[0,];
  // Gamma_[1] ~ normal(1.0,0.1*eps) T[0,];
  // Gamma_[2] ~ normal(1.0,0.1*eps) T[0,];
  c_[1] ~ normal(0.5,0.05*eps) T[0,];
  c_[2] ~ normal(5.0,0.5*eps) T[0,];
  // rmu ~ normal(rh[rid],rsd);
  imu ~ normal(ih[iid],isd);
  ur ~ std_normal();
  for (i in cid)
    cases[i] ~ poisson(ih[i]*eur);
  // PSE
  mpr_a ~ normal(1, 1);
  mpr_b ~ normal(-1, 1);
  phi_a ~ normal(1, 1);
  phi_b ~ normal(-1,1);
  to_vector(mpr) ~ normal(pse[1,mid]*mpr_a+mpr_b, 0.1);
  to_vector(phi) ~ normal(pse[2,pid]*phi_a+phi_b, 0.3);
}

generated quantities {
  matrix[10, nt] yt_ = yt;
  vector[nt] pp_imu;
  vector[nt] pp_mpr;
  vector[nt] pp_phi;
  int pp_cases[nt];
  for (t in 1:nt) {
    int t_ = min({t,nt-1});
    pp_imu[t] = normal_rng(ih[t_],ih[t_]/4);
    pp_mpr[t] = normal_rng(mpr_a*pse[1,t]+mpr_b, 0.1);
    pp_phi[t] = normal_rng(phi_a*pse[2,t]+phi_b, 0.3);
    pp_cases[t] = poisson_rng(ih[t_]*eur);
  }
}

