functions {
  row_vector rvdiff(row_vector x) {
    return x[2:cols(x)] - x[1:cols(x)-1];
  }
  
  matrix run(int nt, real gamma, vector c, vector Gamma, vector ic, real alpha, real J, int nii) {
    real dt = 1.0;
    real sqrt_dt = sqrt(dt);
    real Jb = 0.0;
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
      real Jb_ = Jb;
      if (t == 33) Jb = J; // 34 for Germany
      if (t == 90) Jb = J/2;
      if (t == 150) Jb = 0.0;
      // sigmoidal reset
      if (Jb < Jb_)
	y[,5] = atanh(tanh(y[,5]) * square(Jb/Jb_));
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
	    real dbeta	= -(beta - beta0 + Jb*beta) - gamma*(phi-phi0)*(beta - beta0);
	    real dphi	= -(phi - phi0 + h*Jb)/tau - c[g]*tanh(2*u)*(Jb>0);
	    real du       = -Gamma[g]*(beta - beta0)/sigma; // u0 zero

	    k[i][g,] = [dS, dI, dbeta, dphi, du];
	  }
	}
	y += (dt_/6)*(k[1] + 2*k[2] + 2*k[3] + k[4]); // rk4 update
      }
      for (g in 1:2)
	for (j in 1:3)
	  if (y[g,j] < 0) {
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
  real imu[nr]; // i(t) mean
  real isd[nr]; // i(t) ci/4
  int nm;
  int mid[nm];
  real mpr[nm]; // mobility proporation reduction
  int np;
  int pid[np];
  real phi[np]; // phi values
  real I0;
  int nii; // steps per day to take
  real eps; // loosen priors
}

transformed data {
  int nt = 210;
}

parameters {
  real<lower=0> gamma;
  real<lower=0> Gamma_;
  vector<lower=0>[2] c;
  real alpha;
  real J;
  real mpr_a;
  real mpr_b;
  real phi_a;
  real phi_b;
}

transformed parameters {
  vector[4] ic = [I0, rmu[1]*alpha, 1, 0]'; // I0, beta, phi u
  vector[2] Gamma = [Gamma_, Gamma_]'; // [Gamma_/3, Gamma_*3]';
  matrix[10, nt] yt = run(nt, gamma, c, Gamma, ic, alpha, J, nii);
  row_vector[nt] rh = (yt[5,]+yt[6,]) / alpha * 0.5;
  row_vector[nt] ih = (yt[3,]+yt[4,]) * 86e6;
  matrix[3,nt] pse = yt2pse(yt);
}

model {
  // SIR
  J ~ normal(4.0,0.1*eps) T[0,];
  alpha ~ normal(0.1,0.01*eps) T[0,];
  gamma ~ normal(0.1, 0.01*eps) T[0,];
  Gamma_ ~ normal(1.0,0.1*eps) T[0,];
  //Gamma[2] ~ normal(3.0,0.1*eps) T[0,];
  c[1] ~ normal(0.5,0.05*eps) T[0,];
  c[2] ~ normal(5.0,0.5*eps) T[0,];
  rmu ~ normal(rh[rid],rsd);
  imu ~ normal(ih[rid],isd);
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
  vector[nt] pp_rmu;
  vector[nt] pp_mpr;
  vector[nt] pp_phi;
  for (t in 1:nt) {
    pp_rmu[t] = normal_rng(rh[t],rsd[min({nr,t})]);
    pp_mpr[t] = normal_rng(mpr_a*pse[1,t]+mpr_b, 0.1);
    pp_phi[t] = normal_rng(phi_a*pse[2,t]+phi_b, 0.3);
  }
}

