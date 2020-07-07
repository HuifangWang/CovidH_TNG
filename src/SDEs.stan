data {
  int nt; // 65 days starting March 1st
  vector[nt] r; // ICL R0(t) estimate
  vector[nt] m; // Google mobility retail & recr / 100%
  int nC;
  int iC[nC];
  vector[nC] C;
  real Cl[2];
}

transformed data {
  real alpha = 1.0 / 14.0;
  real h = 0.1;
  real sigma = 30.0;
}

parameters {
  real g;
  real G;
  real c;
  real p0;
  //real Cl[2];
}

transformed parameters {
  // init vectors at decl to avoid NaNs
  vector[nt] b = rep_vector(r[1]*alpha,nt);
  vector[nt] p = rep_vector(p0,nt);
  vector[nt] u = rep_vector(0,nt);
  // step from t to t+1
  for (t in 1:nt-1) {
    real Js = 0.4/2.0*(1-cos((2*pi()*(t - 15 + 29 + 31))/365.0));
    real db = (b[1] + (0.4*m[t]+0.05) - Js - g*(p[t] - p[1])) - b[t];
    real dp = -(p[t] - p[1] + h*m[t]) - c*tanh(2*u[t]);
    real du = -G*(b[t] - b[1]);
    b[t+1] = max([0.0, b[t] + db]);
    p[t+1] = p[t] + dp;
    u[t+1] = u[t] + du/sigma;
  }
  vector[nt-1] Ch = Cl[1]*(p[2:] - p[1:nt-1]) + Cl[2];

}

model {
  // priors
  g ~ normal(0,1); // 0.3, g=0 is null hypothesis
  G ~ normal(1,0.1); // 1
  c ~ normal(1,0.1); // 1
  p0 ~ normal(0.5,0.1);
  //Cl[1] ~ normal(10,2);
  //Cl[2] ~ normal(0.5,0.5);
  // data
  r ~ normal(b/alpha, 0.4);
  C ~ normal(Ch[iC],0.01);
}