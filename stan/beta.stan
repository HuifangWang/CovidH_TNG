data {
  int nt; // 60 days starting March 10th
  vector[nt] r; // RKI R 7 day now-cast
  vector[nt] m; // Google mobility retail & recr / 100%
  int iP[8];
  vector[8] P;
}

transformed data {
  real alpha = 1.0 / 14.0;
  real a0 = (max(r)-min(r))/(max(m)-min(m));
}

parameters {
  vector[nt-1] z;
  real a;
  real l;
  real s;
}

transformed parameters {
  vector[nt] b = r*alpha;
  for (t in 1:nt-1) {
    real db = (b[1] + a*m[t]) - l*b[t];
    b[t+1] = b[t] + db + s*z[t];
  }
}

model {
  // priors
  to_vector(z) ~ std_normal();
  l ~ lognormal(1,1);
  a ~ std_normal();
  s ~ lognormal(0,1);
  // condition on data
  r ~ normal(b/alpha,s);
  // P ~ normal(p[iP],eps/49);
}
