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
  real g;
  real l;
  real s;
}

transformed parameters {
  vector[nt] b = r*alpha;
  vector[nt] p_ = rep_vector(b[1],nt);
  vector[nt] p;
  vector[nt] u = rep_vector(0,nt);
  for (t in 1:nt-1) {
    real Js = t*1.0/nt*0.13; // NY 40%, 2mo
    real db = (b[1] + a*m[t] - Js + g*u[t]) - l*b[t];
    real dp = (b[t] - p_[t] + u[t]*0.9)/3.0;
    real du = (b[1] - b[t] - u[t])/21.0;
    b[t+1] = b[t] + db + s*z[t];
    p_[t+1] = p_[t] + dp;
    u[t+1] = u[t] + du;
  }
  p = (p_ - b[1]) * 35;
}

model {
  // priors
  to_vector(z) ~ std_normal();
  l ~ lognormal(1,1);
  a ~ std_normal();
  u ~ std_normal();
  s ~ lognormal(0,1);
  g ~ lognormal(0.1,1);
  // condition on data
  r ~ normal(b/alpha,s);
  P ~ normal(p[iP],sqrt(s));
}

generated quantities {
  vector[nt] gq_r;
  vector[8] gq_P;
  for (t in 1:nt) gq_r[t] = normal_rng(b[t]/alpha,s);
  for (t in 1: 8) gq_P[t] = normal_rng(p[iP[t]], sqrt(s));
}