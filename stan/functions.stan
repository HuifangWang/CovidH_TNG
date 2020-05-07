// variant in which only (I-S) is stochastically forced with colored noise
// This implements the new equations & log(I) 05/05/2020

// TODO add genaralized g term

functions {
  matrix mat_rng(int m, int n) {
    matrix[m, n] out;
    for (i in 1:n) for (j in 1:m) out[j, i] = normal_rng(0.0, 1.0);
    return out;
  }
  
  real fu(real Jb, real u) {
    return Jb==0.0 ? 0.0 : tanh(fabs(u));
  }
  
  matrix color_dw(matrix dw, vector ar) {
    int p = rows(ar) - 1;
    matrix[rows(dw), cols(dw)] out = dw;
    for (t in (p+1):cols(dw))
      out[,t] = ar[1] + dw[,t-p:t-1] * ar[2:p+1];
    return out;
  }
  
  void reject_nan_mat(matrix m) {
    for (c in 1:cols(m))
      for (r in 1:rows(m))
        if(is_nan(m[r,c]))
          reject("mat[", r, ",", c, "] is nan!");
  }
  
  matrix sub_pse(matrix yt) {
    matrix[3,cols(yt)/7] out;
    for (t in 1:cols(out)) {
      int w = t*7;
      out[1,t] = (yt[5,w] + yt[6,w])/2;
      out[2,t] = (yt[7,w] + yt[8,w])/2;
      out[3,t] = (yt[9,w] + yt[10,w])/2;
    }
    return out;
  }
  
  matrix run(real gamma, vector Gamma, real noise, matrix dwI_, vector arI) {
    int nt = cols(dwI_) + 1;
    matrix[2, nt - 1] dwI = color_dw(dwI_, arI);
    real dt = 1.0;
    real sqrt_dt = sqrt(dt);
    real dt_10 = dt / 10.0;
    real Jb = 0.0;
    real u_mul;
    real Lambda = 1.0 / 80.0 / 365.0 / 2.0;
    real mu = 1.0 / 80.0/ 365.0;
    real alpha = 0.2;
    real phi0 = 1.0;
    real sigma = 30.0;
    real u0 = 0.1;
    real beta0 = 0.4;
    real h = 0.1;
    real tau = 1.0;
    real c = 1.0;
    real I0 = 0.0001 / 2;
    real S0 = (1.0/2) - I0;
    matrix[2, 5] y = to_matrix({S0, S0, log(I0), log(I0), 0.3, 0.3, 1,  1, 0.01, 0.01}, 2, 5); // TODO initial conditions
    matrix[2, 5] k[4];
    matrix[2, 5] y_;
    matrix[10, nt] yt;
    yt[,1] = to_vector(y);
    for (t in 2:nt) {
      real Jb_mul = Jb; // hacky Jb[t]/Jb[t-1]
      // set Jb
      if (t > 30) Jb = 1.0;
      if (t > 90) Jb = 0.75;
      if (t > 105) Jb = 0.5;
      if (t > 120) Jb = 0.25;
      if (t > 150) Jb = 0.0;
      Jb_mul = square(Jb/Jb_mul);
      u_mul = (t == 90 || t == 105 || t == 120 || t==150) ? Jb_mul : 1.0;
      if (is_nan(u_mul)) reject("u_mul nan");
      // if (is_nan(Jb_mul)) reject("Jb_mul nan");
      // 4 evals for RK4
      for (i in 1:4) {
        if (i == 1)
          y_ = y;
        else if (i == 2 || i == 3)
          y_ += dt/2*k[i-1];
        else // i == 4
          y_ += dt*k[i-1]; 
        for (g in 1:2) {
          k[i][g,1] = Lambda - y_[g,3]*sum(exp(y_[,2]))*y_[g,1] - mu*y_[g,1];
          k[i][g,2] = y_[g,3]*sum(y_[,1])-alpha-mu;
          k[i][g,3] = -(y_[g,3]-beta0)*(1.0 + gamma*(y_[g,4]-phi0)) - Jb*y_[g,3];
          k[i][g,4] = ((-(y_[g,4]-(phi0-h*Jb)))-c*fu(Jb, y_[g,5]))/tau;
          k[i][g,5] = 1/sigma*(- Gamma[g]*(y_[g,3]-beta0)+u0*(Jb>0));
        }
        // print("t ", t, " k i ", i, k[i]);
        reject_nan_mat(k[i]);
      }
      // print("y0 ", y);
      y += (dt/6)*(k[1] + 2*k[2] + 2*k[3] + k[4]); // rk4
      // print("y1 ", y);
      reject_nan_mat(y);
      y[,2] += sqrt_dt * noise * dwI[,t-1]; // noise on log(I)
      // print("y2 ", y);
      y[,5] *= u_mul; // update u on Jb change
      // print("y3 ", y);
      yt[,t] = to_vector(y); // save
    }
    return yt;
  }

  // observe I as sum of under reported infections
  vector obs_I(matrix yt, real ur) {
    row_vector[cols(yt)] sum_I = exp(yt[3,]) + exp(yt[4,]);
    return sum_I'/ur;
  }
}
