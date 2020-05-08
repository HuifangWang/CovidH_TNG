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

  int count_nan_mat(matrix m) {
    int n = 0;
    for (c in 1:cols(m))
      for (r in 1:rows(m))
        if(is_nan(m[r,c]))
          n += 1;
    return n;
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
  
  matrix run(vector gamma, vector Gamma, real noise, matrix dwI_, vector arI) {
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
	  { // generalized g
	    // VJ "g0=g2=0, g1, g3 and g4 non-zero ...  g1=g3*Beta0, and g4=-gamma."
	    // SP "so g4>0  g2<0 g0 could be both signs, g3<0  and g1<0"
	    // gamma[] here is std_normal, we xfm here:
	    real beta_dev = y_[g,3] - beta0;
	    real phi_dev = y_[g,4] - phi0;
	    real g0 = gamma[1];
	    real g1 = gamma[2]*Jb;
	    real g2 = gamma[3]*phi_dev;
	    real g3 = gamma[4]*Jb*beta_dev;
	    real g4 = gamma[5]*phi_dev*beta_dev;
	    k[i][g,3] = -beta_dev + g0 + g1 + g2 + g3 + g4;
	  }
          // k[i][g,3] = -(y_[g,3]-beta0)*(1.0 + gamma*(y_[g,4]-phi0)) - Jb*y_[g,3];
          k[i][g,4] = ((-(y_[g,4]-(phi0-h*Jb)))-c*fu(Jb, y_[g,5]))/tau;
          k[i][g,5] = 1/sigma*(- Gamma[g]*(y_[g,3]-beta0)+u0*(Jb>0));
        }
        // print("t ", t, " k i ", i, k[i]);
        // reject_nan_mat(k[i]);
      }
      // print("y0 ", y);
      y += (dt/6)*(k[1] + 2*k[2] + 2*k[3] + k[4]); // rk4
      // print("y1 ", y);
      // reject_nan_mat(y);
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

  vector normal_vec_rng(int n) {
    vector[n] out;
    for (i in 1:n)
      out[i] = normal_rng(0, 1);
    return out;
  }

  vector normal_pos_vec_rng(int n) {
    vector[n] out = normal_vec_rng(n);
    for (i in 1:n)
      out[i] = fabs(out[i]);
    return out;
  }

  vector gamma_xfm(vector gz) {
    // SP "so g4>0  g2<0 g0 could be both signs, g3<0  and g1<0"
    // return [g0, g1, g2, g3, g4]
    vector[4] eg = exp(gz[2:5]) * 0.1;
    eg[1:3] *= -1;
    return append_row([gz[1]*0.1]', eg);
  }

  vector Gamma_xfm(vector Gz) {
    return [0.1, 0.5]' .* exp(Gz);
  }

  // cf reduce_sum
  real sbc_lp_part(int[] rs, int r0, int r1, vector[] gammah_z, vector[] Gammah_z, matrix[] dwIh, vector[] arIh, int no, real noise, vector[] soI, matrix[] pse, real ur, real ur_sd, real pseh_sd) {
    real lp = 0;
    for (r in r0:r1) {
      if (rs[r] == -1) // nans in yt
	continue;
      // ex transformed parameters
      vector[5] gammah = gamma_xfm(gammah_z[r]);
      vector[2] Gammah = Gamma_xfm(Gammah_z[r]);
      matrix[10,no] yth = run(gammah, Gammah, noise, dwIh[r], arIh[r]);
      matrix[3,no/7] pseh = sub_pse(yth);
      row_vector[no] soIh = exp(yth[3,]) + exp(yth[4,]);
      vector[no] urh = soIh' ./ soI[r][1:no];
      // print("gammah", gammah);
      reject_nan_mat([urh']);
      // ex model
      lp += normal_lpdf(gammah_z[r]				| 0, 1);
      lp += normal_lpdf(Gammah_z[r]				| 0, 1);
      lp += normal_lpdf(to_vector(dwIh[r])			| 0, 1);
      lp += normal_lpdf(to_vector(arIh[r])			| 0, 0.1);
      lp += normal_lpdf(to_vector(urh)				| ur, ur_sd);
      lp += normal_lpdf(to_vector(pse[r][,1:cols(pseh)])	| to_vector(pseh), pseh_sd);
    }
    return lp;
  }
}
