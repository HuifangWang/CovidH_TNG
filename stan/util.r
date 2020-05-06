library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

pq <- function(y, lim, col="black", isfirst=FALSE) {
  q <- apply(y, 2, function(x) { quantile(x, c(0.025, 0.5, 0.975), na.rm=T) });
  if (isfirst) {
    plot(q[1,], type="l", ylim=lim, lty=1, col=col);
  } else {
    lines(q[1,], type="l", ylim=lim, lty=1, col=col);
  }
  # lines(q[2,], type="l", ylim=lim, lty=1, col=col);
  lines(q[3,], type="l", ylim=lim, lty=1, col=col);
}

use_stan_funcs <- function(fname) {
  expose_stan_functions(stanc_builder(fname))
}