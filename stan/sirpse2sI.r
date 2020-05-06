source('util.r')
library(tibble)
library(ggplot2)
use_stan_funcs("sirpse2sI.stan")

# setup simulation
noise = 0.5;
nt = 210; # total days
no = 90; # observed days
gamma = 0.2;
Gamma = c(0.5, 0.1);
dwI = mat_rng(2, nt - 1);
arI = mat_rng(5,1)[,1] / 10;
yt = run(gamma, Gamma, noise, dwI, arI);
ur = 2.0; # under report
lur = log(ur);
lur_sd = lur;
oI = yt[3:4,] + lur;
soI = log(exp(oI[1,]) + exp(oI[2,]));

# priors on parameters
gamma_mu = 0.2;
gamma_sd = 1;
Gamma_mu = c(0.1, 0.5);
Gamma_sd = c(1, 1);

# run fit
data = list(
  noise=noise,
  nt=nt,
  no=no,
  soI=soI,
  lur=lur,
  lur_sd=lur_sd,
  gamma_mu=gamma_mu,
  gamma_sd=gamma_sd,
  Gamma_mu=Gamma_mu,
  Gamma_sd=Gamma_sd,
  pseh_sd=0.3,
  use_pse=1,
  pse=sub_pse(yt)
)
fitI<-stan("sirpse2sI.stan", chains=4, data=data, control=list(max_treedepth=3))
pairs(fitI, pars=c("gammah", "Gammah"))
summary(fitI, pars=c("gammah", "Gammah"))

# plot results
pdln <- function(mu,sd) { x<-seq(0,3,len=100);lines(x,dnorm(x,mean=mu,sd=sd)) }
ef <- extract(fitI)
Ip <- exp(ef$Ip)
sIp <- Ip[,1,] + Ip[,2,];
psep <- ef$psep;

par(mfrow=c(3,3))
hist(ef$gamma, probability = T, main=expression(gamma));pdln(gamma_mu, gamma_sd);abline(v=gamma, col="red", lwd=3);
hist(ef$Gamma[,1], probability = T, main=expression(Gamma[1]));pdln(Gamma_mu[1], Gamma_sd[1]);abline(v=Gamma[1], col="red", lwd=3);
text(0.75,1.6,"True value")
text(0.85,0.9,"Posterior")
text(1.3,0.5,"Prior")
hist(ef$Gamma[,2], probability = T, main=expression(Gamma[2]));pdln(Gamma_mu[2], Gamma_sd[2]);abline(v=Gamma[2], col="red", lwd=3);

Ilim=c(0, 0.2)
pq(Ip[,1,], lim=Ilim, isfirst = TRUE); abline(v=90, col="blue");
lines(exp(oI[1,]), type="l", col="red"); title(expression(I[1]))
text(40,0.15,"Observed")
text(150,0.15,"Predicted")
xlab("Days from 02/14")

pq(Ip[,2,], lim=Ilim, isfirst = TRUE); abline(v=90, col="blue");
lines(exp(oI[2,]), type="l", col="red"); title(expression(I[2]))
xlab("Days from 02/14")

pq(sIp, lim=Ilim, isfirst = TRUE); abline(v=90, col="blue");
lines(exp(soI), type="l", col="red"); title(expression(I[1]+I[2]))
xlab("Days from 02/14")

plim=c(-0.1,3)
psetitles=list(
  expression(beta[1]+beta[2]),
  expression(phi[1]+phi[2]),
  expression(u[1]+u[2])
)
pse <- sub_pse(yt)
for (i in 1:3) {
  pq(psep[,i,], lim=plim, isfirst = TRUE);
  lines(pse[i,], type="l", col="red"); title(psetitles[i]);
  abline(v=90/7, col="blue");
  xlab("Weeks from 02/14")
}
