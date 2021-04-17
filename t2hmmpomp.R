require(extraDistr)

hmm_step <- function (X, p0, p1, s0, s1, delta.t, ...)
{
  if (X == 0) {
    X = rbinom(1, 1, p0)
  } else {
    X = rbinom(1, 1, p1)
  }

  c(X = X)
}

hmm_rinit <- function (eta, ...) {
  c(X = as.numeric(rbinom(1, 1, eta)))
}

hmm_dmeas <- function (Y, X, mu0, mu1, s0, s1, log, ...) {
  if (X == 0) {
    dlst(Y, 2, mu0, s0, log)
  } else {
    dlst(Y, 2, mu1, s1, log)
  }
}

hmm_rmeas <- function (X, mu0, mu1, s0, s1, ...) {
  if (X == 0) {
    c(Y=rlst(1, 2, mu0, s0))
  } else {
    c(Y=rlst(1, 2, mu1, s1))
  }
}

source('toydata.R')
data= generate_data()
ts = data.frame(t = 1:length(data$Y), Y = data$Y)
require(pomp)

ts %>% 
  pomp(
    times="t",t0=0,
    rprocess=euler(hmm_step,delta.t=1/12),
    rinit=hmm_rinit,
    rmeasure=hmm_rmeas,
    dmeasure=hmm_dmeas,
    statenames=c("X"),
    paramnames=c("p0","p1","mu0","mu1", "s0", "s1", "eta"),
    partrans=parameter_trans(
      log=c("s0", "s1"),
      logit=c("p0","p1"))
  ) -> hmm

hmm %>% simulate(
  params=c(mu0=0,mu1=3,p0=0.2,p1=0.8,s0 = 0.4, s1 = 0.4, eta = 0.5),
  nsim=20,format="data.frame",include.data=TRUE
) -> sims

sims %>%
  ggplot(aes(x=t,y=Y,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)
