require(pomp)
require(ggplot2)
require(tidyverse)
require(extraDistr)

hmm_step <- Csnippet("
  if (X == 0) {
    X = rbinom(1, p0);
  } else {
    X = rbinom(1, p1);
  }
")

hmm_rinit <- Csnippet("
  X = 0; 
")

hmm_dmeas <- Csnippet("
  if (X == 0) {
    lik = dlst(Y, 2, mu0, s0, give_log);
  } else {
    lik = dlst(Y, 2, mu1, s1, give_log);
  }
")

hmm_rmeas <- Csnippet("
  if (X == 0) {
    Y = rlst(1, 2, mu0, s0);
  } else {
    Y= rlst(1, 2, mu1, s1);
  }
")

source("toydata.R")
data = generate_data(length = 800)
ts = data.frame(t = 1:length(data$Y), Y = data$Y)

hmm <- pomp(data = ts,
  times="t",t0=0,
  rprocess=euler(hmm_step,delta.t=1/12),
  rinit=hmm_rinit,
  rmeasure=hmm_rmeas,
  dmeasure=hmm_dmeas,
  statenames=c("X"),
  paramnames=c("p0","p1","mu0","mu1", "s0", "s1", "eta"),
  partrans=parameter_trans(
    log=c("s0", "s1"),
    logit=c("p0","p1")))

require(ggplot2)
ggplot(ts) + geom_line(aes(t, Y)) + geom_point(aes(t,Y))

hmm %>% simulate(
  params=c(mu0=3,mu1=5,p0=0.2,p1=0.7, s0 = 1, s1 = 1, eta = 0.5),
  nsim=3,format="data.frame",include.data=TRUE
) -> sims

sims %>%
  ggplot(aes(x=t,y=Y,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)