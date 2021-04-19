require(pomp)

data_ts <- read.csv("xrayts.csv")

hmm_step <- Csnippet("
  if (X == 0) {
    X = rbinom(1, p0);
  } else {
    X = rbinom(1, p1);
  }
")

hmm_rinit <- Csnippet("
  X = rbinom(1, eta); 
")

hmm_dmeas <- Csnippet("
  if (X == 0) {
    lik = dnorm(y, mu0, sigma0, give_log);
  } else {
    lik = dnorm(y, mu1, sigma1, give_log);
  }
")

hmm_rmeas <- Csnippet("
  if (X == 0) {
    y = rnorm(mu0, sigma0);
  } else {
    y= rnorm(mu1, sigma1);
  }
")

params=c(mu0=-6,mu1=-4.5,p0=0.2,p1=0.7, sigma0 = 0.2, sigma1 = 0.2, eta = 0.5)

hmm <- pomp(data = data_ts,
            times="t",t0=0,
            rprocess=discrete_time(hmm_step),
            rinit=hmm_rinit,
            rmeasure=hmm_rmeas,
            dmeasure=hmm_dmeas,
            statenames=c("X"),
            paramnames=c("p0","p1","mu0","mu1", "sigma0", "sigma1", "eta"),
            partrans=parameter_trans(
              log=c("sigma0", "sigma1"),
              logit=c("p0","p1")
            ))


hmm %>% simulate(
  params=params,
  nsim=1,format="data.frame",include.data=TRUE
) -> sims

sims %>%
  ggplot(aes(x=t,y=y,group=.id,color=.id=="data"))+ 
  geom_line()+  
  scale_colour_grey(start = 0, end = 0.9)+ theme_bw()


library(foreach)
library(doParallel)
cores <- detectCores()  
registerDoParallel(cores)
library(doRNG)
registerDoRNG(625904618)
require(tidyverse)


tic <- Sys.time()
foreach(i=1:10,.combine=c) %dopar% {
  library(pomp)
  hmm %>% pfilter(params=params,Np=5000)
} -> pf

pf %>% logLik() %>% logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()

pf[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) %>%
  write_csv("hmm_params.csv")

tic <- Sys.time()
registerDoRNG(482947940)
bake(file="results/ghmm2_local_search.rds",{
  foreach(i=1:20,.combine=c) %dopar% {
    library(pomp)
    library(tidyverse)
    hmm %>%
      mif2(
        params=params,
        Np=200, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(mu0 = 0.01, mu1 = 0.005, p0 = 0.01, p1 = 0.01, 
                    sigma0 = 0.01, sigma1 = 0.01, eta = ivp(0.01))
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local
toc <- Sys.time()

t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")

tic <- Sys.time()
registerDoRNG(900242057)
bake(file="results//ghmm2_lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    evals <- replicate(10, logLik(pfilter(mf,Np=2000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
toc <- Sys.time()

t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+mu0+mu1+sigma0+sigma1+p0+01,data=results,pch=16)

read_csv("hmm_params.csv") %>%
  bind_rows(results) %>%
  arrange(-loglik) %>%
  write_csv("hmm_params.csv")

runif_design(
  lower=c(p0=0.1,p1=0.5,mu0=-7, mu1=-6, sigma0 = 0.01, sigma1 = 0.01, eta=0.01),
  upper=c(p0=0.5,p1=1.0,mu0=-5, mu1=-4, sigma0 = 0.7, sigma1 = 0.7, eta=1.0),
  nseq=50
) -> guesses

mf1 <- mifs_local[[1]]

bake(file="results/ghmm2_global_search.rds",{
  registerDoRNG(1270401678)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=c(unlist(guess))) %>%
      mif2(Nmif=100) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=200) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results

t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

read_csv("hmm_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("hmm_params.csv")

read_csv("hmm_params.csv") %>%
  filter(loglik>max(loglik)-50) %>%
  bind_rows(guesses) %>%
  mutate(type=if_else(is.na(loglik),"guess","result")) %>%
  arrange(type) -> all


pairs(~loglik+p0+p1+mu0+mu1+sigma0+sigma1, data=all,
      col=ifelse(all$type=="guess",grey(0.5),"red"),pch=16)