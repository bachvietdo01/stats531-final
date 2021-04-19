require(pomp)
require(ggplot2)
require(tidyverse)

data_ts <- read.csv("xrayts.csv")

ghmmar_filt <- Csnippet("
  if (X == 0) {
    X = rbinom(1, p0);
  } else {
    X = rbinom(1, p1);
  }
  
  mu0 = a0 + b0 * Y_state;
  mu1 = a1 + b1 * Y_state;
  
  Y_state = covaryt;
")

ghmmar_sim <- Csnippet("
  if (X == 0) {
    X = rbinom(1, p0);
  } else {
    X = rbinom(1, p1);
  }
  
  mu0 = a0 + b0 * Y_state;
  mu1 = a1 + b1 * Y_state;
  
  if (X == 0) {
    Y_state = rnorm(mu0, s0);
  } else {
    Y_state = rnorm(mu1, s1);
  }
")

ghmmar_rinit <- Csnippet("
  X = x0;
  mu0 = mu0_init; 
  mu1 = mu1_init;
  
  if (X == 0) {
    Y_state = mu0;
  } else {
    Y_state = mu1;
  }
")

ghmmar_dmeas <- Csnippet("
  if (X == 0) {
    lik = dnorm(y, mu0, s0, give_log);
  } else {
    lik = dnorm(y, mu1, s1, give_log);
  }
")

ghmmar_rmeas <- Csnippet("
  y = Y_state;
")

hmm.filt <- pomp(data=data_ts,
                 statenames=c("X","Y_state", "mu0", "mu1"),
                 paramnames=c("a0","b0","a1","b1","s0","s1","p0","p1",
                              "x0","mu0_init","mu1_init"),
                 times="t",
                 t0=0,
                 covar=covariate_table(
                   time=0:length(data_ts$y),
                   covaryt=c(0,data_ts$y),
                   times="time"),
                 rmeasure=ghmmar_rmeas,
                 dmeasure=ghmmar_dmeas,
                 rprocess=discrete_time(ghmmar_filt,
                                        delta.t=1),
                 rinit=ghmmar_rinit,
                 partrans=parameter_trans(
                   log=c("s0", "s1"),
                   logit=c("p0","p1"))
)

hmm.sim <- pomp(hmm.filt, 
                statenames=c("X","Y_state", "mu0", "mu1"),
                paramnames=c("a0","b0","a1","b1","s0","s1","p0","p1", 
                             "x0", "mu0_init", "mu1_init"),
                rprocess=discrete_time(
                  step.fun=ghmmar_sim,delta.t=1))

# initial values
params=c(p0=0.2,p1=0.7, a0 = -6, b0 = 0.1, a1 = -4.5, b1 = 0.1, s0 = 0.4, 
         s1 = 0.4)
fixed_params = c(x0 = 0, mu0_init = -6, mu1_init = -4.5)


sim1.sim <- simulate(hmm.sim,seed=1,params=c(params, fixed_params))

sim1.filt <- pomp(sim1.sim, 
                  covar=covariate_table(
                    time=c(timezero(sim1.sim),time(sim1.sim)),
                    covaryt=c(obs(sim1.sim),NA),
                    times="time"),
                  statenames=c("X","Y_state", "mu0", "mu1"),
                  paramnames=c("a0","b0","a1","b1","s0","s1","p0","p1", 
                               "x0", "mu0_init", "mu1_init"),
                  rprocess=discrete_time(
                    step.fun=Csnippet(ghmmar_filt),delta.t=1)
)

hmm.sim %>% simulate(
  params=c(p0=0.2,p1=0.7, a0 = -6, b0 = 0.1, a1 = -4.5, b1 = 0.1, s0 = 0.4, 
           s1 = 0.4, x0 = 0, mu0_init = -6, mu1_init = -4.5),
  nsim=2,format="data.frame",include.data=TRUE
) -> sims

sims %>%
ggplot(aes(x=t,y=y,group=.id,color=.id=="data"))+ 
  geom_line()+  theme_bw()

# perform main computation
library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()  
registerDoParallel(cores)
library(doRNG)
registerDoRNG(34118892)

run_level <- 3
Np <-           switch(run_level, 1e3, 1e4, 2e4)
Nmif <-         switch(run_level,  10, 100, 200)
Nreps_eval <-   switch(run_level,   4,  10,  20)
Nreps_local <-  switch(run_level,  10,  20,  20)
Nreps_global <- switch(run_level,  10,  20, 100)

# compute log likelihood
tic <- Sys.time()
stew(file=sprintf("results/arghmm2-%d.rda",run_level),{
  t.pf1 <- system.time(
    pf1 <- foreach(i=1:Nreps_eval,
                   .packages='pomp') %dopar% pfilter(sim1.filt,Np=Np, params = c(params, fixed_params)))
})
(L.pf1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))
toc <- Sys.time()
# write result to csv file
pf1[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik=L.pf1[1],loglik.se=L.pf1[2]) %>%
  write_csv("arghmm2_params.csv")


arghmm.rw <- rw.sd(a0 = 0.1, 
                   a1 = 0.1,
                   b0 = 0.01,
                   b1 = 0.01,
                   p0 = 0.01, 
                   p1 = 0.01, 
                   s0 = 0.01, 
                   s1 = 0.01)

# perform local search
stew(file=sprintf("results/arghmm2_local_search-%d.rda",run_level),{
  t.if1 <- system.time({
    if1 <- foreach(i=1:Nreps_local,
                   .packages='pomp', .combine=c) %dopar% mif2(hmm.filt,
                                                              params=c(params,
                                                                  fixed_params),
                                                              Np=Np,
                                                              Nmif=Nmif,
                                                              cooling.fraction.50=0.5,
                                                              rw.sd = arghmm.rw)
    L.if1 <- foreach(i=1:Nreps_local,
                     .packages='pomp', .combine=rbind) %dopar% logmeanexp(
                       replicate(Nreps_eval, logLik(pfilter(hmm.filt, params=coef(if1[[i]]),Np=Np))), se=TRUE)
  })
})
toc <- Sys.time()

# visualize local search results
if1 %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")

r.if1 <- data.frame(logLik=L.if1[,1],logLik_se=L.if1[,2],
                    t(sapply(if1,coef)))

if (run_level>1) write.table(r.if1,file="arghmm2_params.csv",
                             append=TRUE,col.names=FALSE,row.names=FALSE,
                             sep=",")

arghmm2_box <- rbind(
  a0 = c(-6.5, -4.0),
  a1 = c(-7, -3.0),
  b0 = c(-0.3, 0.3),
  b1 = c(-0.5, 0.3),
  p0 = c(0.01,0.25),
  p1 = c(0.5,0.8),
  s0 = c(0.01, 0.2),
  s1 = c(0.3, 0.6)
)

tic <- Sys.time()
stew(file=sprintf("results/arghmm2box_eval-%d.rda",run_level),{
  t.box <- system.time({
    if.box <- foreach(i=1:Nreps_global,
                      .packages='pomp',.combine=c) %dopar% mif2(if1[[1]],
    params=c(unlist(apply(arghmm2_box,1,function(x)runif(1,x[1], x[2]))),fixed_params))
    L.box <- foreach(i=1:Nreps_global,
                     .packages='pomp',.combine=rbind) %dopar% {
                       logmeanexp(replicate(Nreps_eval, logLik(pfilter(
                         hmm.filt,params=coef(if.box[[i]]),Np=Np))), 
                         se=TRUE)}
  })
})
toc <- Sys.time()

r.box <- data.frame(logLik=L.box[,1],logLik_se=L.box[,2],
                    t(sapply(if.box,coef)))
if(run_level>1) write.table(r.box,file="arghmm2_params.csv",
                            append=TRUE,col.names=FALSE,row.names=FALSE,
                            sep = ",")
summary(r.box$logLik,digits=5)

if.box %>%
  traces() %>%
  melt()  %>%
  droplevels() %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=L1))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  guides(color=FALSE)


if.box %>% 
  as("data.frame") %>% 
  tidyr::gather(variable,value,-t,-.id) %>%
  dplyr::filter(variable == c("cond.logLik","ess")) %>%
  ggplot(aes(x=t,y=value,group=.id,color=.id))+
  geom_line(color='antiquewhite4')+
  facet_wrap(~variable,scales="free_y",ncol=1)+
  guides(color=FALSE)


pairs(~logLik+a0+b0+a1+b1+s0+s1+p0+p1,
      data=subset(r.box,logLik>max(logLik)-20))

require(GGally)
r.box %>%
  dplyr::filter(is.na(logLik) | logLik>max(logLik,na.rm=TRUE)-25) %>%
  ggpairs(columns = c("logLik", "a0", "b0", "a1", "b1", "s0", "s1", "p0", "p1"),
          upper = list(continuous = wrap("points", color = "antiquewhite4")),
          lower = list(continuous = wrap(ggally_cor, alignPercent = 0.8, 
                                         color = "black")),
          diag = list(continuous = wrap("densityDiag", color = "antiquewhite4")))
