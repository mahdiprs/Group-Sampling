#########################################################################
# This script performs the calculations for Section 5.5 of the paper
# (simulation study)
# 26/10/2021
#########################################################################

#########################################################################
# Define functions for calculating densities, likelihoods and other
# quantities
#########################################################################


dty.unclust.oneval <- function(ty,b,Nbar,alpha,beta,R=1000){
  # returns P[ty] (where ty must be a single value)
  #  where ty|p ~ bin(b,1-(1-p)^Nbar) and p~Beta(alpha,beta)
  #if(ty==0) return(exp(lbeta(alpha,beta+Nbar*b)-lbeta(alpha,beta)))
  if(ty==0) return(beta.ratio(alpha=alpha,beta=beta,bplus=Nbar*b))
  if(ty==1) return(b*beta.diff.ratio(alpha=alpha,beta=beta,
                                     add1=Nbar*(b-1),add2=Nbar*b))
  if(ty>1){
    pvals <- qbeta(c(1:R)/(R+1),alpha,beta)
    return(
      choose(b,ty)*mean((1-(1-pvals)^Nbar)^ty*(1-pvals)^(Nbar*(b-ty)))
    )
  }
}

dty.unclust <- function(ty,b,Nbar,alpha,beta,R=1000){
  # returns a vector P[ty] (for each value in ty)
  #  where ty|p ~ bin(b,1-(1-p)^Nbar) and p~Beta(alpha,beta)
  # ty, b and Nbar may be vectors
  if(length(b)==1) b <- rep(b,length(ty))
  if(length(Nbar)==1) Nbar <- rep(Nbar,length(ty))
  out <- rep(NA,length(ty))
  for(j in c(1:length(ty))) out[j]<-dty.unclust.oneval(ty[j],b[j],Nbar[j],alpha,beta,R=R)
  out
}

dty.clust <- function(ty,b,Nbar,alpha,beta,theta,R=1000){
  # returns a vector P[ty] (for each value in ty)
  #  where tyi|pi ~ bin(b,1-(1-pi)^Nbar), pi|p~Beta(theta,theta*(1/p-1)) and
  #     p~Beta(alpha,beta)
  # returns P[ty] for each value in ty
  pvals <- qbeta(c(1:R)/(R+1),alpha,beta)
  log.one.minus.phi.vals <- log1p(-pvals)
  for(i in c(1:(Nbar-1))){
    log.one.minus.phi.vals <- log.one.minus.phi.vals + log1p(-theta/(theta/pvals+i))
  }
  log.phi.vals <- Rmpfr::log1mexp(-log.one.minus.phi.vals)
  prob.ty <- rep(NA,length(ty))
  phi.vals <- exp(log.phi.vals)
  one.minus.phi.vals <- exp(log.one.minus.phi.vals)
  for(k in c(1:length(ty))){
    prob.ty[k] <- mean(choose(b,ty[k])*phi.vals^ty[k]*one.minus.phi.vals^(b-ty[k]))
  }
  return(prob.ty)
}

dty <- function(ty,b,Nbar,alpha,beta,theta=Inf,R=1000){
  if(theta==Inf) return(dty.unclust(ty=ty,b=b,Nbar=Nbar,
                                    alpha=alpha,beta=beta,R=R))
  if(theta!=Inf) return(dty.clust(ty=ty,b=b,Nbar=Nbar,theta=theta,
                                  alpha=alpha,beta=beta,R=R))
}

negll <- function(par,ty,b,Nbar,freq,theta=Inf,R=1000){
  # returns the negative of the log-likelihood of a vector of values of ty
  #  where ty|p ~ bin(b,1-(1-p)^Nbar) and p~Beta(alpha,beta)
  # ty, b and Nbar may be vectors
  # par=(log(alpha),log(beta))
  # freq is the number of times each value of ty appears in dataset
  #  (if freq left blank, it is assumed that each value in ty represents a
  #    single observation)
  # (theta=Inf represents no clustering within groups,
  #   theta=0 represents perfect clustering within groups)
  # (note that alpha,beta control between-population variation i.e.
  #   clustering within consignments)
  if(missing(freq)) freq <- rep(1,length(ty))
  alpha <- exp(par[1])
  beta <- exp(par[2])
  -sum(freq*log(dty(ty=ty,b=b,Nbar=Nbar,alpha=alpha,beta=beta,theta=theta)))
}

posterior.Tx <- function(alpha,beta,ty,b,Nbar,B,R=1e5,topcat=10,mu,shape){
  # returns the posterior distribution of Tx given an observed value of ty
  # for the model without clustering
  # calculates P[Tx=k]=ty for each k=1,...,(Nbar*B)
  #   (values >= topcat are grouped together)
  # the prior may be supplied either as (alpha,beta) or as (mu,shape)
  if(missing(alpha)){
    alpha<-mu*shape
    beta<-(1-mu)*shape
    if(shape==Inf) sim.p <- rep(mu,R)
  } else shape<-alpha+beta
  if(shape<Inf) sim.p <- rbeta(R,alpha+ty,beta+Nbar*(b-ty))
  wt <- ((1-(1-sim.p)^Nbar)/sim.p)^ty
  tX.nonsamp <- rbinom(R,size=(B-b)*Nbar,prob=sim.p)
  if(ty==0) tX.y1 <- 0 else{
    Xi.y1 <- matrix(NA,R,ty)
    probzero.yi <- (1-sim.p)^Nbar
    for(i in c(1:ty)){
      U <- runif(R,min=probzero.yi,max=1)
      Xi.y1[,i] <- qbinom(U,size=Nbar,prob=sim.p)
    }
    tX.y1 <- apply(Xi.y1,1,sum)
  }
  tX <- tX.nonsamp+tX.y1
  tX <- pmin(tX,topcat)
  tX.fac <- factor(tX,levels=c(0:topcat),
                   labels=c(as.character(c(0:(topcat-1))),paste0(topcat,"+")))
  out <- tapply(wt,tX.fac,FUN="sum")/sum(wt)
  out[is.na(out)] <- 0
  return(out)
}

dsn.Tx.clustered <- function(alpha,beta,theta,b,Nbar,B,R=1e5,topcat=10){
  # returns the posterior distribution of Tx given an observed value of ty
  # for the model with clustering within groups, controlled by parameter theta
  # (different numerical method used than in dsn.Tx.clustered)
  # calculates P[Tx=k]=ty for each k=1,...,(Nbar*B)
  #   (values >= topcat are grouped together)
  # the prior must be supplied as (alpha,beta,theta)
  # (theta=Inf represents no clustering within groups,
  #   theta=0 represents perfect clustering within groups)
  # (note that alpha,beta control between-population variation i.e.
  #   clustering within consignments)
  sim.p <- rbeta(R,alpha,beta)
  sim.p.group <- sim.Xi <- matrix(NA,nrow=R,ncol=B)
  for(r in c(1:R)){
    if((theta>0)&(theta<Inf)) sim.p.group[r,] <- rbeta(B,theta,theta*(1/sim.p[r]-1))
    if(theta==0) sim.p.group[r,] <- 1*(runif(B)<=sim.p[r])
    if(theta==Inf) sim.p.group[r,] <- sim.p[r]
    sim.Xi[r,] <- rbinom(B,Nbar,sim.p.group[r,])
  }
  tx.samp <- apply(sim.Xi[,1:b],1,sum)
  tx.nonsamp <- apply(sim.Xi[,(b+1):B],1,sum)
  Tx <- apply(sim.Xi,1,sum)
  ty.samp <- apply(sim.Xi[,1:b]>=1,1,sum)
  Tx.post <- Tx[ty.samp==0]
  Tx <- pmin(Tx,topcat)
  Tx.post <- pmin(Tx.post,topcat)
  Tx.fac <- factor(Tx,levels=c(0:topcat),
                   labels=c(as.character(c(0:(topcat-1))),paste0(topcat,"+")))
  Tx.post.fac <- factor(Tx.post,levels=c(0:topcat),
                        labels=c(as.character(c(0:(topcat-1))),paste0(topcat,"+")))
  prior <- table(Tx.fac)/length(Tx.fac)
  posterior <- table(Tx.post.fac)/length(Tx.post.fac)
  list(prior=prior,posterior=posterior)
}

summary.Tx.clustered <- function(alpha,beta,theta,b,Nbar,B,R=1e5){
  # returns summary information on the prior distribution of Tx and leakage
  # (i.e. not conditional on ty)
  # population assumed to be rejected if ty>=1
  sim.p <- rbeta(R,alpha,beta)
  sim.p.group <- sim.Xi <- matrix(NA,nrow=R,ncol=B)
  for(r in c(1:R)){
    if((theta>0)&(theta<Inf)) sim.p.group[r,] <- rbeta(B,theta,theta*(1/sim.p[r]-1))
    if(theta==0) sim.p.group[r,] <- 1*(runif(B)<=sim.p[r])
    if(theta==Inf) sim.p.group[r,] <- sim.p[r]
    sim.Xi[r,] <- rbinom(B,Nbar,sim.p.group[r,])
  }
  tx.samp <- apply(sim.Xi[,1:b],1,sum)
  tx.nonsamp <- apply(sim.Xi[,(b+1):B],1,sum)
  Tx <- apply(sim.Xi,1,sum)
  ty.samp <- apply(sim.Xi[,1:b]>=1,1,sum)
  Tx.post <- Tx[ty.samp==0]
  L <- Tx*(tx.samp==0)
  c(mean(Tx),mean(tx.samp>0),mean(L),mean(L>0),mean(L[tx.samp==0]),mean(L[tx.samp==0]>0))
}

beta.ratio <- function(alpha,beta,bplus,mu,tot,log=FALSE){
  # calculates beta(alpha,beta+bplus)/beta(alpha,beta)
  # where bplus is a positive integer
  # works reasonably well even if alpha<<beta or alpha<<1
  #  if mu and tot are supplied rather than alpha and beta,
  #  then alpha=mu*tot and beta=mu*tot, and calculation still works if tot=0
  if(!missing(alpha)){
    tot<-alpha+beta
    mu<-alpha/(alpha+beta)
  }
  if(length(mu)==1) mu<-rep(mu,length(tot)) else{
    if(length(tot)==1) tot<-rep(tot,length(mu))
  }
  log.out<-rep(NA,length(mu))
  for(i in c(1:length(mu))){
    log.out[i]<-log1p(-mu[i])+sum(log1p(-mu[i]/(1+c(1:(bplus-1))/tot[i])))
  }
  if(log) return(log.out) else return(exp(log.out))
}

beta.diff.ratio <- function(alpha,beta,add1,add2,mu,tot,log=FALSE){
  # calculates (beta(alpha,beta+add1)-beta(alpha,beta+add2))/beta(alpha,beta)
  # where add1<add2 is a positive integer
  # works reasonably well even if alpha<<beta or alpha<<1
  #  if mu and tot are supplied rather than alpha and beta,
  #  then alpha=mu*tot and beta=mu*tot, and calculation still works if tot=0
  if(add1>=add2) stop("add1 must be strictly less than add2")
  if(!missing(alpha)){
    tot<-alpha+beta
    mu<-alpha/(alpha+beta)
  }
  if(length(mu)==1) mu<-rep(mu,length(tot)) else{
    if(length(tot)==1) tot<-rep(tot,length(mu))
  }
  log.out<-rep(NA,length(mu))
  for(i in c(1:length(mu))){
    log.out[i] <- log1p(-mu[i])+sum(log1p(-mu[i]/(1+c(1:(add1-1))/tot[i]))) + 
      Rmpfr::log1mexp(-sum(log1p(-mu[i]/(1+c(add1:(add2-1))/tot[i]))))
  }
  if(log) return(log.out) else return(exp(log.out))
}

# MLE

BB.group.model <- function(ty,b,B,Nbar,freq,theta=Inf,R=1000,startval,SE=FALSE){
  if(missing(startval)) startval<-c(0,0)
  if(all(ty==0)) return(list(alpha=0,beta=0,mu=0,rho=0,D=1,E.leak=0,prob.leak=0,
                             se.alpha=0,se.beta=0,se.mu=0,se.rho=0,se.D=0,se.par=c(0,0),
                             log.prob.leak=-Inf,pty0=1))
  o <- optim(startval,negll,ty=ty,freq=freq,b=b,Nbar=Nbar,R=R,theta=theta,
             hessian=SE)
  N <- Nbar*B
  n <- Nbar*b
  alpha <- exp(o$par[1])
  beta <- exp(o$par[2])
  mu <- alpha/(alpha+beta)
  rho <- 1/(alpha+beta+1)
  D <- 1+(Nbar-1)*rho
  
  pty0 <- exp(lbeta(alpha,beta+n)-lbeta(alpha,beta))
  #E.leak <- (N-n)*alpha/(alpha+beta+n)*beta(alpha,beta+n)/beta(alpha,beta)
  E.leak <- (N-n)*alpha/(alpha+beta+n)*beta.ratio(alpha=alpha,beta=beta,bplus=n)
  #prob.leak <- beta(alpha,beta+n)/beta(alpha,beta)-
  #  beta(alpha,beta+N)/beta(alpha,beta)
  log.prob.leak <- sum(log1p(-alpha/(alpha+beta+c(0:(n-1))))) +
    Rmpfr::log1mexp(-sum(log1p(-alpha/(alpha+beta+c(n:(N-1))))))
  prob.leak <- exp(log.prob.leak)
  
  if(SE){
    par.vcov <- solve(o$hessian)
    se.par <- sqrt(diag(par.vcov))
    se.alpha <- alpha*se.par[1]
    se.beta <- beta*se.par[2]
    dmu.dpar <- c(1,-1)*alpha*beta/(alpha+beta)^2
    se.mu <- sqrt(as.numeric(t(dmu.dpar)%*%par.vcov%*%dmu.dpar))
    drho.dpar <- -c(alpha,beta)/(alpha+beta+1)^2
    se.rho <- sqrt(as.numeric(t(drho.dpar)%*%par.vcov%*%drho.dpar))
    se.D <- (Nbar-1)*se.rho
    return(list(optim.results=o,alpha=alpha,beta=beta,mu=mu,rho=rho,D=D,
                se.alpha=se.alpha,se.beta=se.beta,se.mu=se.mu,se.rho=se.rho,
                se.D=se.D,se.par=se.par,
                E.leak=E.leak,prob.leak=prob.leak,log.prob.leak=log.prob.leak,
                pty0=pty0))
  }
  return(list(optim.results=o,alpha=alpha,beta=beta,mu=mu,rho=rho,D=D,
              E.leak=E.leak,prob.leak=prob.leak,log.prob.leak=log.prob.leak,
              pty0=pty0))
}

# test:
BB.group.model(ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,B=1000,Nbar=100,
               theta=Inf,SE=TRUE)

BB.group.model(ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,B=1000,Nbar=100,
               theta=0,SE=TRUE)

#########################################################################
# Function to perform simulation for a given set of parameters
#########################################################################

simfn <- function(k,parvalues,seed,R,summarise=TRUE){
  # set seed, set up parameters
  if(!missing(seed)) set.seed(seed)
  Nbar=parvalues$Nbar[k]
  b=parvalues$b[k] ; B=parvalues$B[k]
  alpha=parvalues$alpha[k] ; mu=parvalues$mu[k]
  beta=alpha*(1-mu)/mu
  rho=1/(alpha+beta+1)
  D=1+(Nbar-1)*rho
  theta=parvalues$theta[k]
  num.consignments=parvalues$num.consignments[k]
  N=Nbar*B ; n=Nbar*b
  if(theta==Inf){
    pty0 <- exp(lbeta(alpha,beta+n)-lbeta(alpha,beta))
    E.leak <- (N-n)*alpha/(alpha+beta+n)*beta(alpha,beta+n)/beta(alpha,beta)
    log.prob.leak <- sum(log1p(-alpha/(alpha+beta+c(0:(n-1))))) +
      Rmpfr::log1mexp(-sum(log1p(-alpha/(alpha+beta+c(n:(N-1))))))
    prob.leak <- exp(log.prob.leak)
  }
  if(theta==0){
    pty0 <- exp(lbeta(alpha,beta+b)-lbeta(alpha,beta))
    E.leak <- (N-n)*alpha/(alpha+beta+b)*beta(alpha,beta+b)/beta(alpha,beta)
    log.prob.leak <- sum(log1p(-alpha/(alpha+beta+c(0:(b-1))))) +
      Rmpfr::log1mexp(-sum(log1p(-alpha/(alpha+beta+c(n:(B-1))))))
    prob.leak <- exp(log.prob.leak)
  }
  if((theta>0)&(theta<Inf)){
    p.i.vals<-qbeta(c(1:1000)/1001,alpha,beta)
    ratio.vals <- beta.ratio(alpha=theta,beta=theta*(1/p.i.vals-1),bplus=Nbar)
    pty0<-mean(ratio.vals^b)
    E.leak<-mean((N-n)*p.i.vals*ratio.vals^b)
    prob.leak<-mean(ratio.vals^b*(1-ratio.vals^(B-b)))
  }
  # object to save results
  #simres <- data.frame(alpha=rep(alpha,R),beta=beta,mu=mu,theta=theta,
  #                     a=NA,b=NA,l=NA,L=NA,p.leak=NA,P.leak=NA)
  simres <- data.frame(a=rep(NA,R))
  # loop
  for(r in c(1:R)){
    all.consignments<-data.frame(i=rep(c(1:num.consignments),each=B),
                                 j=rep(c(1:B),times=num.consignments),
                                 p.i=NA,p.ij=NA)
    all.consignments$p.i <- rep(rbeta(num.consignments,alpha,beta),
                                each=B)
    if((theta<Inf)&(theta>0)) all.consignments$p.ij <-
      rbeta(B*num.consignments,rep(theta,B*num.consignments),
            theta*(1/all.consignments$p.i-1))
    if(theta==Inf) all.consignments$p.ij <- all.consignments$p.i
    if(theta==0) all.consignments$p.ij <- rbinom(B*num.consignments,1,
                                                 all.consignments$p.i)
    all.consignments$X.i<-rbinom(B*num.consignments,Nbar,all.consignments$p.ij)
    ty<-with(all.consignments[all.consignments$j<=b,],tapply(X.i>=1,i,FUN=sum))
    Tx <- with(all.consignments,tapply(X.i,i,FUN=sum))
    tab.ty <- table(ty)
    modfit <- BB.group.model(ty=as.numeric(names(tab.ty)),
                               freq=as.numeric(tab.ty),b=b,B=B,Nbar=Nbar,
                               theta=Inf,SE=TRUE)
    simres[r,c("a","b","mu","rho","D")] <- c(modfit$alpha,modfit$beta,modfit$mu,
                                         modfit$rho,modfit$D)
    simres[r,c("se.a","se.b","se.mu","se.rho","se.D")] <-
      c(modfit$se.alpha,modfit$se.beta,modfit$se.mu,modfit$se.rho,modfit$se.D)
    simres[r,c("se.log.a","se.log.b")] <- modfit$se.par
    simres[r,c("prob.leak","pty0","E.leak")] <- c(modfit$prob.leak,modfit$pty0,modfit$E.leak)
    simres[r,"P.leak"] <- mean((ty==0)*(Tx>0))
    simres[r,"mean.ty"] <- mean(ty)
    simres[r,"shipments.positive"] <- sum(ty>0)
    simres[r,"max.ty"] <- max(ty)
  }
  simres$se.a[is.na(simres$se.a)] <- 0
  simres$se.b[is.na(simres$se.b)] <- 0
  simres$se.mu[is.na(simres$se.mu)] <- 0
  simres$se.D[is.na(simres$se.D)] <- 0
  simres$se.rho[is.na(simres$se.rho)] <- 0
  simres$se.log.a[is.na(simres$se.log.a)] <- 0
  simres$se.log.b[is.na(simres$se.log.b)] <- 0
  sim.summary <- c(k=k,alpha=alpha,beta=beta,mu=mu,theta=theta,rho=rho,R=R,Nbar=Nbar,
                   mean.ty=mean(simres$mean.ty),
                   prop.zero.shipments.pathogen=mean(simres$shipments.positive==0),
                   prop.one.shipments.pathogen=mean(simres$shipments.positive==1),
                   prop.ty01=mean(simres$max.ty<=1),
                   pty0=pty0,E.leak=E.leak,prob.leak=prob.leak,
                   rmse.alpha=sqrt(mean((simres$a-alpha)^2)),
                   rmse.alpha2=sqrt(mean((simres$a[simres$max.ty>=2]-alpha)^2)),
                   #rmse.beta=sqrt(mean((simres$b-beta)^2)),
                   rmse.mu=sqrt(mean((simres$mu-mu)^2)),
                   rmse.mu2=sqrt(mean((simres$mu[simres$max.ty>=2]-mu)^2)),
                   rmse.rho=sqrt(mean(simres$rho-rho)^2),
                   rmse.rho2=sqrt(mean((simres$rho[simres$max.ty>=2]-rho)^2)),
                   mad.mu=1.4826*median(abs(simres$mu-mu)),
                   mad.rho=1.4826*median(abs(simres$rho-rho)),
                   #rmse.log.alpha=sqrt(mean((log(simres$a)-log(alpha))^2)),
                   #rmse.log.beta=sqrt(mean((log(simres$b)-log(beta))^2)),
                   #rmse.D=sqrt(mean((simres$D-D)^2)),
                   rmse.pty0=sqrt(mean((simres$pty0-pty0)^2)),
                   rmse.E.leak=sqrt(mean((simres$E.leak-E.leak)^2)),
                   rmse.prob.leak=sqrt(mean((simres$prob.leak-prob.leak)^2)),
                   bias.alpha=mean(simres$a-alpha),
                   bias.alpha2=mean(simres$a[simres$max.ty>=2]-alpha),
                   #bias.beta=mean(simres$b-beta),
                   bias.mu=mean(simres$mu-mu),
                   bias.mu2=mean(simres$mu[simres$max.ty>=2]-mu),
                   bias.rho=mean(simres$rho-rho),
                   bias.rho2=mean(simres$rho[simres$max.ty>=2]-rho),
                   bias.mu.median=median(simres$mu-mu),
                   bias.rho.median=median(simres$rho-rho),
                   #bias.D=mean(simres$D-D),
                   bias.pty0=mean(simres$pty0-pty0),
                   bias.E.leak=mean(simres$E.leak-E.leak),
                   bias.prob.leak=mean(simres$prob.leak-prob.leak),
                   noncover.alpha=mean(abs((simres$a-alpha)/simres$se.a)>qnorm(0.975)),
                   noncover.alpha2=mean(abs((simres$a-alpha)/simres$se.a)[simres$max.ty>=2]>qnorm(0.975)),
                   #noncover.beta=mean(abs((simres$b-beta)/simres$se.b)>qnorm(0.975)),
                   #noncover.log.alpha=mean(abs((log(simres$a)-log(alpha))/simres$se.log.a)>qnorm(0.975)),
                   #noncover.log.beta=mean(abs((log(simres$b)-log(beta))/simres$se.log.b)>qnorm(0.975)),
                   noncover.mu=mean(abs((simres$mu-mu)/simres$se.mu)>qnorm(0.975)),
                   noncover.mu2=mean(abs((simres$mu-mu)/simres$se.mu)[simres$max.ty>=2]>qnorm(0.975)),
                   noncover.rho=mean(abs((simres$rho-rho)/simres$se.rho)>qnorm(0.975)),
                   noncover.rho2=mean(abs((simres$rho-rho)/simres$se.rho)[simres$max.ty>=2]>qnorm(0.975)),
                   #noncover.D=mean(abs((simres$D-D)/simres$se.D)>qnorm(0.975)),
                   prop.leak=mean(simres$P.leak))
  if(summarise) return(sim.summary) else
    return(list(simres=simres,sim.summary=sim.summary))
}

#########################################################################
# Actually perform the simulation over all combinations of
# parameter settings
# Note that the parallelisation may not work on all operating systems.
# The simulation is time consuming. Progress is output to a file
# called "update.txt" in the working directory.
#########################################################################

pargrid <- expand.grid(Nbar=c(100,250,500),b=100,B=1000,alpha=c(0.01,0.015,0.025,0.1,1),
                          num.consignments=100,
                          mu=c(5e-5,5e-4),theta=c(0,5,50,Inf))

library(foreach)
library(doMC)
registerDoMC(cores=6)

t1 <- proc.time()[3]
K <- nrow(pargrid)
cat("Starting simulation study \n",file="update.txt")
simresults <- foreach(k=1:K,.combine=rbind) %dopar% {
  set.seed(81749174)
  a <- simfn(k=k,parvalues=pargrid,R=1000,summarise=TRUE)
  cat(paste0("finished simulation ",k," \n"),file="update.txt",append=TRUE)
  a
}
t2 <- proc.time()[3]
t2-t1
(simresults <- as.data.frame(simresults))
save.image("all_simulation_results.rdata") # saves all results

#########################################################################
# Produce tables and figures for the paper
# (tables will be output to files in LaTex format)
#########################################################################

tabfn <- function(simresults,file.suffix){
  simresults <- simresults[order(-simresults[,"theta"],simresults[,"alpha"]),]
  simresults$theta[simresults$theta==Inf] <- "$\\infty$"

  tabdata <- simresults
  tabdata <- cbind(format(round(tabdata[,c("alpha","beta")],digits=3),nsmall=3),
                   format(round(100*tabdata[,c("rho")],digits=4),nsmall=4),
                   tabdata$theta,
                   format(round(tabdata[,c("rmse.mu")]*1e5,digits=3),nsmall=3),
                   format(round(100*tabdata[,c("rmse.rho")],digits=4),nsmall=4),
                   format(round(tabdata[,c("bias.mu")]*1e5,digits=3),nsmall=3),
                   format(round(100*tabdata[,c("bias.rho")],digits=4),nsmall=4),
                   format(round(100*tabdata[,c("noncover.mu","noncover.rho")],digits=1),nsmall=1))
  colnames(tabdata) <- c("alpha","beta","rho","theta","rmse.mu","rmse.rho","bias.mu","bias.rho",
                         "noncover.mu","noncover.rho")
  tabdata[c(4,8),ncol(tabdata)] <- paste0(tabdata[c(4,8),ncol(tabdata)],"\\vspace{2mm}")

  tabdata2 <- simresults
  tabdata2 <- cbind(format(round(tabdata2[,c("alpha","beta")],digits=3),nsmall=3),
                    format(round(100*tabdata2[,c("rho")],digits=4),nsmall=4),
                    tabdata2$theta,
                    format(round(tabdata2[,c("rmse.mu2")]*1e5,digits=3),nsmall=3),
                    format(round(100*tabdata2[,c("rmse.rho2")],digits=4),nsmall=4),
                    format(round(tabdata2[,c("bias.mu2")]*1e5,digits=3),nsmall=3),
                    format(round(100*tabdata2[,c("bias.rho2")],digits=4),nsmall=4),
                    format(round(100*tabdata2[,c("noncover.mu2","noncover.rho2")],digits=1),nsmall=1),
                    format(round(100*tabdata2[,"prop.zero.shipments.pathogen"]+tabdata2[,"prop.one.shipments.pathogen"],
                                 digits=1),nsmall=1))
  colnames(tabdata2) <- c("alpha","beta","rho","theta","rmse.mu","rmse.rho","bias.mu","bias.rho",
                          "noncover.mu","noncover.rho","pct.excl.")
  tabdata2[c(4,8),ncol(tabdata)] <- paste0(tabdata2[c(4,8),ncol(tabdata2)],"\\vspace{2mm}")

  tabdata3 <- simresults
  tabdata3 <- cbind(format(round(tabdata3[,c("alpha","beta")],digits=3),nsmall=3),
                    format(round(100*tabdata3[,c("rho")],digits=4),nsmall=4),
                    tabdata3$theta,
                    format(round(tabdata3[,c("E.leak")],digits=3),nsmall=3),
                    format(round(tabdata3[,c("prob.leak")],digits=7),nsmall=7),
                    format(round(tabdata3[,c("rmse.E.leak")],digits=3),nsmall=3),
                    format(round(tabdata3[,c("rmse.prob.leak")],digits=7),nsmall=7),
                    format(round(tabdata3[,c("bias.E.leak")],digits=3),nsmall=3),
                    format(round(tabdata3[,c("bias.prob.leak")],digits=3),nsmall=3))
  colnames(tabdata3) <- c("alpha","beta","rho","theta","E.leak","prob.leak",
                          "rmse.E.leak","rmse.prob.leak","bias.E.leak","bias.prob.leak")
  
  write.table(tabdata[,-2],file=paste0("simtab1",file.suffix,".txt"),
              sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
  write.table(tabdata2[,-2],file=paste0("simtab2",file.suffix,".txt"),
              sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
  write.table(tabdata3[,-2],file=paste0("simtab3",file.suffix,".txt"),
              sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
  
  list(tabdata1=tabdata,tabdata2=tabdata2,tabdat3=tabdata3)
}

tabfn(simresults[with(simresults,(theta %in% c(Inf,5,0)&(mu==5e-5)&(Nbar==100)&
                                    (alpha %in% c(0.01,0.015,0.025,0.1)))),],
      file.suffix="_26_10_2021")

for(Nbar.val in unique(simresults$Nbar)){
  for(mu.val in unique(simresults$mu)){
    tabfn(simresults[with(simresults,(mu==mu.val)&(Nbar==Nbar.val)&(alpha!=1)),],
          file.suffix=paste0("Nbar",Nbar.val,"_mu",mu.val,"_26_10_2021"))
  }
}

