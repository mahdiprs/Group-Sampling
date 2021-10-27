##########################################################################
# Perform calculations and produce figures for case study in
# Sections 5.1 - 5.4 of journal paper
##########################################################################

############################################################
# Various functions to calculate densities etc
############################################################

dty.oneval <- function(ty,b,Nbar,alpha,beta,R=1000){
  # returns P[ty] (where ty must be a single value)
  #  where ty|p ~ bin(b,1-(1-p)^Nbar) and p~Beta(alpha,beta)
  if(ty==0) return(exp(lbeta(alpha,beta+Nbar*b)-lbeta(alpha,beta)))
  else{
    pvals <- qbeta(c(1:R)/(R+1),alpha,beta)
    return(
      choose(b,ty)*mean((1-(1-pvals)^Nbar)^ty*(1-pvals)^(Nbar*(b-ty)))
    )
  }
}

dty <- function(ty,b,Nbar,alpha,beta,R=1000){
  # returns P[ty] 
  #  where ty|p ~ bin(b,1-(1-p)^Nbar) and p~Beta(alpha,beta)
  # ty, b and Nbar may be vectors
  if(length(b)==1) b <- rep(b,length(ty))
  if(length(Nbar)==1) Nbar <- rep(Nbar,length(ty))
  out <- rep(NA,length(ty))
  for(j in c(1:length(ty))) out[j]<-dty.oneval(ty[j],b[j],Nbar[j],alpha,beta,R=R)
  out
}

negll <- function(par,ty,b,Nbar,freq){
  if(missing(freq)) freq <- rep(1,length(ty))
  alpha <- exp(par[1])
  beta <- exp(par[2])
  -sum(freq*log(dty(ty=ty,b=b,Nbar=Nbar,alpha=alpha,beta=beta)))
}


############################################################
# Fit model by maximum likelihood
############################################################

( o <- optim(c(0,0),negll,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100))
exp(o$par)
exp(o$par)/sum(exp(o$par))

############################################################
# Functions to calculate posterior
############################################################

posterior.Tx <- function(alpha,beta,ty,b,Nbar,B,R=1e5,topcat=10,mu,shape){
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
  # population assumed to be rejected if ty>=1
  #
  # start by finding alpha and beta s.t. P[Xi=0] is as in the unclustered case
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

##################################################################
# get likelihood surface over a grid of values of alpha and beta
##################################################################

parvals <- expand.grid(logalpha=seq(-7,-1,0.01),logitmu=-seq(-13,-4,0.01))
n <- 9400
N <- 1e5
for(k in c(1:nrow(parvals))){
  if(k%%1e3==0) cat("finished row ",k," of ",nrow(parvals)," \n")
  A <- exp(parvals$logalpha[k])
  mu <- 1/(1+exp(parvals$logitmu[k]))
  B <- A*(1-mu)/mu
  parvals$deviance[k] <- 2*(negll(par=c(log(A),log(B)),
                                  ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,
                                  Nbar=100)-o$value)
}

parvals$alpha <- exp(parvals$logalpha)
parvals$mu <- 1/(1+exp(parvals$logitmu))
parvals$beta <- parvals$alpha*(1-parvals$mu)/parvals$mu
parvals$prob.leak <- beta(parvals$alpha,parvals$beta+9400)/beta(parvals$alpha,parvals$beta)-
  beta(parvals$alpha,parvals$beta+1e5)/beta(parvals$alpha,parvals$beta)
parvals$E.leak <- (1e5-9400)*parvals$alpha/(parvals$alpha+parvals$beta+9400)*
  beta(parvals$alpha,parvals$beta+9400)/beta(parvals$alpha,parvals$beta)
parvals$E.leak.given.positive <- parvals$E.leak/parvals$prob.leak

# MLE

mle <- c(alpha=exp(o$par[1]),beta=exp(o$par[2]),
         mu=exp(o$par[1])/sum(exp(o$par)),
         E.leak=(1e5-9400)*exp(o$par[1])/(exp(o$par[1])+exp(o$par[2])+9400)*
           beta(exp(o$par[1]),exp(o$par[2])+9400)/beta(exp(o$par[1]),exp(o$par[2])),
         prob.leak=beta(exp(o$par[1]),exp(o$par[2])+9400)/beta(exp(o$par[1]),exp(o$par[2]))-
           beta(exp(o$par[1]),exp(o$par[2])+1e5)/beta(exp(o$par[1]),exp(o$par[2]))
)
(mle <- c(mle,E.leak.given.positive=mle["E.leak"]/mle["prob.leak"]))

with(parvals[parvals$deviance<=qchisq(0.9,1),],
     rbind(alpha=range(alpha),beta=range(beta),mu=range(mu),
           E.leak=range(E.leak),prob.leak=range(prob.leak),
           E.leak.given.positive=range(E.leak.given.positive)))

with(parvals[parvals$deviance<=qchisq(0.95,1),],
     rbind(alpha=range(alpha),beta=range(beta),mu=range(mu),
           E.leak=range(E.leak),prob.leak=range(prob.leak),
           E.leak.given.positive=range(E.leak.given.positive)))

with(parvals,plot(E.leak.given.positive,deviance,ylim=c(0,20)))
abline(h=qchisq(0.9,df=1),col="blue")
abline(h=qchisq(0.95,df=1),col="red")

tab <- tapply(parvals$deviance,round(parvals$prob.leak,digits=4),FUN=min)
plot(as.numeric(names(tab)),tab)
abline(h=qchisq(0.9,1),col="red")
parvals.sub <- parvals[round(parvals$prob.leak,digits=4)==0.0783,]
parvals.sub[which.min(parvals.sub$deviance),] # alpha=0.0395575, beta=713.3698

############################################################
# Calculate posteriors and produce Figure 1
############################################################

g <- function(prob){
  epsilon <- 0.0001
  prob <- epsilon+(1-2*epsilon)*prob
  log(prob/(1-prob))-log(epsilon/(1-epsilon))
}
probvals <- c(0,0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.15,c(2:4)/5)
probvals <- c(probvals,0.5,1-probvals)

prior1 <- posterior.Tx(alpha=exp(o$par[1]),beta=exp(o$par[2]),ty=0,Nbar=100,B=1000,b=0,R=1e5)
post1 <- posterior.Tx(alpha=exp(o$par[1]),beta=exp(o$par[2]),ty=0,Nbar=100,B=1000,b=94,R=1e5) # MLE

upperCI.probL <- max(parvals$prob.leak[parvals$deviance<=qchisq(0.9,1)])
parvals.upperCI.probL <- parvals[round(parvals$prob.leak,digits=4)==round(upperCI.probL,digits=4),]
parvals.upperCI.probL <- parvals.upperCI.probL[which.min(parvals.upperCI.probL$deviance),]
prior2 <- posterior.Tx(alpha=parvals.upperCI.probL$alpha,beta=parvals.upperCI.probL$beta,
                       ty=0,Nbar=100,B=1000,b=0,R=1e5) # worst case for P(L>0)
post2 <- posterior.Tx(alpha=parvals.upperCI.probL$alpha,beta=parvals.upperCI.probL$beta,
                      ty=0,Nbar=100,B=1000,b=94,R=1e5)

upperCI.EL <- max(parvals$E.leak[parvals$deviance<=qchisq(0.9,1)])
parvals.upperCI.EL <- parvals[round(parvals$E.leak,digits=3)==round(upperCI.EL,digits=3),]
parvals.upperCI.EL <- parvals.upperCI.EL[which.min(parvals.upperCI.EL$deviance),]
prior3 <- posterior.Tx(alpha=parvals.upperCI.EL$alpha,beta=parvals.upperCI.EL$beta,
                       ty=0,Nbar=100,B=1000,b=0,R=1e5) # worst case for P(L>0)
post3 <- posterior.Tx(alpha=parvals.upperCI.EL$alpha,beta=parvals.upperCI.EL$beta,
                      ty=0,Nbar=100,B=1000,b=94,R=1e5)

prior4 <- posterior.Tx(mu=1/9401,shape=Inf,ty=0,Nbar=100,B=1000,b=0,R=1e5) 
    # shape=Inf, mean=1/(n+1)
post4 <- posterior.Tx(mu=1/9401,shape=Inf,ty=0,Nbar=100,B=1000,b=94,R=1e5)

plotpost <- function(prior,post,main,cex.leg=0.8){
  subhead<-paste0("Prior $P\\[T_x>0\\]=",format(round(1-prior[1],digits=3),nsmall=3),
                  "$ ; Posterior $P\\[T_x>0|t_y=0\\]=",
                  format(round(1-post[1],digits=3),nsmall=3),"$")
  a <- barplot(g(post),xaxt="n",yaxt="n",ylim=g(c(0,0.99)),
               xlab=latex2exp::TeX("Number of Infected Seeds in Lot $(T_X)$"),
               ylab="probability",main=latex2exp::TeX(main))
  axis(side=1,at=c(-Inf,a),labels=c("",as.character(0:9),"10+"))
  axis(side=2,at=g(probvals),labels=probvals,las=2,cex.axis=0.8)
  #text(a[3],g(0.975),labels=latex2exp::TeX(subhead),adj=0)
  points(a,g(prior),type="b",lty=1,col=2,pch=1)
  legend("top",legend=c(latex2exp::TeX("Prior $P\\[T_X= ... \\]$"),
                        latex2exp::TeX("Posterior $P\\[T_X= ... | t_y=0 \\]"),
                        " ",
                        latex2exp::TeX(paste0("prior leakage prob $P\\[T_X>0\\]=",
                                              format(round(1-prior[1],digits=3),nsmall=3))),
                        latex2exp::TeX(paste0("posterior leakage prob $P\\[T_X>0|t_y=0\\]=",
                                              format(round(1-post[1],digits=3),nsmall=3)))),
         fill=c(NA,"grey",NA,NA,NA),
         border=c(NA,"black",NA,NA,NA),
         pch=c(1,NA,NA,NA,NA),
         lty=c(1,NA,NA,NA,NA),cex=cex.leg,
         col=c("red",NA,NA,NA,NA))
}

postscript(file="cucurbit_posteriors_26_10_2021.eps")
par(mfrow=c(2,2))
plotpost(prior1,post1,("1: Maximum Likelihood Fit"))
plotpost(prior2,post2,"2: upper 90% CI for $P\\[L>0\\]$")
plotpost(prior3,post3,"3: upper 90% CI for $E\\[L\\]$")
plotpost(prior4,post4,"4: Worst Case for $E\\[L\\]$: $p_i=1/(n+1)$")
dev.off()
par(mfrow=c(1,1))


##################################################################
# Calculate MLEs under models with clustering with various 
# values of theta
##################################################################

dty.clust <- function(ty,b,Nbar,alpha,beta,theta,R=1000){
  # returns P[ty] (where ty must be a single value)
  #  where tyi|pi ~ bin(b,1-(1-pi)^Nbar), pi|p~Beta(theta,theta*(1/p-1)) and
  #     p~Beta(alpha,beta)
  pvals <- qbeta(c(1:R)/(R+1),alpha,beta)
  #phi.vals <- Nbar*pvals - 0.5*Nbar*(Nbar-1)*(pvals^2+pvals*(1-pvals)/(theta+1))
  #phi.vals <- 1 - beta(theta,theta*(1/pvals-1)+Nbar)/beta(theta,theta*(1/pvals-1))
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


negll.clust <- function(par,ty,b,Nbar,freq,theta,R=1000){
  if(missing(freq)) freq <- rep(1,length(ty))
  alpha <- exp(par[1])
  beta <- exp(par[2])
  -sum(freq*log(dty.clust(ty=ty,b=b,Nbar=Nbar,alpha=alpha,beta=beta,theta=theta,R=R)))
}

( o.clust1 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/0.1,R=1e5))
( o.clust2 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/0.5,R=1e5))
( o.clust3 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/1,R=1e5))
( o.clust4 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/5,R=1e5))
( o.clust4b <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/20,R=1e5))
( o.clust5 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/50,R=1e5))
( o.clust6 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/200,R=1e5))
( o.clust7 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/1000,R=1e5))
( o.clust8 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=1/10000,R=1e5)) 
( o.clust9 <- optim(o$par,negll.clust,ty=c(0,4,4,7,21),freq=c(68,1,1,1,1),b=94,Nbar=100,theta=0,R=1e5))

##################################################################
# Calculate posteriors under model with clustering within groups
# and produce Figure 2
##################################################################

set.seed(3010449)
clust1 <- dsn.Tx.clustered(alpha=exp(o.clust1$par[1]),beta=exp(o.clust1$par[2]),theta=1/0.1,Nbar=100,B=1000,b=94,R=1e5)
clust2 <- dsn.Tx.clustered(alpha=exp(o.clust2$par[1]),beta=exp(o.clust2$par[2]),theta=1/0.5,Nbar=100,B=1000,b=94,R=1e5)
clust3 <- dsn.Tx.clustered(alpha=exp(o.clust3$par[1]),beta=exp(o.clust3$par[2]),theta=1/1,Nbar=100,B=1000,b=94,R=1e5)
clust4 <- dsn.Tx.clustered(alpha=exp(o.clust4$par[1]),beta=exp(o.clust4$par[2]),theta=1/5,Nbar=100,B=1000,b=94,R=1e5)
clust4b <- dsn.Tx.clustered(alpha=exp(o.clust4b$par[1]),beta=exp(o.clust4b$par[2]),theta=1/20,Nbar=100,B=1000,b=94,R=1e5)
clust5 <- dsn.Tx.clustered(alpha=exp(o.clust5$par[1]),beta=exp(o.clust5$par[2]),theta=1/50,Nbar=100,B=1000,b=94,R=1e5)
clust6 <- dsn.Tx.clustered(alpha=exp(o.clust6$par[1]),beta=exp(o.clust6$par[2]),theta=1/200,Nbar=100,B=1000,b=94,R=1e5)
clust7 <- dsn.Tx.clustered(alpha=exp(o.clust7$par[1]),beta=exp(o.clust7$par[2]),theta=1/1000,Nbar=100,B=1000,b=94,R=1e5)
clust8 <- dsn.Tx.clustered(alpha=exp(o.clust8$par[1]),beta=exp(o.clust8$par[2]),theta=1/10000,Nbar=100,B=1000,b=94,R=1e5)
clust9 <- dsn.Tx.clustered(alpha=exp(o.clust9$par[1]),beta=exp(o.clust9$par[2]),theta=0,Nbar=100,B=1000,b=94,R=1e5)

set.seed(3010449)
clust1.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/0.1,Nbar=100,B=1000,b=94,R=1e5)
clust2.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/0.5,Nbar=100,B=1000,b=94,R=1e5)
clust3.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/1,Nbar=100,B=1000,b=94,R=1e5)
clust4.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/5,Nbar=100,B=1000,b=94,R=1e5)
clust4b.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/20,Nbar=100,B=1000,b=94,R=1e5)
clust5.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/50,Nbar=100,B=1000,b=94,R=1e5)
clust6.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/200,Nbar=100,B=1000,b=94,R=1e5)
clust7.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/1000,Nbar=100,B=1000,b=94,R=1e5)
clust8.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/10000,Nbar=100,B=1000,b=94,R=1e5)
clust9.naive <- dsn.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=0,Nbar=100,B=1000,b=94,R=1e5)

# SIMPLE SUMMARY CLUSTERED PRIORS and POSTERIORS

set.seed(3010449)
tab1 <- rbind(
  summary.Tx.clustered(alpha=exp(o.clust1$par[1]),beta=exp(o.clust1$par[2]),theta=1/0.1,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust2$par[1]),beta=exp(o.clust2$par[2]),theta=1/0.5,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust3$par[1]),beta=exp(o.clust3$par[2]),theta=1/1,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust4$par[1]),beta=exp(o.clust4$par[2]),theta=1/5,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust4b$par[1]),beta=exp(o.clust4b$par[2]),theta=1/20,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust5$par[1]),beta=exp(o.clust5$par[2]),theta=1/50,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust6$par[1]),beta=exp(o.clust6$par[2]),theta=1/200,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust7$par[1]),beta=exp(o.clust7$par[2]),theta=1/1000,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust8$par[1]),beta=exp(o.clust8$par[2]),theta=1/10000,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o.clust9$par[1]),beta=exp(o.clust9$par[2]),theta=0,Nbar=100,B=1000,b=94,R=1e5))

set.seed(3010449)
tab2 <- rbind(
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/0.1,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/0.5,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/1,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/5,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/20,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/50,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/200,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/1000,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=1/10000,Nbar=100,B=1000,b=94,R=1e5),
  summary.Tx.clustered(alpha=exp(o$par[1]),beta=exp(o$par[2]),theta=0,Nbar=100,B=1000,b=94,R=1e5))

tab1
tab2

postscript(file="cucurbit_clustered_posteriors_26_10_2021.eps")
par(mfrow=c(3,2))
plotpost(clust1$prior,clust1$posterior,("1: $ \\theta^{-1} =0.1$"),cex.leg=0.6)
plotpost(clust2$prior,clust2$posterior,("2: $ \\theta^{-1} =0.5$"),cex.leg=0.6)
plotpost(clust3$prior,clust3$posterior,("3: $ \\theta^{-1} =1$"),cex.leg=0.6)
plotpost(clust4$prior,clust4$posterior,("4: $ \\theta^{-1} =5$"),cex.leg=0.6)
plotpost(clust5$prior,clust5$posterior,("5: $ \\theta^{-1} =50$"),cex.leg=0.6)
plotpost(clust9$prior,clust9$posterior,("6: $ \\theta =0$"),cex.leg=0.6)
dev.off()
par(mfrow=c(1,1))

postscript(file="cucurbit_naive_clustered_posteriors_26_10_2021.eps") # not included in paper
par(mfrow=c(3,2))
plotpost(clust1.naive$prior,clust1.naive$posterior,("1: $ \\theta^{-1} =0.1$"),cex.leg=0.6)
plotpost(clust2.naive$prior,clust2.naive$posterior,("2: $ \\theta^{-1} =0.5$"),cex.leg=0.6)
plotpost(clust3.naive$prior,clust3.naive$posterior,("3: $ \\theta^{-1} =1$"),cex.leg=0.6)
plotpost(clust4.naive$prior,clust4.naive$posterior,("4: $ \\theta^{-1} =5$"),cex.leg=0.6)
plotpost(clust5.naive$prior,clust5.naive$posterior,("5: $ \\theta^{-1} =50$"),cex.leg=0.6)
plotpost(clust9.naive$prior,clust9.naive$posterior,("6: $ \\theta =0$"),cex.leg=0.6)
dev.off()
par(mfrow=c(1,1))




