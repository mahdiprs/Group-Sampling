## Only run examples in interactive R sessions
library(shiny)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(ggplot2)
library(reshape2)
library(actuar)
library(mathjaxr)
library(shinycssloaders)
library(extraDistr)
library(shinythemes)
library(sqldf)
library(kSamples)
library(splus2R)
colsize =3

# functions

fp <- function(p, alpha, beta, Nbar, ty, b) {
    fp =p^(alpha - 1)*(1 - (1 - p)^Nbar)^ty*(1 - p)^(beta - 1 + Nbar*(b - ty))
    return(fp)
}

pf <- function(p, alpha, beta, Nbar, ty, b) {
    #
    #marginal = integrate(function(x) {fp(x,alpha, beta, Nbar, ty, b)}, 0,1)
    pf = fp(p, alpha, beta, Nbar, ty, b)#/marginal$value
    return(pf)
}

Tx_given_ty_NC <-function(N,Nbar,b,ty,l,alpha,beta){
    val =0
    if (l >= ty){
        tmpVal =sapply(0:ty, function(r) (-1)^r *exp(lbeta(l+alpha,N+beta-l)+lchoose(ty,r)+lchoose(N-(b-ty+r)*Nbar,l)))
        tmpVal = sum(tmpVal)
        val = tmpVal
    }else{
        val =0
    }

    #print(cat(val,l))
    return(val)
}


Tx_tx_given_ty_NC <-function(N,Nbar,b,ty,l,alpha,beta){
    val =0
        tmpVal =sapply(0:ty, function(r) (-1)^r *exp(lbeta(l+alpha,N+Nbar*r-Nbar*ty+beta-l)+lchoose(ty,r)+lchoose(N-b*Nbar,l)))
        tmpVal = sum(tmpVal)
        val = tmpVal
    return(val)
}


# clustered functions
phi =function(p,theta,Nbar){
    phi =  1 - gamma(theta*(1/p - 1) + Nbar)*gamma(theta/p)/(gamma(theta/p + Nbar)*gamma(theta/p - theta))
    return(phi)
}

fpc = function(p, alpha, beta,theta, Nbar, ty, b){
    fpc =p^(alpha - 1)*(1 - p)^(beta - 1)*phi(theta, p, Nbar)^ty*(1 - phi(theta, p, Nbar))^(b - ty)
    return(fpc)
}

pfc <- function(p, alpha, beta, theta, Nbar, ty, b) {
    #
    #marginal = integrate(function(x) {fpc(x,alpha, beta, theta,Nbar, ty, b)}, 0,1)
    pfc = fpc(p, alpha, beta, theta, Nbar, ty, b)#/marginal$value
    return(pfc)
}

# zero truncated beta binomial distribution pdf
ztbb <- function(k,Nbar,alpha,beta) {
    #
    if (k>0){
        ztbb = dbbinom(k, Nbar, alpha = alpha, beta = beta, log = FALSE)/
            (1-dbbinom(0, Nbar, alpha = alpha, beta = beta, log = FALSE))
    }else{
        ztbb = 0
    }
    return(ztbb)
}

# simple sampling algorithm
sample_fun <- function(n,support, probs) sample(
    support, # the support of pdf
    n, # n draws
    TRUE,  # with replacement
    probs # using these probabilities
)

# adjust cell
widthCell =20
widthbar =270
steps = 6
colsize =3
maxgrowth =3
# annotate setting
hval =-1.5
# number of simulation for Tx
nsimVal = 100000

textInputRow<-function (inputId, label, value = "")
{
    div(style="display:inline-block",
        tags$label(label, `for` = inputId),
        tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

# Server code used for all examples
    server <- function(input, output,session) {
        #p = .1
        ########################################################################
        ## non clustered
        #######################################################################
        data <- eventReactive(input$go, {
            # number of units in the population
            N= input$consignment
            # number of units in the sample tested
            n = input$tested_n# .2 * N
            # beta dist parameters
            alpha = input$alpha
            beta =input$beta
            # number of units in in a single bulk sample
            Nbar = input$Nbar
            ty = input$ty
            # number of bulk samples tested bNbar = n
            b=floor(n/Nbar)
            # Number of bulk samples in the population
            B= floor(N/Nbar)

            # input validation
            validate(
                need(input$consignment>0 & input$tested_n >0, '$N$ and $n$ must be greater than 0'),
                need(input$tested_n >0, '$n$ must be greater than 0'),
                need(input$consignment>input$tested_n, '$N$ must be greater than $n$'),
                need(input$alpha >0 & input$beta >0, '$\\alpha$ and $\\beta$ must be greater than 0'),
                need(input$ty >=0 & ty%%1==0 , '$t_y$ must be an integer greater and equal to 0'),
                need(input$ty <=b , paste0('$t_y$ must be smaller than ',b+1)),
                need(input$Nbar >0 & input$Nbar < input$tested_n , '$\\bar{N}$ must be greater than 0 and less than $n$')
            )
            # number of sample nsim
            nsim = nsimVal
            prior = rbeta(nsim, alpha, beta)
            prior = as.data.frame(prior)

            # posterior
            support <- seq(0.0, 1, 0.0001)
            post_size = length(support)
            #post_size
            posterior <- data.frame(xp = rep(NA, post_size), yp = rep(NA, post_size))
            for (i in c(1:post_size)) {
                posterior[i, 1] <- fp(support[i], alpha, beta, Nbar, ty, b)
                posterior[i, 2] <- support[i]
            }

            posterior = subset(posterior, !is.na(xp) &!is.inf(xp))
            posterior$xp = posterior$xp/sum(posterior$xp)
            Tx = data.frame(x= rep(NA,N),y= rep(NA,N))

            if (ty>0){
                marginal = integrate(function(x) {fp(x,alpha, beta, Nbar, ty, b)}, 0,1)
                for (i in c(1:N)) {
                    Tx[i,1]= Tx_given_ty_NC(N,Nbar,b,ty,i,alpha,beta)/marginal$value
                    Tx[i,2] =i
                }
            }else{
                i<-0
                tmpCnt =0
                marginal = lbeta(alpha, Nbar*b + beta)
                while (i < N) {
                    Tx[i+1,1]= exp(lbeta(alpha + i, N + beta - i) -marginal+lchoose(N-b*Nbar,i))
                    Tx[i+1,2] =i
                    tmpCnt = tmpCnt+Tx[i+1,1]
                    if (tmpCnt>0.9999){
                        i <- N+1
                    } else{
                        i =i+1
                    }
                }
            }
            Tx = subset(Tx, !is.na(x))

            Tx_tx = data.frame(xtx= rep(NA,N),ytx= rep(NA,N))

            if (ty>0){
                marginal = integrate(function(x) {fp(x,alpha, beta, Nbar, ty, b)}, 0,1)
                for (i in c(1:N)) {
                    Tx_tx[i,1]= Tx_tx_given_ty_NC(N,Nbar,b,ty,i,alpha,beta)/marginal$value
                    Tx_tx[i,2] =i
                }
            }else{
                Tx_tx = Tx
                colnames(Tx_tx) <- c('xtx','ytx')
            }

            Tx_tx = subset(Tx_tx, !is.na(xtx))
            data = bind_rows(prior,posterior,Tx,Tx_tx)
            return(data)
        })


        output$plotPriorNC <- renderPlot({
            tmp = data()
            prior = subset(tmp, !is.na(prior))
            # mean
            alpha = input$alpha
            beta =input$beta
            mean = alpha/(alpha+beta)#mean(prior$prior,na.rm = T)
            median = median(prior$prior, na.rm = T)
            sd = sd(prior$prior, na.rm = T)
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)
            if (mean >median){
                trs = mean+3*sd
            }else{
                trs = median+3*sd
            }

            #       plot beta prior
            ggplot(data=subset(prior, prior<trs), aes(x=prior,y=..count../sum(..count..)))+
                geom_histogram(alpha=.5, position="identity")+
                labs(x=paste("prevalence (p)\n", "Mean= ", round(mean,digits = 8),", ", "Median= ", round(median,digits = 8),", ", "sd= ",
                             round(sd,digits = 8),collapse = ", "), y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Prior: (p)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())
        })
        # posterior prevalence
        output$plotPostPreNC <- renderPlot({
            tmp = data()
            posterior = subset(tmp, !is.na(xp))

            mean = sum(posterior$xp*posterior$yp)

            sd =sqrt(sum(posterior$xp* (posterior$yp-mean)^2))

            i<-1
            tmpMed =0
            while (i < nrow(posterior)) {
                tmpMed = tmpMed+posterior$xp[i]
                if (tmpMed>0.5){
                    median = posterior$yp[i]
                    i <- nrow(posterior)+1
                } else{
                    i =i+1
                }
            }
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)
            if (mean >median){
                trs = mean+3*sd
            }else{
                trs = median+3*sd
            }
            ggplot(data=subset(posterior, yp<trs), aes(x=yp,y=xp))+
                geom_col(fill = "grey")+
                labs(x=paste("prevalence (p|ty)\n", "Mean= ", round(mean,digits = 6), ", ",
                             "Median= ",round(median,digits = 6),", ",
                             "sd= ",round(sd,digits = 6)),y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Posterior: (p|ty)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())
        })

        # posterior Tx
        output$plotPostTx <- renderPlot({
            #sample = simulated
            # plot Tx histogram
            tmp = data()
            Tx = subset(tmp, !is.na(x))
            Tx$x = Tx$x/sum(Tx$x)
            #Tx$y = 1:nrow(Tx)
            mean = sum(Tx$x*Tx$y)
            #median = NA #median(Tx$x, na.rm = T)
            sd =sqrt(sum(Tx$x* (Tx$y-mean)^2))

            i<-1
            tmpMed =0
            N = input$consignment
            while (i < N) {
                tmpMed = tmpMed+Tx$x[i]
                if (tmpMed>=0.5){
                    median = Tx$y[i]
                    i <- N+1
                } else{
                    i =i+1
                }
            }

            # summary stats
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)

            if (mean >median){
                trs = mean+3*sd
            }else{
                trs = median+3*sd
            }

            ggplot(data=subset(Tx, y< trs), aes(x=y,y=x))+
                geom_col(fill = "grey")+
                labs(x=paste("Nbr contaminated units (Tx|ty)\n", "Mean= ", round(mean,digits = 6), ", ",
                             "Median= ",round(median,digits = 6),", ",
                             "sd= ",round(sd,digits = 6)),y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Posterior:  (Tx|ty)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())

        })

        # posterior Tx_tx
        output$plotPostTx_tx <- renderPlot({
            # plot Tx histogram
            tmp = data()
            Tx_tx = subset(tmp, !is.na(xtx))
            Tx_tx$xtx = Tx_tx$xtx/sum(Tx_tx$xtx)
            mean = sum(Tx_tx$xtx*Tx_tx$ytx)
            sd =sqrt(sum(Tx_tx$xtx* (Tx_tx$ytx-mean)^2))

            i<-1
            tmpMed =0
            N = input$consignment
            while (i < N) {
                tmpMed = tmpMed+Tx_tx$xtx[i]
                if (tmpMed>0.5){
                    median = Tx_tx$ytx[i]
                    i <- N+1
                } else{
                    i =i+1
                }
            }

            # summary stats
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)

            if (mean >median){
                trs = mean+3*sd
            }else{
                trs = median+3*sd
            }

            ggplot(data=subset(Tx_tx, ytx< trs), aes(x=ytx,y=xtx))+
                geom_col(fill = "grey")+
                labs(x=paste("Nbr contaminated units (Tx-tx|ty)\n", "Mean= ", round(mean,digits = 6), ", ",
                             "Median= ",round(median,digits = 6),", ",
                             "sd= ",round(sd,digits = 6)),y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Posterior: (Tx-tx|ty)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())

        })

        ########################################################################
        ## with clustering
        #######################################################################
        dataC <- eventReactive(input$goC, {
            # number of units in the population
            N= input$consignmentC
            # number of units in the sample tested
            n = input$tested_nC# .2 * N
            # beta dist parameters
            alpha = input$alphaC
            beta =input$betaC
            # Clustering level
            theta = input$theta
            # number of units in in a single bulk sample
            Nbar = input$NbarC
            ty = input$tyC
            # number of bulk samples tested bNbar = n
            b=floor(n/Nbar)
            # Number of bulk samples in the population
            B= floor(N/Nbar)
            # sample observations of yi, i=1..b
            #p = .1

            # input validation
            validate(
                need(input$consignmentC>0 & input$tested_nC >0, '$N$ and $n$ must be greater than 0'),
                need(input$tested_nC >0, '$n$ must be greater than 0'),
                need(input$consignmentC>input$tested_nC, '$N$ must be greater than $n$'),
                need(input$alphaC >0 & input$betaC >0, '$\\alpha$ and $\\beta$ must be greater than 0'),
                need(input$tyC >=0 & input$tyC%%1==0 , '$t_y$ must be an integer greater and equal to 0'),
                need(input$tyC <=b , paste0('$t_y$ must be smaller than ',b+1)),
                need(input$theta >0 , '$\\thetha$ must be greater than 0'),
                need(input$NbarC >0 & input$NbarC < input$tested_nC , '$\\bar{N}$ must be greater than 0 and less than $n$')
            )
            # set support of pdf
            support <- seq(0.000001, 1, 0.000001)

            # Find the probability densities over the support
            probs <- pfc(support, alpha, beta, theta, Nbar, ty, b)
            probs<-probs[!is.na(probs)]
            probs <- probs/sum(probs) # normalize

            # # number of sample nsim
            nsim = 10000
            simulated <- sample_fun(nsim,support,probs)
            posterior = simulated
            posterior = as.data.frame(posterior)
            prior = rbeta(nsim, alpha, beta)
            prior = as.data.frame(prior)
            y = c(rep(1, ty), rep(0, b-ty), rep(1, B-b))
            Txc = data.frame(x= rep(0,length(nsim))) # Tx in the paper
            txcl = data.frame(tx= rep(0,length(nsim))) # tx in the paper

            supportZTBB <- seq(0, Nbar, 1)

            # Find the probability densities over the support
            for (cnt in c(1:nsim)){
                p = simulated[cnt]
                Xis <- 0
                if (ty>0){
                    probsTBB = data.frame(x= rep(NA,Nbar),y= rep(NA,Nbar))
                    for (i in c(0:Nbar)){
                        probsTBB[i+1,1] = ztbb(i, Nbar,theta,theta*(1/p-1))
                        probsTBB[i+1,2] = i
                    }
                    probsTBB<-subset(probsTBB, !is.na(x))
                    probsTBB <- probsTBB$x/sum(probsTBB$x) # normalize

                    tmpVal = sum(sample_fun(ty,supportZTBB,probsTBB))
                    Txc[cnt,1] = sum(rbbinom(B-b, Nbar,theta,theta*(1/p-1)))+
                        tmpVal
                    txcl[cnt,1] = tmpVal
                }else{
                    Txc[cnt,1] = sum(rbbinom(B-b, Nbar,theta,theta*(1/p-1)))
                    txcl[cnt,1] = 0
                }
            }


            data = bind_rows(prior,posterior,Txc,txcl)
            return(data)
        })

        output$plotPriorC <- renderPlot({
            tmp = dataC()
            prior = subset(tmp, !is.na(prior))
            mean = mean(prior$prior,na.rm = T)
            median = median(prior$prior, na.rm = T)
            sd = sd(prior$prior, na.rm = T)
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)
            ##  plot beta prior
            ggplot(data=subset(prior, prior<= max(mean,median)+3*sd), aes(x=prior,y=..count../sum(..count..)))+
                geom_histogram(alpha=.5, position="identity",bins =  nclass.FD(prior$prior))+
                labs(x=paste("prevalence (p)\n", "Mean= ", round(mean,digits = 6),", ", "Median= ",round(median,digits = 6),", ",
                             "sd= ",round(sd,digits = 6)),
                     y="pdf", show.legend = TRUE)+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Prior: (p)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())
        })
        # posterior prevalence
        output$plotPostPreC <- renderPlot({
            tmp = dataC()
            posterior = subset(tmp, !is.na(posterior))
            mean = mean(posterior$posterior,na.rm = T)
            median = median(posterior$posterior, na.rm = T)
            sd = sd(posterior$posterior, na.rm = T)
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)
            #       plot beta prior
            ggplot(data=subset(posterior,posterior<=max(mean,median)+3*sd), aes(x=posterior,y=..count../sum(..count..)))+
                geom_histogram(alpha=.5, position="identity",bins =  nclass.FD(posterior$posterior))+
                labs(x=paste("prevalence (p|ty)\n", "Mean= ", round(mean,digits = 6),", ", "Median= ",round(median,digits = 6),", ",
                             "sd= ",round(sd,digits = 6)), y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Posterior: (p|ty)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())
        })

        # posterior Tx
        output$plotPostTxC <- renderPlot({
            tmp = dataC()
            Tx = subset(tmp, !is.na(x))
            mean = mean(Tx$x,na.rm = T)
            median = median(Tx$x, na.rm = T)
            sd = sd(Tx$x, na.rm = T)
            stats <- data.frame(Ref = c("Mean", "Median"),
                                vals = c(mean, median),
                                stringsAsFactors = FALSE)

            ggplot(data=subset(Tx,x<= max(median,mean)+3*sd), aes(x=x,y=..count../sum(..count..)))+
                geom_histogram(binwidth = 1,alpha=.5, position="identity",bins =  nclass.FD(Tx$x))+
                labs(x=paste("Nbr contaminated units (Tx|ty)\n", "Mean= ", round(mean,digits = 6),", ",
                             "Median= ",round(median,digits = 6), ", ", "sd= ",round(sd,digits = 6)), y="pdf")+
                geom_vline(aes(xintercept=median,
                               color="median"), linetype="solid",
                           size=1) +
                geom_vline(aes(xintercept=mean,
                               color="mean"), linetype="solid",
                           size=1) +
                scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                ggtitle("Posterior: (Tx|ty)") +
                theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                    legend.justification=c(1,0),
                                                                    legend.position=c(.8, .75),
                                                                    legend.background = element_blank(),
                                                                    legend.key = element_blank())
        })

        # posterior Tx-tx
        output$plotPostTx_txC <- renderPlot({

            tmp = dataC()
            Txc = subset(tmp, !is.na(x))
            tycc = input$tyC
            if (tycc ==0){
                mean = mean(Txc$x,na.rm = T)
                median = median(Txc$x, na.rm = T)
                sd = sd(Txc$x, na.rm = T)
                stats <- data.frame(Ref = c("Mean", "Median"),
                                    vals = c(mean, median),
                                    stringsAsFactors = FALSE)

                ggplot(data=subset(Txc, x<=max(mean,median)+3*sd), aes(x=x,y=..count../sum(..count..)))+
                    geom_histogram(binwidth = 1,alpha=.5, position="identity",bins =  nclass.FD(Txc$x))+
                    labs(x=paste("Nbr contaminated units (Tx-tx|ty)\n", "Mean= ", round(mean,digits = 6),", ",
                                 "Median= ",round(median,digits = 6), ", ", "sd= ",round(sd,digits = 6)), y="pdf")+
                    geom_vline(aes(xintercept=median,
                                   color="median"), linetype="solid",
                               size=1) +
                    geom_vline(aes(xintercept=mean,
                                   color="mean"), linetype="solid",
                               size=1) +
                    scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                    ggtitle("Posterior: (Tx-tx|ty)") +
                    theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                        legend.justification=c(1,0),
                                                                        legend.position=c(.8, .75),
                                                                        legend.background = element_blank(),
                                                                        legend.key = element_blank())
            }else{
                txcl =subset(tmp, !is.na(tx))
                tt = sqldf("select *, count(x) as feq from Txc group by
           x")
                tt$feq = tt$feq/sum(tt$feq)

                tt1 = sqldf("select *, count(tx) as feq1 from txcl group by
           tx")
                tt1$feq1 = tt1$feq1/sum(tt1$feq1)

                x1 <- tt$x
                p1 <- tt$feq
                x2 <- -tt1$tx
                p2 <- tt1$feq1

                Tx_txc = conv(x1,p1,x2,p2)
                Tx_txc = as.data.frame(Tx_txc)
                Tx_txc = subset(Tx_txc,V1>0)
                Tx_txc$V2 = Tx_txc$V2/sum(Tx_txc$V2)

                mean = sum(Tx_txc$V1*Tx_txc$V2)
                sd =sqrt(sum(Tx_txc$V2* (Tx_txc$V1-mean)^2))

                i<-1
                tmpMed =0
                N = input$consignmentC
                while (i < N) {
                    tmpMed = tmpMed+Tx_txc$V2[i]
                    if (tmpMed>0.5){
                        median = Tx_txc$V1[i]
                        i <- N+1
                    } else{
                        i =i+1
                    }
                }

                # summary stats
                stats <- data.frame(Ref = c("Mean", "Median"),
                                    vals = c(mean, median),
                                    stringsAsFactors = FALSE)

                if (mean >median){
                    trs = mean+3*sd
                }else{
                    trs = median+3*sd
                }

                ggplot(data=subset(Tx_txc, V1< trs), aes(x=V1,y=V2))+
                    geom_col(fill = "grey")+
                    labs(x=paste("Nbr contaminated units (Tx-tx|ty)\n", "Mean= ", round(mean,digits = 6), ", ",
                                 "Median= ",round(median,digits = 6),", ",
                                 "sd= ",round(sd,digits = 6)),y="pdf")+
                    geom_vline(aes(xintercept=median,
                                   color="median"), linetype="solid",
                               size=1) +
                    geom_vline(aes(xintercept=mean,
                                   color="mean"), linetype="solid",
                               size=1) +
                    scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))+
                    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                    ggtitle("Posterior: (Tx-tx|ty)") +
                    theme(plot.title = element_text(hjust = 0.5))+theme(legend.title = element_text(size=12, color = "salmon", face="bold"),
                                                                        legend.justification=c(1,0),
                                                                        legend.position=c(.8, .75),
                                                                        legend.background = element_blank(),
                                                                        legend.key = element_blank())
            }


        })
    }

    # All cells at 300 pixels wide, with cell padding
    # and a border around everything
    ui <- navbarPage(theme = shinytheme("superhero"),
                     "Group sampling inference tool",
        tabPanel("Without clustering",
        fluidRow(withMathJax(),
          tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({tex2jax: {
                    inlineMath: [ ['$','$'], ['\\(','\\)'] ]}});
                </script>")),
        div(style = "display:inline-block"),responsive = TRUE, br(),
        splitLayout(
        style = "border: 1px solid silver;",
        cellWidths = c("25%", "25%","25%","25%"),
        cellArgs = list(style = "padding: 6px"),
        withSpinner(plotOutput("plotPriorNC", height = 400, width = 400),type=4),
        withSpinner(plotOutput("plotPostPreNC", height = 400, width = 400),type=4),
        withSpinner(plotOutput("plotPostTx", height = 400, width = 400),type=4),
        withSpinner(plotOutput("plotPostTx_tx", height = 400, width = 400),type=4)
            ),splitLayout(
                style = "border: 1px solid silver;",
                cellWidths = c("64%","36%"),
                cellArgs = list(style = "padding: 12px"),
        column(colsize, h5("Model inputs"),
               fluidRow(
                   column(6,
                          style = list("padding-right: 5px;"),
                          numericInput("consignment", label = "$N$", value = 10000, step = 500)
                   ),
                   column(6,
                          style = list("padding-left: 5px;"),
                          numericInput("tested_n", label = "$n$", value = 2000, step = 0.5)
                   )
               ),
               fluidRow(
                   column(6,
                          style = list("padding-right: 5px;"),
                          numericInput("alpha", label = "$\\alpha$", value = 1, step = 0.05)
                   ),
                   column(6,
                          style = list("padding-left: 5px;"),
                          numericInput("beta", label = "$\\beta$", value = 100, step = 0.01)
                   )
               ),
               fluidRow(
                   column(6,
                          style = list("padding-right: 5px;"),
                          numericInput("Nbar", label = "$\\bar{N}$", value = 100, step = 0.5)
                   ),
                   column(6,
                          style = list("padding-left: 5px;"),
                          numericInput("ty", label = "$t_y$", value = 3, step = .5)
                   ), fluidRow(
                       actionButton("go", "Update"), br(),br(),
                       fluidRow(
                           h5("Set the parameters then ‘’Update’’ to generate the plots"))
                   ))
                  ),
        column(colsize, fluidRow(column(12,
                                h3("Model parameter description"),
                                "$N$: number of units in the population", br(), "$n$: number of units in the sample tested", br(),
                                "$\\alpha$: Shape parameter of beta distribution $($prior$)$",br(),
                                "$\\beta$: Scale parameter of beta distribution $($prior$)$",br(),
                                "$\\bar{N}$: number of units in a single bulk-sample $($group$)$", br(),
                                "$t_y$: number of infected bulk-samples $($groups$)$ in the sample tested", br(),
                                "$t_x$: number of infected units in the sample tested", br(),
                                "$T_x$: number of infected units in the population", br(),
                                "$p$: probability that a single drawn unit is infected", br()
        ))))
        )),tabPanel("With clustering",
                    fluidRow(withMathJax(),
                             tags$div(HTML("<script type='text/x-mathjax-config'>
                MathJax.Hub.Config({tex2jax: {
                    inlineMath: [ ['$','$'], ['\\(','\\)'] ]}});
                </script>")),
                             div(style = "display:inline-block"),responsive = TRUE, br(),
                             splitLayout(
                                 style = "border: 1px solid silver;",
                                 cellWidths = c("25%", "25%","25%","25%"),
                                 cellArgs = list(style = "padding: 6px"),
                                 withSpinner(plotOutput("plotPriorC", height = 400, width = 400),type=4),
                                 withSpinner(plotOutput("plotPostPreC", height = 400, width = 400),type=4),
                                 withSpinner(plotOutput("plotPostTxC", height = 400, width = 400),type=4),
                                 withSpinner(plotOutput("plotPostTx_txC", height = 400, width = 400),type=4)
                             ),splitLayout(
                                 style = "border: 1px solid silver;",
                                 cellWidths = c("64%","36%"),
                                 cellArgs = list(style = "padding: 6px"),
                                 column(colsize, h5("Model inputs"),
                                        fluidRow(
                                            column(6,
                                                   style = list("padding-right: 5px;"),
                                                   numericInput("consignmentC", label = "$N$", value = 10000, step = 500)
                                            ),
                                            column(6,
                                                   style = list("padding-left: 5px;"),
                                                   numericInput("tested_nC", label = "$n$", value = 2000, step = 0.5)
                                            )
                                        ),
                                        fluidRow(
                                            column(6,
                                                   style = list("padding-right: 5px;"),
                                                   numericInput("alphaC", label = "$\\alpha$", value = 1, step = 0.05)
                                            ),
                                            column(6,
                                                   style = list("padding-left: 5px;"),
                                                   numericInput("betaC", label = "$\\beta$", value = 100, step = 0.01)
                                            )
                                        ),
                                        fluidRow(
                                            column(6,
                                                   style = list("padding-right: 5px;"),
                                                   numericInput("NbarC", label = "$\\bar{N}$", value = 100, step = 0.5)
                                            ),
                                            column(6,
                                                   style = list("padding-left: 5px;"),
                                                   numericInput("tyC", label = "$t_y$", value = 3, step = .5)
                                            ),column(6,
                                                     style = list("padding-right: 5px;"),
                                                     numericInput("theta", label = "$\\theta$", value = 3, step = .5), br(),
                                            ),
                                            fluidRow(
                                                actionButton("goC", "Update"), br(),br(),
                                                fluidRow(
                                                    h5("Set the parameters then ‘’Update’’ to generate the plots"))
                                            ))
                                 ),
                                 column(colsize, fluidRow(column(12,
                                                                 h3("Model parameter description"),
                                                                 "$N$: number of units in the population", br(), "$n$: number of units in the sample tested", br(),
                                                                 "$\\alpha$: Shape parameter of beta distribution $($prior$)$",br(),
                                                                 "$\\beta$: Scale parameter of beta distribution $($prior$)$",br(),
                                                                 "$\\bar{N}$: number of units in a single bulk-sample $($group$)$", br(),
                                                                 "$t_y$: number of infected bulk-samples $($groups$)$ in the sample tested", br(),
                                                                 "$t_x$: number of infected units in the sample tested", br(),
                                                                 "$T_x$: number of infected units in the population", br(),
                                                                 "$p$: probability that a single drawn unit is infected", br(),
                                                                 "$\\theta$: level of clustering"
                                 ))))
                    ))
        )
    shinyApp(ui, server)

