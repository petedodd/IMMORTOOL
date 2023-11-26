## the non-shiny functions

##' Function to simulate a cohort
##'
##' This function will simulate a cohort with relevant times-to-event under a null assumption for the treatment effect. Typical use returns summary statistics.
##' 
##' @title ITBstats
##' @param N The size of the cohort for simulation (unrelated to the experimental study cohort size)
##' @param rtt.exposure random time-to-exposure generator
##' @param rtt.death random time-to-death generator
##' @param rtt.ltfu random time-to-ltfu generator
##' @param returnraw Logical (default FALSE) controlling return type. If true, the raw cohort will be returned. Otherwise, a list of summary statistics will be returned
##' @return Depends on [returnraw]
##' @author Pete Dodd
##' @export
##' @import data.table
##' @import gmodels
ITBstats <- function(N=1e4,           #simulation cohort size
                     Tstop=90,        #end time TODO include below w/slider
                     rtt.exposure,    #random time-to-exposure
                     rtt.death,       #random time-to-death
                     rtt.ltfu,        #random time-to-ltfu
                     returnraw=FALSE, #return TZ
                     ...
                     ){
  ## TODO set seed
  ## create cohort
  TZ <- data.table::data.table(t.e=rtt.exposure(N),
                   t.d=rtt.death(N),
                   t.l=rtt.ltfu(N))

  TZ[,exposed:=ifelse(t.e<pmin(t.d,t.l,Tstop),TRUE,FALSE)]
  TZ[,died:=ifelse(t.d<pmin(t.l,Tstop),TRUE,FALSE)]
  TZ[1:floor(N/2),exposed:=FALSE] #ensure some non-exposed

  ## print?
  gmodels::CrossTable(TZ$exposed,TZ$died)
  ## use person time to death by group to calculate death rate

  ## a) person time from 0
  TZ[,PTa:=pmin(t.d,t.l,Tstop)]
  NPT <- TZ[,.(deaths=sum(died),PT=sum(PTa)),by=exposed]
  NPT[,rate:=deaths/PT]
  (RR2 <- NPT[exposed==TRUE]$rate / NPT[exposed==FALSE]$rate)

  suissa.p <- TZ[exposed==TRUE,sum(t.e)/sum(PTa)] #Suissa p fraction of T1 that us U

  ## b) person time from exposure
  TZ[,PTb:=ifelse(exposed==TRUE,PTa-t.e,PTa)]
  NPT2 <- TZ[,.(deaths=sum(died),PT=sum(PTb)),by=exposed]
  NPT2[,rate:=deaths/PT]
  (RR3 <- NPT2[exposed==TRUE]$rate / NPT2[exposed==FALSE]$rate)

  ## return
  if(!returnraw){
    list(table.a=NPT,RR.a=RR2,table.b=NPT2,RR.b=RR3,
         F.e=TZ[,mean(exposed)],F.d=TZ[,mean(died)], #NOTE ltfu=1-death
         suissa.k=NPT[exposed==TRUE]$PT / NPT[exposed==TRUE]$PT,
         suissa.p=suissa.p,
         SElnIRR1.a = sqrt(N*sum(1/NPT$deaths)),
         SElnIRR1.b = sqrt(N*sum(1/NPT2$deaths))
         )
  } else {
    return(TZ)
  }
}

##' Function to visualize input distributions
##'
##' This plots the exposure (treatment), mortality, and LTFU time-to-event distributions
##'
##' @title makeDistPlot
##' @param inputs a list of parameters containing the Weibull parameters for exposure (k.e, L.e), death (k.d, L.d), and ltfu (k.l, L.l)
##' @param k.e (tmp)
##' @param k.d  (tmp)
##' @param k.l  (tmp)
##' @param N The size of the cohort for simulation (unrelated to the experimental study cohort size)
##' @param rtt.exposure random time-to-exposure generator
##' @param rtt.death random time-to-death generator
##' @param rtt.ltfu random time-to-ltfu generator
##' @param returnraw Logical (default FALSE) controlling return type. If true, the raw cohort will be returned. Otherwise, a list of summary statistics will be returned
##' @return A ggplot2 graph
##' @author Pete Dodd
##' @export
##' @import ggplot2
makeDistPlot <- function(input){
  ggplot2::ggplot() +
    ggplot2::xlim(0,input$T.max) +
    ggplot2::geom_function(ggplot2::aes(colour="exposure"),fun=dweibull,
                           args=list(shape=input$k.e,scale=input$L.e)) +
    ggplot2::geom_function(ggplot2::aes(colour="death"),fun=dweibull,
                           args=list(shape=input$k.d,scale=input$L.d)) +
    ggplot2::geom_function(ggplot2::aes(colour="LTBU"),fun=dweibull,
                           args=list(shape=input$k.l,scale=input$L.l)) +
    ggplot2::xlab('Time') + ggplot2::ylab('') +
    ggplot2::theme(legend.title=ggplot2::element_blank(),legend.position='top')
}


##' Make plot of effect and CIs by sample size
##'
##' This is the effect by various analyses, with approximate confidence intervals, plotted against sample size.
##'
##' @title make Effect by Sample size Plot
##' @param ANS A list of summary statistics returned by [ITBstats]
##' @return A ggplot2 graph
##' @author Pete Dodd
##' @export
##' @import data.table
##' @import ggplot2
makeEffectBySamplePlot <- function(ANS){
  ## make plot data
  pfd <- data.table::data.table(n=c(10,20,30,50,80,100,2:10*100))
  pfd[,RR.a:=ANS$RR.a]; pfd[,RR.b:=ANS$RR.b]
  pfd <- data.table::melt(pfd,id='n')
  pfd[,c('qnt','type'):=data.table::tstrsplit(variable,split='\\.')]
  pfd[,SElnIRR1:=ANS$SElnIRR1.a] #same everywhere
  pfd[,RRhi:=exp(log(value) + 1.96*SElnIRR1/sqrt(n))]
  pfd[,RRlo:=exp(log(value) - 1.96*SElnIRR1/sqrt(n))]
  pfd[,analysis:=type]
  ggplot2::ggplot(pfd,ggplot2::aes(x=n,y=value,ymin=RRlo,ymax=RRhi,col=analysis,group=analysis))+
    ggplot2::geom_line()+## geom_pointrange()+
    ggplot2::geom_line(ggplot2::aes(y=RRhi),lty=2)+
    ggplot2::geom_line(ggplot2::aes(y=RRlo),lty=2)+
    ggplot2::xlab('Sample size')+ggplot2::ylab('Rate ratio')+
    ggplot2::geom_hline(yintercept=1,col='black',lty=3)+
    ggplot2::theme(legend.position='top')
}

## Weibull scaled to have max height = 1
dweibull1 <- function(x,shape,scale){
  kk <- (shape-1)/shape
  dweibull(x=x,shape=shape,scale=scale)/
    ((shape/scale) * (kk^kk) * exp(-kk) )
}

##' Plot survival and treatment distributions.
##'
##' The plot returned is the survival (from death) function with a scaled time-to-treatment distribution superimposed. Note: the time-to-treatment distribution will differ from the time-to-treatment distribution in the data due to competition with mortality. The competition with mortality is accounted for in fitting and simulation, but not in this plot.
##' @title Plot survival and treatment distributions.
##' @param input a list of parameters including Weibull treatment distribution parameters (k.e, L.e), Weibull mortality distribution parameters (k.d, L.d), and the maximum time (T.max).
##' @return A plot with survival and treatment times
##' @author Pete Dodd
##' @export
makeTMplot <- function(input){
  ggplot2::ggplot() +
    ggplot2::xlim(0,input$T.max) +
    ggplot2::geom_function(ggplot2::aes(colour="exposure times"),fun=dweibull1,
                           args=list(shape=input$k.e,scale=input$L.e),n=1e3) +
    ggplot2::geom_function(ggplot2::aes(colour="survival function"),
                           fun=function(x) exp(-(x/input$L.d)^input$k.d),n=1e3) +
    ggplot2::xlab('Time') + ggplot2::ylab('') +
    ggplot2::theme(legend.title=ggplot2::element_blank(),legend.position='top')
}



## --- some smaller non-exported utilities ----

## == now for death times:

## ' Errors in cumulative mortality
## '
## ' This function is a utility for fitting a Weibull to the observed mortality in a control arm.
## '
## ' @title Residual for observed mortality at time T
## ' @param T A 2xN matrix: the first column being observation times; the second column being fractional mortality at each time in the (control) cohort
## ' @param km A proposed Weibull shape parameger
## ' @param Lm A proposed Weibull scale parameter
## ' @return The error as a difference
## ' @author Pete Dodd
CFRET <- function(T,km,Lm) 1-exp(-(T[1]/Lm)^km) - T[2] #first arg time, second corresponding mort'y

## ' SSE mortality errors
## '
## ' This is a utility function for fitting a Weibull to observed mortality in a control arm..
## '
## ' @title SSE for fitting a Weibull to observed mortality
## ' @param M A 2xN matrix: the first column being observation times; the second column being fractional mortality at each time in the (control) cohort
## ' @param x A 2-vector equal to (log(Weibull shape), log(Weibull scale))
## ' @return The SSE for the Weibull and these data
## ' @author Pete Dodd
morterr <- function(M,x){
  x <- exp(x)
  tmp <- apply(M,1,function(y) (CFRET(y,x[1],x[2]))^2) #vector of squared errors
  sum(tmp)
}


##' For fitting to mortality data.
##'
##' This takes data from a control arm on cumulative mortality and fits a Weibull distribution to it so as to minimize the sum-of-squares error.
##'
##' @title Get the best-fit parameters for mortality
##' @param M A 2xN matrix: the first column being observation times; the second column being fractional mortality at each time in the (control) cohort
##' @return A list with a logical flag to indicate convergence, and the best-fit Weibull shape and scale parameters.
##' @author Pete Dodd
##' @export
getMortParz <- function(M){
  out <- optim(par=c(0,0),fn=function(x)morterr(M,x))
  ans <- list(k.d=exp(out$par[1]),L.d=exp(out$par[2]),converged=TRUE)
  if(abs(out$convergence)>0) ans$converged <- FALSE
  if(!ans$converged) warning('Mortality parameter fitting has not converged!')
  ans
}



## == now for exposure times:

## density for staying alive and being exposed at time x
dns <- function(x,ke,le,km,lm) dweibull(x,shape=ke,scale=le) * exp(-(x/lm)^km)
## normalization for the above
norm <- function(T,ke,le,km,lm) integrate(function(x) dns(x,ke,le,km,lm),lower=0,upper=T)$value
## pdf for exposure time conditional on receiving exposure
normq <- function(T,ke,le,km,lm) norm(T,ke,le,km,lm)/norm(Inf,ke,le,km,lm)

## error for treatment/exposure, as above
experr <- function(M,x){
  x <- exp(x)
  tmp <- apply(M,1,function(y) (normq(y[1],x[1],x[2],x[3],x[4])-y[2])^2) #vector of squared errors
  sum(tmp)
}



##' For fitting to treatment data.
##'
##' This takes data on cumulative percentiles of cumulative treatment among those treated, and fits a Weibull distribution to it so as to minimize the sum-of-squares error, taking into account the competing hazard of death.
##'
##' @title Get the best-fit parameters for treatment
##' @param M A 2xN matrix: the first column being observation times; the second column being cumulative fraction treated at each time in the cohort of those who are ultimately treated
##' @param km The best-fit Weibull shape parameter for mortality
##' @param lm The best-fit Weibull scale parameter for mortality
##' @return A list with a logical flag to indicate convergence, and the best-fit Weibull shape and scale parameters.
##' @author Pete Dodd
##' @export
getTxParz <- function(M,km,lm){
  y <- log(c(km,lm))
  out <- optim(par=c(0,0),fn=function(x) experr(M,c(x,y)))
  ans <- list(k.e=exp(out$par[1]),L.e=exp(out$par[2]),converged=TRUE)
  if(abs(out$convergence)>0) ans$converged <- FALSE
  if(!ans$converged) warning('Treatment parameter fitting has not converged!')
  ans
}
