## the non-shiny functions
##' Function to simulate a cohort
##'
##' This function will simulate a cohort with relevant times-to-event under a null assumption for the treatment effect. Typical use returns summary statistics.
##'
##' @title ITBstats
##' @param N The size of the cohort for simulation (unrelated to the experimental study cohort size)
##' @param Tstop The time at which observation stops
##' @param Tlandmark The time used for a landmark analysis
##' @param Texc The time used for exclude early events & time analysis
##' @param rtt.exposure random time-to-exposure generator
##' @param rtt.death random time-to-death generator
##' @param rtt.ltfu random time-to-ltfu generator
##' @param returnraw Logical (default FALSE) controlling return type. If true, a list comprising the raw cohort and the raw landmark cohort will be returned. Otherwise, a list of summary statistics will be returned
##' @return Depends on [returnraw].
##'
##' Note that if returnraw==FALSE, the following labels are used for analysis variants:
##' a) person time from 0 (standard); b) person time from exposure; c) landmark analysis: drop those dead/ltfu before landmark AND exposed after landmark -> control; d) excluding early events & reset clock
##'
##' @author Pete Dodd
##' @export
##' @import data.table
##' @import gmodels
ITBstats <- function(N=1e4,           #simulation cohort size
                     Tstop=90,        #end time
                     Tlandmark=1,     #time used in landmark analysis
                     Texc=0,     #time for exclude early deaths analysis
                     rtt.exposure,    #random time-to-exposure
                     rtt.death,       #random time-to-death
                     rtt.ltfu,        #random time-to-ltfu
                     returnraw=FALSE, #return TZ
                     ...
                     ){
  ## create cohort
  TZ <- data.table::data.table(t.e=rtt.exposure(N),
                   t.d=rtt.death(N),
                   t.l=rtt.ltfu(N))
  ## NOTE these exposed/dead categories are used in analyses a & b & d, reset in landmark (c)
  TZ[,exposed:=ifelse(t.e<pmin(t.d,t.l,Tstop),TRUE,FALSE)]
  TZ[,died:=ifelse(t.d<pmin(t.l,Tstop),TRUE,FALSE)]

  ## print
  gmodels::CrossTable(TZ$exposed,TZ$died)

  ## now use person time to death by group to calculate death rate

  ## a) person time from 0
  TZ[,PTa:=pmin(t.d,t.l,Tstop)]
  NPT <- TZ[,.(deaths=sum(died),PT=sum(PTa)),by=exposed]
  NPT[,rate:=deaths/PT]
  RR2 <- NPT[exposed==TRUE]$rate / NPT[exposed==FALSE]$rate
  frac.d.control.a <- TZ[died == TRUE, mean(!exposed)] # frac deaths in control

  ## Suissa stats
  suissa.p <- TZ[exposed==TRUE,sum(t.e)/sum(PTa)] #Suissa p fraction of T1 that us U

  ## b) person time from exposure
  TZ[,PTb:=ifelse(exposed==TRUE,PTa-t.e,PTa)]
  NPT2 <- TZ[,.(deaths=sum(died),PT=sum(PTb)),by=exposed]
  NPT2[,rate:=deaths/PT]
  RR3 <- NPT2[exposed==TRUE]$rate / NPT2[exposed==FALSE]$rate
  frac.d.control.b <- TZ[died == TRUE, mean(!exposed)] # frac deaths in control

  ## c) landmark analysis: drop those dead/ltfu before landmark AND exposed after landmark -> control
  TZL <- TZ[!(t.d < Tlandmark | t.l < Tlandmark)]
  cat('Landmark dropping ', nrow(TZ)-nrow(TZL),' patients from ',nrow(TZ),'\n')
  TZL[,exposed:=ifelse(t.e<pmin(Tlandmark,t.d,t.l,Tstop),TRUE,FALSE)] #should be only t.e<Tlandmark
  TZL[,died:=ifelse(t.d<pmin(t.l,Tstop),TRUE,FALSE)]
  TZL[,PT:=pmin(t.d,t.l,Tstop)-Tlandmark]
  NPTL <- TZL[,.(deaths=sum(died),PT=sum(PT)),by=exposed]
  NPTL[,rate:=deaths/PT]
  RRL <- NPTL[exposed==TRUE]$rate / NPTL[exposed==FALSE]$rate
  frac.d.control.c <- TZL[died == TRUE, mean(!exposed)] # frac deaths in control

  ## d) excluding early events & reset clock
  TZE <- TZ[!(t.d < Texc | t.l < Texc)]
  cat('Exclude early events dropping ', nrow(TZ)-nrow(TZE),' patients from ',nrow(TZ),'\n')
  TZE[,PT:=pmin(t.d,t.l,Tstop)-Texc]            #NOTE resetting clock
  NPTE <- TZE[,.(deaths=sum(died),PT=sum(PT)),by=exposed]
  NPTE[,rate:=deaths/PT]
  RRE <- NPTE[exposed==TRUE]$rate / NPTE[exposed==FALSE]$rate
  frac.d.control.d <- TZE[died == TRUE, mean(!exposed)] # frac deaths in control

  ## return
  if(!returnraw){
    list(table.a=NPT,RR.a=RR2,frac.d.control.a=frac.d.control.a,  #person time from 0
         table.b=NPT2,RR.b=RR3,frac.d.control.b=frac.d.control.b, #person time from exposure
         table.c=NPTL,RR.c=RRL,frac.d.control.c=frac.d.control.c, #landmark
         table.d=NPTE,RR.d=RRE,frac.d.control.d=frac.d.control.d, #exclude early events
         F.e=TZ[,mean(exposed)],F.d=TZ[,mean(died)], #NOTE ltfu=1-death
         suissa.k=NPT[exposed==TRUE]$PT / NPT[exposed==TRUE]$PT,
         suissa.p=suissa.p
         ## TODO document this: not sure relevant to onward calx
         ## https://sphweb.bumc.bu.edu/otlt/mph-modules/ep/ep713_randomerror/ep713_randomerror4.html
         )
  } else {
    return(list(cohort=TZ,landmark.cohort=TZL,excearlyevent.cohort=TZE))
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
                           args=list(shape=input$k.e,scale=input$L.e),n=1e3) +
    ggplot2::geom_function(ggplot2::aes(colour="death"),fun=dweibull,
                           args=list(shape=input$k.d,scale=input$L.d),n=1e3) +
    ggplot2::geom_function(ggplot2::aes(colour="LTFU"),fun=dweibull,
                           args=list(shape=input$k.l,scale=input$L.l),n=1e3) +
    ggplot2::xlab('Time') + ggplot2::ylab('Hazard') +
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



##' Plot survival and treatment distributions.
##'
##' This is a utility to check fitted survival distributions against data for plausibility.
##' 
##' @title makeFitCheckPlot
##' @param input a list of parameters including Weibull treatment distribution parameters (k.e, L.e), Weibull mortality distribution parameters (k.d, L.d), and the maximum time (T.max), denominator type, and the mortality and treatment data.
##' @return A plot
##' @author Pete Dodd
##' @export
makeFitCheckPlot <- function(input){
  tz <- seq(from=0,to=input$Tmax,by=0.1)
  ttl <- ifelse(input$denominator=='cohort',
                'Treatment denominator = whole cohort',
                'Treatment denominator = those ultimately treated')
  plot(input$mortality.times, input$mortality.fracs,
       type='b',xlim=c(0,input$Tmax),ylim=c(0,1),
       xlab='Time',ylab='Fraction',main=ttl)
  lines(tz,1-exp(-(tz/input$L.d)^input$k.d),lty=2)
  lines(input$treatment.times, input$treatment.fracs, type = "b",col=2)
  CE <- tz #cumulative exposure
  if(input$denominator=='cohort'){
    for(i in 1:length(CE)) CE[i] <- norm(tz[i],input$k.e,input$L.e,input$k.d,input$L.d)
  } else{
    for(i in 1:length(CE)) CE[i] <- normq(tz[i],input$k.e,input$L.e,input$k.d,input$L.d)
  }
  ## CE <- 1 - exp(-(tz / input$L.e)^input$k.e)
  lines(tz, CE,col=2,lty=2)
  legend('topleft',## 1, 0.9,
         legend = c("death", "treatment"),
         col = c("black","red"),lty=1,pch=1
         )
}

##' Fit Weibull to times and fractions
##'
##' This uses an OLS fit to transformed Weibull parameters to fit. See Weibull distribution Wikipedia page for explanation.
##' 
##' @title Fit Weilbull
##' @param t vector of time points
##' @param F vector of fractions with event
##' @return vector (k,L) of Weibull shape/scale parameters
##' @author Pete Dodd
##' @export
Yfit <- function(t, F) {
  Y <- log(-log(1 - F))
  T <- log(t)
  mdl <- lm(data = data.frame(T, Y), Y ~ .)
  K <- unname(coef(mdl))
  c(k = K[2], L = exp(-K[1] / K[2]))
}


##' For fitting to treatment data.
##'
##' This takes data on cumulative percentiles of cumulative treatment among those treated, and fits a Weibull distribution to it so as to minimize the sum-of-squares error, taking into account the competing hazard of death.
##'
##' @title Get the best-fit parameters for treatment
##' @param M A 2xN matrix: the first column being observation times; the second column being cumulative fraction treated at each time in the cohort of those who are ultimately treated
##' @param mort.parms The best-fit Weibull shape and scale parameter for mortality as vector
##' @param denominator Whether fractions have the whole cohort as denominator (default: 'cohort') or only those ultimately 'exposed'
##' @return A list with a logical flag to indicate convergence, and the best-fit Weibull shape and scale parameters.
##' @author Pete Dodd
##' @export
getTxParz <- function(M,mort.parms,denominator='cohort'){
  y <- log(mort.parms)
  initial.guess <- log(Yfit(M[,1],M[,2])) #non-competing, denominator=cohort version
  if(any(!is.finite(initial.guess))) initial.guess <- c(0,0)
  out <- optim(par=initial.guess,fn=function(x) experr(M,c(x,y),denominator=denominator))
  ans <- list(k.e=exp(out$par[1]),L.e=exp(out$par[2]),converged=TRUE)
  if(abs(out$convergence)>0) ans$converged <- FALSE
  if(!ans$converged) warning('Treatment parameter fitting has not converged!')
  ans
}

## density for staying alive and being exposed at time x (competing hazards)
dns <- function(x,ke,le,km,lm) dweibull(x,shape=ke,scale=le) * exp(-(x/lm)^km)
## normalization for the above: ie cumulative density
norm <- function(T,ke,le,km,lm){
  if(ke<1 & T < le/200){ #safety: avoid numerical integration over blow-up
    1-exp(-(T/le)^ke) #non-competing version which is close
  } else{
    integrate(function(x) dns(x,ke,le,km,lm),lower=0,upper=T)$value
  }
}
## pdf for exposure time conditional on receiving exposure
normq <- function(T,ke,le,km,lm) norm(T,ke,le,km,lm)/norm(Inf,ke,le,km,lm)
## error for treatment/exposure, as above
experr <- function(M,x,denominator='cohort'){
  x <- exp(x)
  if( !(denominator=='cohort' | denominator=='exposed') )
    stop("denominator must either by 'cohort' or 'exposed'")
  if(denominator=='cohort')
    tmp <- apply(M,1,function(y) (norm(y[1],x[1],x[2],x[3],x[4])-y[2])^2) #vector of squared errors
  if(denominator=='exposed')
      tmp <- apply(M,1,function(y) (normq(y[1],x[1],x[2],x[3],x[4])-y[2])^2) #vector of squared errors
  sum(tmp)
}



##' For computing a factor for approximate IRR CIs
##'
##' For computing a factor for approximate IRR CIs.
##' @title CIfactor
##' @param N the number of deaths observed
##' @param frac.deaths.control the fraction of deaths that are expected in the control arm (returned by ITBstats)
##' @return A factor for computation of approximate CIs. High = mid * factor; Low = mid / factor.
##' @author Pete Dodd
##' @export
CIfactor <- function(N, frac.deaths.control) {
  a <- N * frac.deaths.control
  b <- N * (1 - frac.deaths.control)
  exp(1.96 * sqrt(1 / a + 1 / b)) # for approx CIs
}


##' A shortcut function to generate all results
##'
##' This is a shortcut function to fit to mortality and treatment data and return results. This uses a cohort of 10,000 by default.
##' @title Make Results List
##' @param mortality.times a vector of mortality data times
##' @param mortality.fracs a vector of mortality fractions among all patients corresponding to times
##' @param treatment.times a vector of treatment data times
##' @param treatment.fracs a vector of corresponding treatment fractions among all patients
##' @param N the number of deaths observed (for CIs)
##' @param simulation.cohort.size the size of the cohort used in the simulation
##' @param Tmax the maximum time horizon
##' @param Tlandmark the time used in landmark analysis
##' @param Tearly the time used in the drop-early-events & time analysis
##' @param plotfit Plot the Weibull fits?
##' @param denominator Whether fractions have the whole cohort as denominator (default: 'cohort') or only those ultimately 'exposed'
##' @param ... any extras
##' @param frac.deaths.control the fraction of deaths that are expected in the control arm (returned by ITBstats)
##' @return A list of results from ITBstats together with the number of deaths
##' @author Pete Dodd
##' @export
makeResultList <- function(mortality.times, mortality.fracs,
                           treatment.times, treatment.fracs,
                           N,simulation.cohort.size=1e5,
                           Tmax = 30, Tlandmark = 1, Tearly = 1,
                           plotfit=FALSE,denominator='cohort',
                           ...) {

  ## mortality data and fit
  mortality.parms <- Yfit(mortality.times, mortality.fracs)

  ## treatment data and fit
  treatment.data <- cbind(treatment.times, treatment.fracs)
  treatment.parms <- getTxParz(treatment.data, mortality.parms, denominator=denominator)

  ## combine
  input <- list(k.d=mortality.parms[1], L.d=mortality.parms[2],
             k.e=treatment.parms$k.e, L.e=treatment.parms$L.e)
  input$T.max <- Tmax

  if(plotfit){
    tz <- seq(from=0,to=Tmax,by=0.1)
    ttl <- ifelse(denominator=='cohort',
                  'Treatment denominator = whole cohort',
                  'Treatment denominator = those ultimately treated')
    plot(mortality.times, mortality.fracs,
         type='b',xlim=c(0,Tmax),ylim=c(0,1),
         xlab='Time',ylab='Fraction',main=ttl)
    lines(tz,1-exp(-(tz/input$L.d)^input$k.d),lty=2)
    lines(treatment.times, treatment.fracs, type = "b",col=2)
    CE <- tz #cumulative exposure
    if(denominator=='cohort'){
      for(i in 1:length(CE)) CE[i] <- norm(tz[i],input$k.e,input$L.e,input$k.d,input$L.d)
    } else{
      for(i in 1:length(CE)) CE[i] <- normq(tz[i],input$k.e,input$L.e,input$k.d,input$L.d)
    }
    ## CE <- 1 - exp(-(tz / input$L.e)^input$k.e)
    lines(tz, CE,col=2,lty=2)
    legend('topleft',## 1, 0.9,
      legend = c("death", "treatment"),
      col = c("black","red"),lty=1,pch=1
    )
  }

  # run cohort
  ans <- ITBstats(
    N = simulation.cohort.size,
    Tstop = Tmax, Tlandmark = Tlandmark, Texc = Tearly,
    rtt.exposure = function(n) rweibull(n, input$k.e, input$L.e),
    rtt.death = function(n) rweibull(n, input$k.d, input$L.d),
    rtt.ltfu = function(n) rweibull(n, 1, 36500),...
  )
  ans$N <- N # add in total observed deaths

  ## return
  ans
}


##' For generating a results table from a list of individual study results
##'
##' Generates a results table including CIs for analyses a)-d) of ITB stats as well as the fraction of deaths expected to occur in the control arm.
##' @title Make a table of results
##' @param list.of.results a named list of results from makeResultList, ie from ITBstats with N=#deaths appended
##' @return a data frame of results for all elements in the list
##' @author Pete Dodd
##' @export
makeCItable <- function(list.of.results){
  ## names
  nmz <- names(list.of.results)
  ## numbers
  Nz <- unlist(lapply(list.of.results,function(X)X[['N']]))
  ## mid-points
  Az <- unlist(lapply(list.of.results,function(X)X[['RR.a']]))
  Bz <- unlist(lapply(list.of.results,function(X)X[['RR.b']]))
  Cz <- unlist(lapply(list.of.results,function(X)X[['RR.c']]))
  Dz <- unlist(lapply(list.of.results,function(X)X[['RR.d']]))
  ## death fractions
  faz <- unlist(lapply(list.of.results,function(X)X[['frac.d.control.a']]))
  fbz <- unlist(lapply(list.of.results,function(X)X[['frac.d.control.b']]))
  fcz <- unlist(lapply(list.of.results,function(X)X[['frac.d.control.c']]))
  fdz <- unlist(lapply(list.of.results,function(X)X[['frac.d.control.d']]))
  ## CI factors
  FAZ <- CIfactor(Nz,faz)
  FBZ <- CIfactor(Nz,fbz)
  FCZ <- CIfactor(Nz,fcz)
  FDZ <- CIfactor(Nz,fdz)
  ## formatted CIs
  A <- paste0(signif(Az,2)," (",signif(Az/FAZ,2),"-",signif(Az*FAZ,2),")") #format CIs
  B <- paste0(signif(Bz,2)," (",signif(Bz/FBZ,2),"-",signif(Bz*FBZ,2),")") #format CIs
  C <- paste0(signif(Cz,2)," (",signif(Cz/FCZ,2),"-",signif(Cz*FCZ,2),")") #format CIs
  D <- paste0(signif(Dz,2)," (",signif(Dz/FDZ,2),"-",signif(Dz*FDZ,2),")") #format CIs
  ## answer
  tab <- data.frame(id=nmz,N=Nz,
                    RR.a=Az,RR.b=Bz,RR.c=Cz,RR.d=Dz,
                    frac.d.control.a=faz,frac.d.control.b=fbz,
                    frac.d.control.c=fcz,frac.d.control.d=fdz,
                    F.a=FAZ,F.b=FBZ,F.c=FCZ,F.d=FDZ,
                    CI.a=A,CI.b=B,CI.c=C,CI.d=D)
  rownames(tab) <- NULL
  tab
}
