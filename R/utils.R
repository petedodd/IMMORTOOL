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
makeDistPlot <- function(inputs, #TODO bring other parms under inputs
                         k.e,k.d,k.l){
  ggplot2::ggplot() +
    ggplot2::xlim(0,2*max(input$L.e,input$L.d,input$L.l)) +
    ggplot2::geom_function(aes(colour="exposure"),fun=dweibull,args=list(shape=k.e,scale=input$L.e)) +
    ggplot2::geom_function(aes(colour="death"),fun=dweibull,args=list(shape=k.d,scale=input$L.d)) +
    ggplot2::geom_function(aes(colour="LTBU"),fun=dweibull,args=list(shape=k.l,scale=input$L.l)) +
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
  ggplot2::ggplot(pfd,aes(x=n,y=value,ymin=RRlo,ymax=RRhi,col=analysis,group=analysis))+
    ggplot2::geom_line()+## geom_pointrange()+
    ggplot2::geom_line(aes(y=RRhi),lty=2)+
    ggplot2::geom_line(aes(y=RRlo),lty=2)+
    ggplot2::xlab('Sample size')+ggplot2::ylab('Rate ratio')+
    ggplot2::geom_hline(yintercept=1,col='black',lty=3)+
    ggplot2::theme(legend.position='top')
}
