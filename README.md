# IMMORTOOL <img src="man/figures/logo.png" align="right" height="139" alt="" />
Tool to explore the potential for immortal time bias


## Installation & use ##


If you have the `devtools` package installed, you can install `IMMORTOOL` with the following command from R:

```R
devtools::install_github('petedodd/IMMORTOOL',dependencies=TRUE)
```

The library can then be loaded ahead of use. The shiny app can be started as follows:
```R
library(IMMORTOOL)
run_app()
```

Depending on your system, this may open the app in a system web browser. You may need to stop the app from R by typing Cntrl-C or otherwise exiting the function.

Apart from the shiny app, a number of utility functions are provided that can be used without the GUI: see below and additional documentation.

## Examples


### van der Vaart et al.

```R
## mortality data and fit
mortality.data <- cbind(c(7,30,90),c(25/298,76/298,103/298)) #col1 = times; col2 = cumulative deaths
mortality.parms <- getMortParz(mortality.data)

##treatment data and fit
treatment.data <- cbind(c(6,9,12),c(0.25,0.5,0.75)) # median & IQR timings for those treated
treatment.parms <- getTxParz(treatment.data, mortality.parms$k.d, mortality.parms$L.d)

##plot
input <- c(mortality.parms,treatment.parms)
input$T.max <- 90
makeTMplot(input) #NOTE unlike fitting & simulation, the treatment timing plot excludes competing mortality

## run cohort
ans.vdV <- ITBstats(N=1e4,Tstop=90,Tlandmark=1,
                    rtt.exposure = function(n) rweibull(n,input$k.e,input$L.e),
                    rtt.death = function(n) rweibull(n,input$k.d,input$L.d),
                    rtt.ltfu = function(n) rweibull(n,1,36500))
N.vdV <- list(a=103,b=44) #event counts (deaths)
```


### Jones & Fowler (using Kumar et al.)

```R
## mortality data and fit
mortality.data <- cbind(c(14,28,90),c(0.107,0.143,0.173))
mortality.parms <- getMortParz(mortality.data)

##treatment data and fit
treatment.data <- cbind(c(1),c((540-37)/540)) # 37 of treated > 1d 
treatment.parms <- getTxParz(treatment.data, mortality.parms$k.d, mortality.parms$L.d)

##plot
input <- c(mortality.parms,treatment.parms)
input$T.max <- 90
makeTMplot(input)

## run cohort
ans.JnF <- ITBstats(N=1e4,Tstop=90,Tlandmark=1,
                    rtt.exposure = function(n) rweibull(n,input$k.e,input$L.e),
                    rtt.death = function(n) rweibull(n,input$k.d,input$L.d),
                    rtt.ltfu = function(n) rweibull(n,1,36500))
N.JnF <- list(a=12,b=105) #event counts

```


### Kaul et al.

```R
## mortality data and fit
mortality.data <- cbind(c(7,30),c(0.5,0.66)) # 2 points is not ideal
mortality.parms <- getMortParz(mortality.data)

##treatment data and fit
treatment.data <- cbind(c(1),c(0.5)) # no data: assume median 1 day
treatment.parms <- getTxParz(treatment.data, mortality.parms$k.d, mortality.parms$L.d)

##plot
input <- c(mortality.parms,treatment.parms)
input$T.max <- 90
makeTMplot(input)

## run cohort
ans.Kaul <- ITBstats(N=1e4,Tstop=30,Tlandmark=1,
                    rtt.exposure = function(n) rweibull(n,input$k.e,input$L.e),
                    rtt.death = function(n) rweibull(n,input$k.d,input$L.d),
                    rtt.ltfu = function(n) rweibull(n,1,36500))

N.Kaul <- list(a=11,b=14)

```

### Combined:

Looking at these together:

```R
## combine answers:
CIfac <- function(L) exp(1.96*sqrt(1/L$a+1/L$b)) #for approx CIs
F <- unlist(lapply(list(N.vdV,N.JnF,N.Kaul),CIfac))
A <- c(ans.vdV$RR.a,ans.JnF$RR.a,ans.Kaul$RR.a)
B <- c(ans.vdV$RR.b,ans.JnF$RR.b,ans.Kaul$RR.b)
A <- paste0(signif(A,2)," (",signif(A/F,2),"-",signif(A*F,2),")") #format CIs
B <- paste0(signif(B,2)," (",signif(B/F,2),"-",signif(B*F,2),")")

## table:
ANS <- data.frame(author=c('vdv','JnF','Kaul'),A = A,B = B)

```



### License ###

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

