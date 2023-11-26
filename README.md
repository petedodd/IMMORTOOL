# IMMORTOOL <img src="man/figures/logo.png" align="right" height="139" alt="" />
Tool to explore the potential for immortal time bias


## Installation & use ##


If you have the `devtools` package installed, you can install `IMMORTOOL` with the following commend from R:

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
makeTMplot(input)

```


### Jones & Fowler (using Kumar et al)

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

```


### Kaul

```R
## mortality data and fit
mortality.data <- cbind(c(7,30),c(0.5,0.66)) # 2 points is not ideal
mortality.parms <- getMortParz(mortality.data)

##treatment data and fit
treatment.data <- cbind(c(1),c(0.5) # no data: assume median 1 day
treatment.parms <- getTxParz(treatment.data, mortality.parms$k.d, mortality.parms$L.d)

##plot
input <- c(mortality.parms,treatment.parms)
input$T.max <- 90
makeTMplot(input)

```


### License ###

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

