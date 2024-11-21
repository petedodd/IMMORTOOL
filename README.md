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

## Use online ##

If you want to use without installing, you can access a shiny-live version of this app at:

https://petedodd.github.io/IMMORTOOL-live/

Please note: this will take some moments to load. It may be necessary to access this site using Firefox or chrome browsers.


## Examples of command line use ##

It may be preferable to use the functions provided in this package directly from R, rather than running the shiny app. This is more flexible and reproducible. To access examples of use, please build the vignette on installation, e.g.

```R
devtools::install_github('petedodd/IMMORTOOL',build_vignettes=TRUE)
```
after which the vignette can be accessed in the usual way from R:

```R
vignette('IMMORTOOL-use-from-R')
```

Alternatively, you can access the vignette Rmd file in the vignettes folder of this repository, and create the vignette from R using:

```R
 rmarkdown::render('IMMORTOOL-use-from-R.Rmd')
```
which will take around a minute.






### License ###

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

