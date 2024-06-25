#' The application server-side
#'
#' @param input,output,session Internal parameters or {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  ## Your application server logic

  ## ---------- create cohort ---------------
  ITBoutput <- reactive({
    ITBstats(N=1e4, Tstop=input$T.max,
             rtt.exposure = function(n) rweibull(n,input$k.e,input$L.e),
             rtt.death = function(n) rweibull(n,input$k.d,input$L.d),
             rtt.ltfu = function(n) rweibull(n,1,36500))
  })

  ##  --------- fitting --------------------
  ## making boxes for fitting
  output$morttimeinput <- renderUI({
    mortnum <- as.integer(input$mortnum)
    lapply(1:mortnum, function(i) {
      numericInput(inputId = paste0("morttime", i), label = paste("Time ", i),
                   value=i, min = 0, max = 90) #TODO change max
    })
  })
  output$mortfracinput <- renderUI({
    mortnum <- as.integer(input$mortnum)
    lapply(1:mortnum, function(i) {
      numericInput(inputId = paste0("mortfrac", i), label = paste("Fraction dead ", i),
                   value=i/10, min = 0, max = 1)
    })
  })

  ## TTE fitting
  output$TTEtimeinput <- renderUI({
    TTEnum <- as.integer(input$TTEnum)
    lapply(1:TTEnum, function(i) {
      numericInput(inputId = paste0("TTEtime", i), label = paste("Time ", i),
                   value=i, min = 0, max = 90) #TODO change max
    })
  })
  output$TTEfracinput <- renderUI({
    TTEnum <- as.integer(input$TTEnum)
    lapply(1:TTEnum, function(i) {
      numericInput(inputId = paste0("TTEfrac", i), label = paste("Time-to-treatment quantile ", i),
                   value=i/10, min = 0, max = 1)
    })
  })

  ## capture the mortality/treatment times & fracs
  fractimes <- reactive({
    ## mortality fit first
    mt <- mf <- list()
    for(i in 1:input$mortnum){
      mt[[i]] <- input[[paste0('morttime',i)]]
      mf[[i]] <- input[[paste0('mortfrac',i)]]
    }
    mt <- unlist(mt)
    mf <- unlist(mf)
    if(any(!c(mt,mf)>0)) stop('Mortality data not >0!')
    if(!all(mt==cummax(mt))) stop('Mortality times do not increase!')
    if(!all(mf==cummax(mf))) stop('Mortality fractions do not increase!')
    ## treatment
    tt <- tf <- list()
    for(i in 1:input$TTEnum){
      tt[[i]] <- input[[paste0('TTEtime',i)]]
      tf[[i]] <- input[[paste0('TTEfrac',i)]]
    }
    tt <- unlist(tt)
    tf <- unlist(tf)
    if(any(!c(tt,tf)>0)) stop('Mortality data not >0!')
    if(!all(tt==cummax(tt))) stop('Treatment times do not increase!')
    if(!all(tf==cummax(tf))) stop('Treatment fractions do not increase!')
    list(mt=mt,mf=mf,tt=tt,tf=tf)
  })


  ## fit action TODO
  fitVals <- eventReactive(input$do, {
  ## fitVals <- reactive({
    FT <- fractimes()
    mt <- FT$mt; mf <- FT$mf
    tt <- FT$tt; tf <- FT$tf
    ## mortality fit
    D <- Yfit(mt, mf)
    ## then treatment fit
    if(!any(!is.finite(D))){
      TD <- cbind(tt,tf)
      T <- getTxParz(TD,D)
    }
    if(!T$converged){
      txt <- 'Unfortunately, while the mortality converged, the treatment fit did not converge!'
    } else {
      ## output
      txt <- paste0('Suggested parameters based on this are:\n Mortality (shape,scale)=(',
                    D['k'],' , ',D['L'],' )\n Treatment (shape,scale)=(',
                    T$k.e,' , ',T$L.e,' )\n')
    }
    list(txt = txt,
         k.d = D['k'],
         L.d = D['L'],
         k.e = ifelse(!T$converged,NULL,T$k.e),
         L.e = ifelse(!T$converged,NULL,T$L.e)
         )
  })

  ## text to report results
  output$fits <- renderText({
    FV <- fitVals()
    FV$txt
  })


  ## ---- other stuff -------
  output$FitPlot <- renderPlot({
    makeFitCheckPlot(input,fractimes(),fitVals())
  })

  output$dists <- renderTable({
    k <- 2
    L.l <- 36500
    data.table(event=c('exposure','death','ltfu'),
               mean=c(input$L.e*gamma(1+1/k),input$L.d*gamma(1+1/k),36500*gamma(1+1/k)),
               sd=c(input$L.e,input$L.d,36500)*sqrt(gamma(1+2/k)-gamma(1+1/2)^2),
               `2.5% percentile`=qweibull(0.025,shape=2,scale=c(input$L.e,input$L.d,36500)),
               `25% percentile`=qweibull(0.25,shape=2,scale=c(input$L.e,input$L.d,36500)),
               `50% percentile`=qweibull(0.5,shape=2,scale=c(input$L.e,input$L.d,36500)),
               `75% percentile`=qweibull(0.75,shape=2,scale=c(input$L.e,input$L.d,36500)),
               `97.5% percentile`=qweibull(0.975,shape=2,scale=c(input$L.e,input$L.d,36500))
               )
  })
  
  ## output$distPlot <- renderPlot({
  ##   makeDistPlot(input)
  ## })

  ## output$TMPlot <- renderPlot({
  ##   makeTMplot(input)
  ## })

  
  ## ## new approach designed to allow persistence of data between tabs:
  ## v <- reactiveValues(data=NULL) 
 
  ## ## fit action TODO
  ## observeEvent(input$do, {
  ##   ## mortality fit first
  ##   mt <- mf <- list()
  ##   for(i in 1:input$mortnum){
  ##     mt[[i]] <- input[[paste0('morttime',i)]]
  ##     mf[[i]] <- input[[paste0('mortfrac',i)]]
  ##   }
  ##   mt <- unlist(mt)
  ##   mf <- unlist(mf)
  ##   MD <- cbind(mt,mf)
  ##   if(any(!MD>0)) stop('Mortality data not >0!')
  ##   if(!all(mt==cummax(mt))) stop('Mortality times do not increase!')
  ##   if(!all(mf==cummax(mf))) stop('Mortality fractions do not increase!')
  ##   ## D <- getMortParz(MD) #fit TODO old version
  ##   D <- Yfit(mt, mf) #mortality.parms TODO denominator button
  ##   ## then treatment fit
  ##   if(!any(!is.finite(D))){
  ##     tt <- tf <- list()
  ##     for(i in 1:input$TTEnum){
  ##       tt[[i]] <- input[[paste0('TTEtime',i)]]
  ##       tf[[i]] <- input[[paste0('TTEfrac',i)]]
  ##     }
  ##     tt <- unlist(tt)
  ##     tf <- unlist(tf)
  ##     TD <- cbind(tt,tf)
  ##     if(any(!TD>0)) stop('TTE data not >0!')
  ##     if(!all(tt==cummax(tt))) stop('Treatment times do not increase!')
  ##     if(!all(tf==cummax(tf))) stop('Treatment fractions do not increase!')
  ##     ## T <- getTxParz(TD,D$k.d,D$L.d) #old version TODO denominator type
  ##     T <- getTxParz(TD,D)
  ##   }
  ##   ## if(!D$converged){
  ##   ##   txt <- 'Unfortunately, the mortality fit did not converge!'
  ##   ## } else 
  ##   if(!T$converged){
  ##     txt <- 'Unfortunately, while the mortality converged, the treatment fit did not converge!'
  ##   } else {
  ##     ## output
  ##     txt <- paste0('Suggested parameters based on this are:\n Mortality (shape,scale)=(',
  ##                   D['k'],' , ',D['L'],' )\n Treatment (shape,scale)=(',
  ##                   T$k.e,' , ',T$L.e,' )\n')
  ##   }
  ##   ## capture outputs
  ##   v$k.e <- T$k.e; v$L.e <- T$L.e
  ##   v$k.d <- D['k']; v$L.d <- D['L']
  ##   v$txt <- txt
  ## })

  ## output$fits <- renderText({
  ##   if(is.null(v$txt)) return()
  ##   v$txt
  ## })

  ## output$FitPlot <- renderPlot({
  ##   if(is.null(v$txt)) return()
  ##   makeFitCheckPlot(input)
  ## })

  ## output$resultPlot <- renderPlot({
  ##   ANS <- ITBoutput()
  ##   makeEffectBySamplePlot(ANS)
  ## })

  ## output$tableA1 <- renderTable({
  ##   ANS <- ITBoutput()
  ##   ANS$table.a
  ## })
  ## output$tableA2 <- renderTable({
  ##   ANS <- ITBoutput()
  ##   data.table(exposed=ANS$F.e,died=ANS$F.d,RR=ANS$RR.a)
  ## })
  ## output$tableB1 <- renderTable({
  ##   ANS <- ITBoutput()
  ##   ANS$table.b
  ## })
  ## output$tableB2 <- renderTable({
  ##   ANS <- ITBoutput()
  ##   data.table(exposed=ANS$F.e,died=ANS$F.d,RR=ANS$RR.b)
  ## })


}
