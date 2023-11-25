#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  ## Your application server logic
  ## TODO wrap these into smaller functions
  ## create data
  k.e <- k.d <- k.l <- 2
  ## L.e <- 1
  ## L.d <- 2
  ## L.l <- 5
  ITBoutput <- reactive({
    ITBstats(N=1e4,
             rtt.exposure = function(n) rweibull(n,k.e,input$L.e),
             rtt.death = function(n) rweibull(n,k.d,input$L.d),
             rtt.ltfu = function(n) rweibull(n,k.l,input$L.l))
  })


  output$distPlot <- renderPlot({
    makeDistPlots(inputs,k.e,k.d,k.l) #TODO
  })

  output$resultPlot <- renderPlot({
    ANS <- ITBoutput()
    makeEffectBySamplePlot(ANS)
  })

  output$tableA1 <- renderTable({
    ANS <- ITBoutput()
    ANS$table.a
  })
  output$tableA2 <- renderTable({
    ANS <- ITBoutput()
    data.table(exposed=ANS$F.e,died=ANS$F.d,RR=ANS$RR.a)
  })
  output$tableB1 <- renderTable({
    ANS <- ITBoutput()
    ANS$table.b
  })
  output$tableB2 <- renderTable({
    ANS <- ITBoutput()
    data.table(exposed=ANS$F.e,died=ANS$F.d,RR=ANS$RR.b)
  })
  output$dists <- renderTable({
    k <- 2
    data.table(event=c('exposure','death','ltfu'),
               mean=c(input$L.e*gamma(1+1/k),input$L.d*gamma(1+1/k),input$L.l*gamma(1+1/k)),
               sd=c(input$L.e,input$L.d,input$L.l)*sqrt(gamma(1+2/k)-gamma(1+1/2)^2),
               `2.5% percentile`=qweibull(0.025,shape=2,scale=c(input$L.e,input$L.d,input$L.l)),
               `25% percentile`=qweibull(0.25,shape=2,scale=c(input$L.e,input$L.d,input$L.l)),
               `50% percentile`=qweibull(0.5,shape=2,scale=c(input$L.e,input$L.d,input$L.l)),
               `75% percentile`=qweibull(0.75,shape=2,scale=c(input$L.e,input$L.d,input$L.l)),
               `97.5% percentile`=qweibull(0.975,shape=2,scale=c(input$L.e,input$L.d,input$L.l))
               )
  })

}
