#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  ## Your application server logic

  ITBoutput <- reactive({
    ITBstats(N=1e4, Tstop=input$T.max,
             rtt.exposure = function(n) rweibull(n,input$k.e,input$L.e),
             rtt.death = function(n) rweibull(n,input$k.d,input$L.d),
             rtt.ltfu = function(n) rweibull(n,1,36500))
  })


  output$distPlot <- renderPlot({
    makeDistPlot(input)
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

}
