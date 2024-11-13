#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    fluidPage(
      h1("IMMORTOOL:"),
      h2("Exploring the potential for immortal time bias"),

      ## Main panel for displaying outputs ----
        mainPanel(
          tabsetPanel(type = "tabs",
                      tabPanel(
                        title = "README",
                        value = "S1",
                        h1("Background"),
                        h2("Immortal time bias"),
                        tags$img(src = 'www/illustrations.png',
                                 width = '1000px'## ,
   ## style = 'position: absolute; position: absolute; width: 1024px; height: 768px;'
                                 ),
                        p("Immortal time bias is illustrated in the above diagram. By definition, the outcome cannot occur before the exposure. "),
                        h2("Analysis A"),
                        tags$img(src = 'www/analysisA.png',width = '1000px'),
                        p("The most basic approach is simply to consider all time in those ultimately exposed as exposed. This misclassifies some unexposed time as exposed. In IMMORTOOL output the rate ratio for this analysis is labelled 'a'."),
                        h2("Analysis B"),
                        p("Another approach excludes those with early deaths before a specified time. In IMMORTOOL output the rate ratio for this analysis is labelled 'b'."),
                        h2("Analysis C"),
                        p("As above, this approach excludes those with early deaths before a specified time, but also resets the clock to measure person time from this specified time. In IMMORTOOL output the rate ratio for this analysis is labelled 'c'."),
                        h2("Analysis D"),
                        ## tags$img(src = 'www/analysisD.png',width = '1000px'),
                        p("Landmark analysis attempts to mitigate the effects of immortal time bias in two ways: 1) Those dying before landmark time are excluded; 2) Those exposed after landmark are considered unexposed. For Landmark analysis exposure status at the landmark time is definitive; whereas in analysis c, exposure status and the end time is definitive. In IMMORTOOL output the rate ratio for this analysis is labelled 'd'."),
                        p("NOTE: the Landmark time input is used for analyses b-d"),
                        p(" ")
                        ## code omitted
                      ),
                      tabPanel("Fitting distributions",
                               value="S1",
                               h2("Adjust how many data points there are for mortality/treatment, then press button to fit."),
                               h3("What maximum time should we use for plots and simulation?"),
                               numericInput("Tmax", "Time horizon",
                                            value = 90, min = 1, max = 365),
                               h3("Is the denominator for the treatment data the whole cohort or those ultimately exposed?"),
                               radioButtons("denominator", "Denominator type:",
                                            c("Whole cohort" = "cohort",
                                              "Only those ultimately treated" = "exposed")),
                               h3("Times and fractions should increase (mortality is cumulative)."),
                               fluidRow(
                                 column(width=6,
                                        numericInput("mortnum", "Number of mortality datapoints",
                                                     value = 3, min = 2, max = 5),
                                        uiOutput("morttimeinput"),
                                        uiOutput("mortfracinput")),
                                 column(width=6,
                                        numericInput("TTEnum", "Number of treatment time data points",
                                                     value = 2, min = 1, max = 5),
                                        uiOutput("TTEtimeinput"),
                                        uiOutput("TTEfracinput"))
                               ),
                               actionButton("do", "Try to fit!"),
                               ## textOutput('fits'),
                               h3("(See next tab for fit)")
                               ),
                      tabPanel("Time-to-event distributions", plotOutput(outputId = "FitPlot")),
                      tabPanel("Run simulation",
                               h2("Sliders on fitted values by default. Adjust, then press go!\n"),
                               textOutput("fits"),
                               radioButtons("fitsorsliders", "Use fit or slider values?",
                                            c("Slider values" = "svl",
                                              "Fit values" = "fvl")),
                               h2("\n"),
                               actionButton("runsim", "Run simulation!"),
                               h3("(See next tab for results)"),
                               h2("\n"),
                               ## Input: Slider for End time ----
                               sliderInput(inputId = "T.max",
                                           label = "Follow-up time:",
                                           width='100%',
                                           min = 1,
                                           max = 100,
                                           value = 30),

                               ## Input: Slider for exposure ----
                               sliderInput(inputId = "L.e",
                                           label = "Exposure scale parameter:",
                                           width='100%',
                                           min = 0.001,
                                           max = 500,
                                           value = 10.6),
                               sliderInput(inputId = "k.e",
                                           label = "Exposure shape parameter:",
                                           width='100%',
                                           min = -1,
                                           max = 20,
                                           value = 2.3),

                               ## Input: Slider for death ----
                               sliderInput(inputId = "L.d",
                                           label = "Death scale parameter:",
                                           width='100%',
                                           min = 0.001,
                                           max = 500,
                                           value = 428),
                               sliderInput(inputId = "k.d",
                                           label = "Death shape parameter:",
                                           width='100%',
                                           min = -1,
                                           max = 20,
                                           value = 0.5),

                               ## Input: Time for landmark analysis ---
                               sliderInput(inputId = "T.landmark",
                                           label = "Time for landmark analysis:",
                                           width='100%',
                                           min = 1,
                                           max = 100,
                                           value = 3)

                               ),
                      tabPanel("Simulation results", textOutput("inputtext"), tableOutput("tableA1"))

                      )
        )
      )

  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "IMMORTOOL"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
