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

      ## Sidebar layout with input and output definitions ----
      sidebarLayout(

        ## Sidebar panel for inputs ----
        sidebarPanel(

          ## Input: Slider for End time ----
          sliderInput(inputId = "T.max",
                      label = "Follow-up time:",
                      min = 1,
                      max = 100,
                      value = 90),

          ## Input: Slider for exposure ----
          sliderInput(inputId = "L.e",
                      label = "Exposure scale parameter:",
                      min = 0.001,
                      max = 500,
                      value = 10.6),
          sliderInput(inputId = "k.e",
                      label = "Exposure shape parameter:",
                      min = -1,
                      max = 20,
                      value = 2.3),

          ## Input: Slider for death ----
          sliderInput(inputId = "L.d",
                      label = "Death scale parameter:",
                      min = 0.001,
                      max = 500,
                      value = 428),
          sliderInput(inputId = "k.d",
                      label = "Death shape parameter:",
                      min = -1,
                      max = 20,
                      value = 0.5),

          ## Input: Time for landmark analysis ---
          sliderInput(inputId = "T.landmark",
                      label = "Time for landmark analysis:",
                      min = 1,
                      max = 100,
                      value = 3)

          ## ## Input: Slider for death ----
          ## sliderInput(inputId = "L.l",
          ##             label = "LTFU scale parameter:",
          ##             min = 2,
          ##             max = 1000,
          ##             value = 5),
          ## sliderInput(inputId = "k.l",
          ##             label = "LTFU shape parameter:",
          ##             min = -1,
          ##             max = 20,
          ##             value = 1)
        ),

        ## Main panel for displaying outputs ----
        mainPanel(
          tabsetPanel(type = "tabs",
                      tabPanel(
                        title = "README",
                        value = "S1",
                        h1("Background"),
                        h2("Explanations of what is happening and how to interpret"),
                        tags$img(src = 'www/illustrations.png',
                                 width = '1000px'## ,
                                 ## style = 'position: absolute; position: absolute; width: 1024px; height: 768px;'
                                 ),
                        p(" Yet more text TODO")
                                        # code omitted
                      ),
                      tabPanel("Fitting distributions",
                               value="S1",
                               h2("Adjust how many data points there are for mortality/treatment, then press button to fit."),
                               h3("Times and fractions should increase (mortality is cumulative)."),
                               fluidRow(
                                 column(width=6,
                                        numericInput("mortnum", "Number of mortality datapoints",
                                                     value = 3, min = 2, max = 5),
                                        uiOutput("morttimeinput"),
                                        uiOutput("mortfracinput")),
                                 column(width=6,
                                        numericInput("TTEnum", "Number of treatment time data points",
                                                     value = 1, min = 1, max = 5),
                                        uiOutput("TTEtimeinput"),
                                        uiOutput("TTEfracinput"))
                               ),
                               actionButton("do", "Try to fit!"),
                               textOutput('fits')
                               ),
                      ## tabPanel("Time-to-event distributions", plotOutput(outputId = "distPlot")),
                      tabPanel("Time-to-event distributions", plotOutput(outputId = "TMPlot")),
                      tabPanel("Distribution stats", tableOutput('dists')),
                      tabPanel("Output Plot", plotOutput(outputId = "resultPlot"))#,
                      ## tabPanel("Rates, method a", tableOutput('tableA1')),
                      ## tabPanel("Rate ratio, method a", tableOutput('tableA2')),
                      ## tabPanel("Rates, method b", tableOutput('tableB1')),
                      ## tabPanel("Rate ratio, method b", tableOutput('tableB2'))
                      )
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
