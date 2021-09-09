
#' Shiny UI part of the shiny app
#'
#' @return A shiny UI
initialize_ui = function() {
  
  # Define UI for application that draws a histogram
  shiny::navbarPage(
    title = "Synergy DDR",
    theme = bslib::bs_theme(version = 4, bootswatch = "materia"),
    shiny::tabPanel("Welcome",
             shiny::fluidPage(
               shiny::fluidRow(
                 
                 shiny::column(width = 10,
                               shiny::wellPanel(style = "background-color: #fcfcfc;",
                                                shiny::h3("Welcome"),
                                                app_welcome_paragraph()),
                               shiny::br(),
                               shiny::wellPanel(style = "background-color: #fcfcfc;",
                                                shiny::h3("Abstract"),
                                                abstract(),
                                                shiny::p(),
                                                shiny::img(src="modeloutline.png",
                                                           style="width: 35%;display: block;margin-left: auto;margin-right: auto;"))
                 ),
                 shiny::column(width = 2,
                               shiny::wellPanel(style = "background-color: #fcfcfc;",
                                                shiny::h3("NEWS"),
                                                app_news_list()))
               )
             )),
    shiny::navbarMenu(
      "Model",
      shiny::tabPanel(
        "Input data",
        shiny::fluidPage(
          shiny::fluidRow(
            shiny::column(
              width=12,
              nested_plot_ui(
                "input", x_breaks_name="Feature",
                info_title = shiny::h4("Input Feature Overview"),
                info_text = shiny::HTML(
                  "This plot provides an overview of the different feature",
                  "categories and the number of features within each category",
                  "that are use for our surrogate model. A click on one of the",
                  "categories opens another plot with a more detailed view on",
                  "the selected category. Click on 'Show data' to see individual",
                  "statistics of a feature or feature group of interest.")
              ))
          )
        )),
      shiny::tabPanel("Prediction",
               shiny::fluidPage(
                 prediction_ui("aoc_bliss_prediction")
               )),
      shiny::tabPanel("SHAP analysis",
               shiny::fluidPage(
                 shiny::fluidRow(
                   shiny::column(
                     width=6,
                     nested_plot_ui("aoc", x_breaks_name="Feature",
                                    info_title = shiny::h4("Efficacy"),
                                    info_text = shiny::HTML(
                                      "Overview of features and their",
                                      "individual contribution to the",
                                      "overall prediction. A click",
                                      "on one of the categories opens",
                                      "another plot with a detailed",
                                      "view of the selected category.",
                                      "Click on 'Show data' to see",
                                      "individual SHAP values of a",
                                      "feature or feature group within",
                                      "the five cross validatation sets."))
                     
                   ),
                   shiny::column(
                     width=6,
                     nested_plot_ui("bliss", x_breaks_name="Feature",
                                    info_title = shiny::h4("Synergy"),
                                    info_text = shiny::HTML(
                                      "Overview of features and their",
                                      "individual contribution to the",
                                      "overall prediction. A click",
                                      "on one of the categories opens",
                                      "another plot with a detailed",
                                      "view of the selected category.",
                                      "Click on 'Show data' to see",
                                      "individual SHAP values of a",
                                      "feature or feature group within",
                                      "the five cross validatation sets."))
                   )
                 )
               )
      )),
    shiny::tabPanel("Predict",
             score_data_ui("score_data")
    ),
    shiny::navbarMenu("Help",
               shiny::tabPanel("FAQ",
                        shiny::tags$iframe(src="faq.html",
                                    style="height: 100vh; width: 100%",
                                    frameBorder="0")
               ),
               shiny::tabPanel("Contact")),
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    shinyjs::useShinyjs()
  )
  
}