#' @title 
#' Welcome text
#' 
#' @description 
#' Welcome text for landing page
#' 
app_welcome_paragraph = function() {
  shiny::HTML("This shiny app provides access to data of our recent publication.",
              "The main goal of the app is to provide you insights into the",
              "input data of our model (Model > Input data), the predicted",
              "synergy and efficacy scores (Model > Prediction) and the",
              "downstream evauation of feature contributions using a SHAP",
              "analysis (Model > Shap analysis).")
}
