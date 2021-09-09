#' @title 
#' UI Component of collabsible box
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#'  serve function.
#' @param ui Shiny ui components to be shown in the collapsible box.
#' 
collapsible_box_ui = function(id, ui) {
  ns = shiny::NS(id)
  
  shiny::tagList(
    shiny::actionLink(ns("collapse"), "", class="btn_collabsible"),
    shiny::div(id = ns("collapse_elements"), class="collapse_subentry",
        ui
    )
  )
}

#' @title 
#' Server component for a collapsible HTML box
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#'  UI function.
#' @param collapsible_title A character that is used as a title for the box
#' 
collapsible_box_server = function(id, collapsible_title="Test") {
  shiny::moduleServer(
    id, function(input, output, session) {

      # Reactive that saves the box state
      open = shiny::reactiveVal(FALSE)

      # Observer for the box header, shows/hides box content
      shiny::observeEvent(input$collapse, ignoreNULL = FALSE, {

        if (!open()){
          html = shiny::div(class="collapse_entry", collapsible_title, style="text-align:left;",
                     shiny::span(shiny::icon("chevron-down"), style="float:right;"))
          shiny::updateActionLink(session, "collapse", label = as.character(html))
          shinyjs::hideElement("collapse_elements", anim = TRUE, animType = "slide")
        } else {
          html = shiny::div(class="collapse_entry", collapsible_title, style="text-align:left;",
                            shiny::span(shiny::icon("chevron-up"), style="float:right;"))
          shiny::updateActionLink(session, "collapse", label = as.character(html))
          shinyjs::showElement("collapse_elements", anim = TRUE, animType = "slide")
        }
        open(!open())
      })
    }
  )
}
