
#' @title 
#' UI component for displaying a plot and table output for synergy and efficacy
#' scores
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' server function.
#' 
prediction_ui = function(id) {
  ns = shiny::NS(id)
  
  # Create a Row with two columns for the table and plot output
  shiny::fluidRow(
    shiny::column(width = 7,
                  shiny::div(class="collapse_subentry",
                      DT::dataTableOutput(ns("mytable"), width = "100%"),
                      shiny::h4("Efficacy and Synergy Table"),
                      shiny::HTML("Table of treatment combinations for the",
                                  "selected cell line. Please use the table",
                                  "header for sorting the table and the search",
                                  "input field in the top right to search for",
                                  "keywords. A click on one or multiple rows",
                                  "in the table will highlight the respective",
                                  "treatment combination in the plot on the",
                                  "right side."))
                  ),
    shiny::column(width = 5,
                  shiny::div(class="collapse_subentry",
                      shiny::div(id = "container",
                                 shiny::div(style="float:left;width:100px;",
                                            shiny::actionButton(ns("prvBtn"),
                                                                label = "",
                                                                icon = shiny::icon("caret-left"))),
                             shiny::div(style="float:right;width:100px;",
                                        shiny::actionButton(ns("nextBtn"),label = "",
                                                            icon = shiny::icon("caret-right"))),
                             shiny::div(style="margin:0 auto;width:200px",
                                        shinyWidgets::pickerInput(
                                          ns("aoc_bliss_sample"),
                                          width = "100%",
                                          choices = character(0),
                                          selected = character(0)))
                             ),
                  ggiraph::girafeOutput(ns("aoc_bliss")),
                  shiny::h4("Efficacy and Synergy Plot"),
                  shiny::HTML("Plot of predicted synergy and efficacy scores",
                              "for the selected cell line. Use the input field",
                              "and buttons above the plot to screen different",
                              "cell lines. The list of available cells is ordered",
                              "descreasingly by the correlation coefficient of",
                              "the synergy and efficacy scores.")
                  ))
  )
}

#' @title 
#' Server component for displaying a plot and table output for synergy and 
#' efficacy scores
#' 
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' server function.
#' @param view A character of the fst file listing predictions
#' 
prediction_server = function(id, view= list()) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      # Observer for an update of the prediction reactive and updates the
      # sample selector above the plot
      shiny::observeEvent(efficacy_stat(), {
        models = efficacy_stat() %>%
          dplyr::arrange(dplyr::desc(estimate)) %>%
          dplyr::pull(.identifier_sample_name)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = models)
      })

      # Reactive for identifying top correlated synergy and efficacy scores
      # also handles previous and next action buttons
      efficacy_stat = shiny::reactive({
        shiny::req(!is.null(view))
        
        # calculate correlation for each sample and order
        stat = load_fst(fst_name = view,
                        cols = c(".prediction_aoc", ".prediction_bliss",
                                 ".identifier_sample_name")) %>%
          tidyr::nest(data=-.identifier_sample_name) %>%
          dplyr::mutate(stat = purrr::map(data, 
                                          ~ stats::cor.test(.x$.prediction_aoc,
                                                            .x$.prediction_bliss,
                                                            exact = FALSE,
                                                            method = "spearman")),
                        stat = purrr::map(stat, broom::tidy)) %>%
          tidyr::unnest(stat) %>%
          dplyr::select(-data) %>%
          dplyr::arrange(-estimate)

        shinyjs::disable("prvBtn")
        if (nrow(stat) == 1) {
          shinyjs::disable("nextBtn")
        }

        stat
      }) %>%
        shiny::bindCache(view)

      # Observer for sample update (input above the plot) which handles the index
      # for accessing data for the plot and table and disables/enables 
      # next/previous button
      shiny::observeEvent(input$aoc_bliss_sample, {

        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$.identifier_sample_name == input$aoc_bliss_sample)

        # depending on the possition in the sample list, deactivate/activate
        # buttons above the plot
        shinyjs::disable("prvBtn")
        if (curr_idx > 1)
          shinyjs::enable("prvBtn")

        shinyjs::disable("nextBtn")
        if (curr_idx < nrow(stat))
          shinyjs::enable("nextBtn")
      })

      # Observer for the previous input button which updates the index of the
      # selected sample and the picker input above the plot
      shiny::observeEvent(input$prvBtn, {

        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$.identifier_sample_name == input$aoc_bliss_sample)
        prev_idx = min(1, curr_idx - 1)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = as.character(stat$.identifier_sample_name),
                                        selected = as.character(stat$.identifier_sample_name[prev_idx]))
      })
      
      # Observer for the next input button which updates the index of the
      # selected sample and the picker input above the plot
      shiny::observeEvent(input$nextBtn, {

        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$.identifier_sample_name == input$aoc_bliss_sample)
        next_idx = min(nrow(stat), curr_idx + 1)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = as.character(stat$.identifier_sample_name),
                                        selected = as.character(stat$.identifier_sample_name[next_idx]))
      })

      # Reactive value to sync plot and table data
      selected_rows = shiny::reactiveVal()

      # Observer for the selection of a point in the scatter plot
      shiny::observeEvent(input$aoc_bliss_selected, ignoreNULL = FALSE, {

        shiny::req(!is.na(predictions()))

        curr_sel = predictions() %>%
          dplyr::filter(id %in% !!selected_rows(),
                        .identifier_sample_name == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id)

        deselected = setdiff(curr_sel, input$aoc_bliss_selected)
        if (length(deselected) > 0)
          selected_rows(setdiff(selected_rows(), deselected))

        shiny::req(!is.null(input$aoc_bliss_selected))
        selected_rows(as.numeric(unique(c(selected_rows(), input$aoc_bliss_selected))))
      })

      # Reactive for reading efficacy and synergy predictions for the sample
      # that was selected above the plot and format for the table output
      predictions = shiny::reactive({
        shiny::req(!is.null(input$aoc_bliss_sample))
        preds = load_fst(fst_name = view,
                         cols     = c(".identifier_sample_name",
                                      ".metadata_moa_1", ".metadata_moa_2",
                                      ".metadata_treatment_1",
                                      ".metadata_treatment_2",
                                      ".prediction_bliss",
                                      ".prediction_aoc",
                                      ".metadata_cancer_type",
                                      ".metadata_cancer_subtype"),
                         distinct = TRUE,
                         filter   = list(.identifier_sample_name = input$aoc_bliss_sample))

        preds %>%
          dplyr::ungroup() %>%
          dplyr::mutate(`Cell Line` = .identifier_sample_name,
                        AoC   = .prediction_aoc,
                        Bliss   = .prediction_bliss) %>%
          dplyr::rename(`MoA 1` = .metadata_moa_1,
                        `MoA 2` = .metadata_moa_2,
                        `Treatment 1` = .metadata_treatment_1,
                        `Treatment 2` = .metadata_treatment_2,
                        `Cancer Type` = .metadata_cancer_type,
                        `Cancer Subtype` = .metadata_cancer_subtype) %>%
          dplyr::mutate(id = dplyr::row_number())

      }) %>%
        shiny::bindCache(view,
                         input$aoc_bliss_sample)

      
      # Table output of prediction scores
      output$mytable <- DT::renderDataTable({
        shiny::req(!is.null(predictions()))

        predictions() %>%
          dplyr::mutate(dplyr::across(tidyselect:::where(is.numeric),
                                      ~ round(.x, digits = 2))) %>%
          dplyr::mutate(dplyr::across(tidyselect:::where(is.character),
                                      as.factor)) %>%
          dplyr::select(-dplyr::starts_with(c(".id", ".meta", ".pred")), -id) %>%
          DT::datatable(rownames = FALSE,
                        options = list(scrollX = TRUE, autoWidth = TRUE,
                                       initComplete = DT::JS(
                                         "function(settings, json) {",
                                         "$('tbody').css({'word-break': 'break-all'});",
                                         "}"),
                                       columnDefs = list(list(
                                           targets = "_all",
                                           render = DT::JS(
                                             "function(data, type, row, meta) {",
                                             "return type === 'display' && data != null && data.length > 10 ?",
                                             "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                             "}")
                                         )),
                                       class = "display"))
      })

      # variable needed to identify selected rows
      proxy = DT::dataTableProxy('mytable')

      # Observer for selected rows in the table which updates the selected
      # data, which are synced with the plot
      shiny::observeEvent(input$mytable_rows_selected, ignoreNULL = FALSE, {

        shiny::req(!is.na(predictions()))

        curr_sel = predictions() %>%
          dplyr::filter(id %in% !!selected_rows(),
                        .identifier_sample_name == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id)

        deselected = setdiff(curr_sel, input$mytable_rows_selected)
        if (length(deselected) > 0)
          selected_rows(setdiff(selected_rows(), deselected))

        shiny::req(!is.null(input$mytable_rows_selected))
        selected_rows(as.numeric(unique(c(selected_rows(), input$mytable_rows_selected))))
      })

      # Update selected rows which are selected in the plot
      shiny::observeEvent(selected_rows(),{
        proxy %>%
          DT::selectRows(as.numeric(selected_rows()))
      })

      # Create a plot output
      output$aoc_bliss = ggiraph::renderggiraph({

        shiny::req(!is.null(predictions()),
                   !is.null(input$aoc_bliss_sample),
                   input$aoc_bliss_sample %in% predictions()$.identifier_sample_name)

        plt_df = predictions() %>%
          dplyr::filter(!is.na(.prediction_aoc), !is.na(.prediction_bliss)) %>%
          dplyr::filter(.identifier_sample_name == input$aoc_bliss_sample) %>%
          dplyr::mutate(tool_text = glue::glue("MoA1:\t\t{`MoA 1`}",
                                               "MoA2:\t\t{`MoA 2`}",
                                               "Treatment1:\t{`Treatment 1`}",
                                               "Treatment2:\t{`Treatment 2`}",
                                               .sep = "\n"))
        plt = ggplot2::ggplot(plt_df,
                              ggplot2::aes(x=.prediction_aoc,
                                           y=.prediction_bliss)) +
          ggiraph::geom_point_interactive(ggplot2::aes(tooltip=tool_text,
                                                       data_id = id)) +
          ggplot2::theme_bw(base_size = 16) +
          ggplot2::labs(subtitle = glue::glue("{input$aoc_bliss_sample} (n={nrow(plt_df)})"),
                        x = "Efficacy [AoC]",
                        y = "Synergy [Bliss]") +
          ggplot2::xlim(c(-.2, 1.5)) +
          ggplot2::ylim(c(-.2, 1.5)) +
          ggplot2::theme(aspect.ratio = 1)

        curr_sel = predictions() %>%
          dplyr::filter(id %in% !!selected_rows(),
                        .identifier_sample_name == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id) %>%
          as.character()

        girafe(ggobj = plt,
               options = list(opts_selection(type = "multiple",
                                             selected = curr_sel)))
      })
    })
}
