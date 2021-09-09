#' @title 
#' UI component for displaying a plot and table output for synergy and efficacy
#' scores of an uploaded file
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' server function.
#' 
score_data_ui = function(id) {
  ns = shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      shiny::column(width = 12,
                    shiny::div(class="collapse_subentry",
                        shiny::h4("Scoring User Dataset"),
                        shiny::HTML("This section allows you to load molecular",
                                    "data of your cell lines of interest and",
                                    "predict synergy and efficacy scores for",
                                    "the analysed drug combinations. Please",
                                    "select a file using the 'Select' button",
                                    "and click on 'Score'. The file must be",
                                    "formated in the same format as the example file",
                                    "Please keep the number of cell lines at a",
                                    "minimum, them runtime depends on your system."),
                        shiny::div(style="margin-top: 40px;"),
                        shiny::fileInput(
                          ns("file1"),
                          "Select .tsv file with molecular data",
                          multiple = FALSE,
                          buttonLabel = "Select",
                          placeholder="Select file",
                          accept = c("text/table/tsv",
                                     "text/tab-separated-values,text/plain",
                                     ".tsv")),
                        shiny::downloadLink(ns("example_file_download"),
                                            "Example File")
                    )

      )
    ),
    shiny::p(),
    shiny::uiOutput(ns("result_ui"))
  )
}

#' @title 
#' Server component for displaying a plot and table output for synergy and 
#' efficacy scores for an uploaded file
#' 
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' server function.
#' 
score_data_server = function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      example_df = shiny::reactive({
        system.file("extdata/example_input.tsv",
                    package = "guanlabddrdrugcombination") %>%
          read.table(header = TRUE, sep = "\t")
      })

      output$example_file_download <- shiny::downloadHandler(
          filename = "example_molecular_data.tsv",
          content = function(file) {
            write.table(example_df(), file, sep = "\t", row.names = FALSE)
            }
        )

      valid_df = shiny::reactive({
        file       = input$file1
        extension  = tools::file_ext(file$datapath)

        shiny::validate(shiny::need(!is.null(file),
                                    "Please upload data"))
        shiny::validate(shiny::need(extension == "tsv",
                                    "Please upload a tsv file"))

        df = readr::read_tsv(file$datapath, col_names = TRUE,
                             show_col_types = FALSE)

        # verify correct columns
        genes = unique(c(genes_aoc, genes_bliss))
        missing_genes = setdiff(genes, colnames(df))
        shiny::validate(shiny::need(length(missing_genes) == 0,
                                    glue::glue("Could not find all required",
                                               "genes in uploaded data. Please",
                                               "add: {paste(missing_genes,",
                                               "collapse=', ')}", .sep = " ")))

        shiny::validate(shiny::need("sample" %in% colnames(df),
                                    paste("Could not find sample names in",
                                          "uploaded file. Please add a sample",
                                          " column to your data")))

        df
      })

      prediction = shiny::reactive({

        shiny::validate(shiny::need(valid_df(),
                                    "No valid data uploaded"))

        predict_optimal_drug_combination(model_path = system.file("models",
                                                                  package = "guanlabddrdrugcombination"),
                                         mol_df=valid_df())
      })

      output$result_ui = shiny::renderUI({

        ns = session$ns
        shiny::validate(shiny::need(prediction(),
                                    "Invalid prediction, please contact the app author"))
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
      })


      shiny::observeEvent(efficacy_stat(), {

        models = efficacy_stat() %>%
          dplyr::arrange(dplyr::desc(estimate)) %>%
          dplyr::pull(sample)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = models)
      })

      efficacy_stat = shiny::reactive({
        stat = prediction()$prediction %>%
          tidyr::nest(data=-sample) %>%
          dplyr::mutate(stat = purrr::map(data, ~ cor.test(.x$AOC,
                                                           .x$Bliss,
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
      })

      shiny::observeEvent(input$aoc_bliss_sample, {

        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$sample == input$aoc_bliss_sample)

        shinyjs::disable("prvBtn")
        if (curr_idx > 1)
          shinyjs::enable("prvBtn")

        shinyjs::disable("nextBtn")
        if (curr_idx < nrow(stat))
          shinyjs::enable("nextBtn")
      })

      shiny::observeEvent(input$prvBtn, {


        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$sample == input$aoc_bliss_sample)
        prev_idx = min(1, curr_idx - 1)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = as.character(stat$sample),
                                        selected = as.character(stat$sample[prev_idx]))
      })

      shiny::observeEvent(input$nextBtn, {


        shiny::req(!is.null(input$aoc_bliss_sample))

        # find the current selected value
        stat = efficacy_stat()

        curr_idx = which(stat$sample == input$aoc_bliss_sample)
        next_idx = min(nrow(stat), curr_idx + 1)

        shinyWidgets::updatePickerInput(session, "aoc_bliss_sample",
                                        choices = as.character(stat$sample),
                                        selected = as.character(stat$sample[next_idx]))
      })


      selected_rows = shiny::reactiveVal()

      shiny::observeEvent(input$aoc_bliss_selected, ignoreNULL = FALSE, {

        shiny::req(!is.na(predictions()))

        curr_sel = predictions() %>%
          dplyr::filter(id %in% !!selected_rows(),
                        sample == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id)

        deselected = setdiff(curr_sel, input$aoc_bliss_selected)
        if (length(deselected) > 0)
          selected_rows(setdiff(selected_rows(), deselected))

        shiny::req(!is.null(input$aoc_bliss_selected))
        selected_rows(as.numeric(unique(c(selected_rows(), input$aoc_bliss_selected))))
      })


      predictions = shiny::reactive({
        shiny::req(!is.null(input$aoc_bliss_sample))

        preds = prediction()$prediction %>%
          dplyr::filter(sample == input$aoc_bliss_sample)

        preds %>%
          dplyr::ungroup() %>%
          dplyr::mutate(`Cell Line` = sample) %>%
          dplyr::mutate(id = dplyr::row_number())

      })

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
                                         render = JS(
                                           "function(data, type, row, meta) {",
                                           "return type === 'display' && data != null && data.length > 10 ?",
                                           "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                           "}")
                                       )),
                                       class = "display"))
      })

      proxy = DT::dataTableProxy('mytable')

      shiny::observeEvent(input$mytable_rows_selected, ignoreNULL = FALSE, {

        shiny::req(!is.na(predictions()))

        curr_sel = predictions() %>%
          dplyr::filter(id %in% !!selected_rows(),
                        sample == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id)

        deselected = setdiff(curr_sel, input$mytable_rows_selected)
        if (length(deselected) > 0)
          selected_rows(setdiff(selected_rows(), deselected))

        shiny::req(!is.null(input$mytable_rows_selected))
        selected_rows(as.numeric(unique(c(selected_rows(), input$mytable_rows_selected))))
      })

      shiny::observeEvent(selected_rows(),{
        proxy %>%
          DT::selectRows(as.numeric(selected_rows()))
      })

      output$aoc_bliss = ggiraph::renderggiraph({

        shiny::req(!is.null(predictions()),
                   !is.null(input$aoc_bliss_sample),
                   input$aoc_bliss_sample %in% predictions()$sample)

        plt_df = predictions() %>%
          dplyr::filter(!is.na(AOC), !is.na(Bliss)) %>%
          dplyr::filter(sample == input$aoc_bliss_sample) %>%
          dplyr::mutate(tool_text = glue::glue("MoA1:\t\t{`MoA 1`}",
                                               "MoA2:\t\t{`MoA 2`}",
                                               "Treatment1:\t{`Treatment 1`}",
                                               "Treatment2:\t{`Treatment 2`}",
                                               .sep = "\n"))
        plt = ggplot2::ggplot(plt_df,
                              ggplot2::aes(x=AOC,
                                           y=Bliss)) +
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
                        sample == !!input$aoc_bliss_sample) %>%
          dplyr::pull(id) %>%
          as.character()

        girafe(ggobj = plt,
               options = list(opts_selection(type = "multiple",
                                             selected = curr_sel)))
      })

    })
}
