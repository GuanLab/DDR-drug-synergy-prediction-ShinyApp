
#' @title 
#' UI component for displaying an interactive plot that is clickable
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' server function.
#' @param x_breaks_name A character for the label of an input field for selecting
#' above the breaks of the x axis of the plot
#' @param info_text Small description of the plot and its content
#' @param info_title Small title of the plot (wrapped in \link[shiny]{h2})
#' 
nested_plot_ui = function(id, x_breaks_name="Breaks", info_text=shiny::HTML(""),
                          info_title=shiny::h2("")) {
  ns = shiny::NS(id)

  shiny::column(
    width = 12,
    shiny::div(class="collapse_subentry",

        shiny::div(shinyjs::hidden(shiny::actionButton(ns("prev"), "",
                                                       style="margin-bottom: 10px; padding-top: 5px; padding-bottom: 5px; z-index: 100; position: absolute;",
                                                       icon = shiny::icon("arrow-circle-left"))),
            shiny::span(shinyWidgets::pickerInput(ns("x_breaks"), x_breaks_name,
                                                           multiple = TRUE,
                                                           choices = character(0),
                                                           selected = character(0),
                                           options = shinyWidgets::pickerOptions(dropdownAlignRight=TRUE)),
                 style="float: right")),
        shiny::div(style="height: 85px; width: 20; z-index: -1; position: relative; margin-bottom: 10px;"),
        shinyjqui::jqui_resizable(
          ggiraph::ggiraphOutput(ns("categories_plt"), height = "450px", width = "100%")
        ),
        info_title,
        info_text,
        shiny::div(style="margin-bottom: 40px;")),
    collapsible_box_ui(ns("collapse"), DT::DTOutput(ns("categories_tab")))
    )
}

#' @title 
#' Server component for displaying an interactive plot that is clickable
#' 
#' @param id An ID string that corresponds with the ID used to call the module's
#' UI function.
#' @param top_level_name A character naming the name of the x cvategories for
#' the initial plot (used as x title and for the previous button)
#' @param plot_title A optional character for the plot title
#' @param summary An shiny reactive with the details for the plot output
#' @param ignore_click A character vector of x values that should be ignored
#' when clicked on
#'
nested_plot_server = function(id, top_level_name, plot_title=NULL, summary,
                              ignore_click  = c("")) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      clicked_category = shiny::reactiveVal()

      collapsible_box_server("collapse", "Show data")

      last_click = shiny::reactiveVal(top_level_name)
      last_last_click = shiny::reactiveVal()
      level = shiny::reactiveVal(0)
      last_clicked_category = shiny::reactiveVal()
      last_x = shiny::reactiveVal()
      last_last_x = shiny::reactiveVal()
      last_summary = shiny::reactiveVal()
      last_x_variable = shiny::reactiveVal()

      shiny::observeEvent(
        input$categories_plt_selected, {


          x_var = summary()$x

          shiny::req(level() < 2)
          res = NULL

          if (!is.null(input$categories_plt_selected))
            res = summary()$df %>%
              dplyr::filter(!!rlang::sym(x_var) == input$categories_plt_selected)

          ignore_selects = intersect(unique(res[[x_var]]),
                                     ignore_click)
          shiny::req(length(ignore_selects) == 0)

          level(level() + 1)
          if (level() == 1) {
            button_label = top_level_name#"Top Level"
          } else {
            last_clicked_category(last_summary())
            last_last_click(last_click())
            last_last_x(last_x())
            button_label = last_click()
          }
          last_x_variable(last_x())

          shiny::updateActionButton(session = session, inputId = "prev",
                                    label = button_label)
          last_summary(summary()$df)
          last_click(input$categories_plt_selected)
          last_x(x_var)

          clicked_category(unique(as.character(res[[x_var]])))
        })

      shiny::observeEvent(level(), {
        if (level() > 0) {
          shinyjs::show("prev")
        } else {
          shinyjs::hide("prev")
        }
      })

      shiny::observeEvent(input$prev, {
        x_var = last_x_variable()
        level(level() - 1)
        if (level() == 0) {
          res = NULL
          last_click(top_level_name)
          last_x("category")
        } else if (level() == 1) {
          res = last_clicked_category() %>%
            dplyr::filter(!!rlang::sym(x_var) == !!last_last_click())

          last_click(last_last_click())
          last_x(last_last_x())
          last_summary(res)
          res = as.character(res[[x_var]])
          shiny::updateActionButton(session = session, inputId = "prev",
                                    label = top_level_name)#"Top Level")
        } else {
          res = last_clicked_category() %>%
            dplyr::filter(!!rlang::sym(x_var) == !!last_last_click()) %>%
            dplyr::pull(!!x_var)

          last_click(last_last_click())
          last_last_x(last_x())
          shiny::updateActionButton(session = session, inputId = "prev",
                                    label = last_last_click())
        }

        clicked_category(unique(res))
      })


      shiny::observeEvent(summary()$df, {
        shiny::req(!is.null(summary()))

        df        = summary()$df
        x_var     = summary()$x
        y_var     = summary()$y

        choices = df[[x_var]] %>%
          unique() %>%
          as.character()

        selected = choices

        if(length(choices) > 10) {
          selected = df %>%
            dplyr::group_by(!!rlang::sym(x_var))
          if (!is.null(summary()$fill) | !is.null(summary()$color)) {
            selected = selected %>%
              dplyr::summarise(stat = sum(!!rlang::sym(y_var)))
          } else if ("middle" %in% colnames(selected)) {
            selected = selected %>%
              dplyr::ungroup() %>%
              dplyr::rename(stat = middle)
          } else {
            selected = selected %>%
              dplyr::summarise(stat = median(!!rlang::sym(y_var)))
          }
           selected = selected %>%
            dplyr::slice_max(stat, n = 10) %>%
            dplyr::pull(x_var) %>%
            as.character()
        }

        shinyWidgets::updatePickerInput(session,
                                        "x_breaks",
                                        choices  = choices,
                                        selected = selected)
      })

      output$categories_plt = ggiraph::renderggiraph({
        shiny::validate(
          shiny::need(!is.null(summary()), "Loading data"))

        df        = summary()$df
        x_var     = summary()$x
        y_var     = summary()$y
        fill_var  = summary()$fill
        color_var = summary()$color
        y_title   = summary()$y_title
        x_title   = summary()$x_title
        fill_title = summary()$legend_title

        castable = TRUE
        if (!is.character(df[[x_var]]))
          castable = all(df[[x_var]] == as.integer(df[[x_var]]), na.rm = TRUE)

        shiny::req(is.null(input$x_breaks) | any(input$x_breaks %in% df[[x_var]]))

        if (!is.null(input$x_breaks))
          df = df %>%
            dplyr::filter(!!rlang::sym(x_var) %in% input$x_breaks)


        if(is.null(y_var) & all(c("middle", "max", "min", "lower", "upper") %in%
                                colnames(df))) {


          plt = ggplot(df, aes_string(x=x_var, ymin="min", ymax="max",
                         middle="middle", upper="upper",
                         lower="lower"))+
            ggiraph::geom_boxplot_interactive(stat="identity",
                                              aes(tooltip=feature,
                                                  data_id=feature)) +
            ggplot2::theme_bw(base_size = 20) +
            ggplot2::scale_x_discrete(guide = guide_axis(n.dodge = 2))

        } else if ((is.character(df[[x_var]]) | castable) &
            is.double(df[[y_var]])) {
          stat = df %>%
            dplyr::group_by(!!rlang::sym(x_var)) %>%
            dplyr::summarise(med = median(!!rlang::sym(y_var))) %>%
            dplyr::arrange(-med) %>%
            utils::head(n=10)

          df = df %>%
            dplyr::semi_join(stat, by=x_var) %>%
            dplyr::ungroup() %>%
            dplyr::select(!!x_var, !!y_var) %>%
            dplyr::mutate(!!x_var := factor(!!rlang::sym(x_var),
                                            levels=stat[[x_var]]))

          plt = ggplot2::ggplot(df, ggplot2::aes_string(x=x_var, y=y_var)) +
            ggiraph::geom_boxplot_interactive(ggplot2::aes_string(tooltip = x_var,
                                                                  data_id = x_var)) +
            ggplot2::labs(title = plot_title) +
            ggplot2::theme_bw(base_size = 20)

          if (!is.null(y_title))
            plt = plt + ggplot2::ylab(y_title)

          if (!is.null(x_title))
            plt = plt + ggplot2::xlab(x_title)

          sub_label = function(x) stringr::str_sub(x, 1, 27)
          max_label_length = max(nchar(levels(df[[x_var]])))
          if (max_label_length > 15) {
            values = sub_label(levels(df[[x_var]]))
            names(values) = levels(df[[x_var]])
          }

          if (max_label_length > 15) {
            plt = plt +
              ggplot2::scale_x_discrete(labels=values) +
              ggplot2::coord_flip()
          } else if (max_label_length > 8) {
            plt = plt +
              ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
          } else {
          }
         #}  else if (is.numeric(df[[x_var]]) & is.numeric(df[[y_var]])) {
         #
         # plt = df %>%
         #   ggplot2::ggplot(aes_string(x=x_var, y=y_var, color=color_var)) +
         #   ggiraph::geom_point_interactive(ggplot2::aes(tooltip = x_var,
         #                                                data_id = x_var)) +
         #   ggplot2::labs(title = plot_title) +
         #   ggplot2::theme_bw(base_size = 20)
         # 
         # if (!is.null(y_title))
         #   plt = plt + ggplot2::ylab(y_title)
         # 
         # if (!is.null(x_title))
         #   plt = plt + ggplot2::xlab(x_title)
        } else if ((is.character(df[[x_var]]) | is.factor(df[[x_var]])) &
                    (is.integer(df[[y_var]]) | all(df[[y_var]] == as.integer(df[[y_var]])))) {

          plt = df %>%
            ggplot2::ggplot() +
            ggiraph::geom_col_interactive(ggplot2::aes_string(x=x_var, y=y_var,
                                                              data_id=x_var,
                                                              tooltip=x_var,
                                                              fill=fill_var)) +
            ggplot2::theme_bw(base_size=20) +
            viridis::scale_fill_viridis(discrete = TRUE) +
            ggplot2::guides(fill=guide_legend(fill_title))

          sub_label = function(x) stringr::str_sub(x, 1, 27)
          max_label_length = max(nchar(levels(df[[x_var]])))
          if (max_label_length > 15) {
            values = sub_label(levels(df[[x_var]]))
            names(values) = levels(df[[x_var]])
          }

          if (max_label_length > 15) {
            plt = plt +
              ggplot2::scale_x_discrete(labels=values) +
              ggplot2::coord_flip()
          } else if (max_label_length > 8) {
            plt = plt +
              ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
          } else {
          }

          if (max(df[[y_var]]) - min(df[[y_var]]) > 100)
            plt = plt +
              ggplot2::scale_y_log10()

          if (!is.null(y_title))
            plt = plt + ggplot2::ylab(y_title)

          if (!is.null(x_title))
            plt = plt + ggplot2::xlab(x_title)
        }

        width = input$categories_plt__shinyjquiBookmarkState__resizable$width / 72
        if (base::isTRUE(length(width) == 0))
          width = 9.11
        height = input$categories_plt__shinyjquiBookmarkState__resizable$height / 72
        if (base::isTRUE(length(height) == 0))
          height = 4.8

        x = girafe(ggobj = plt,
                   width_svg = width,
                   height_svg = height)

        x = girafe_options(x,
                           opts_selection(type = "single"))
      })

      output$categories_tab = DT::renderDT({

        df        = summary()$df
        x_var     = summary()$x
        y_var     = summary()$y
        fill_var  = summary()$fill
        color_var = summary()$color
        y_title   = summary()$y_title
        x_title   = summary()$x_title
        fill_title = summary()$legend_title

        if (!is.null(y_title)) {
          df = df %>%
            dplyr::rename(!!y_title := !!y_var)
          y_var = y_title
        }

        if (!is.null(fill_title)) {
          df = df %>%
            dplyr::rename(!!fill_title := !!fill_var)
          fill_var = fill_title
        }

        if (!is.null(last_click())) {
          df = df %>%
            dplyr::rename(!!last_click() := !!x_var)
          x_var = last_click()
        }

        if (is.character(df[[x_var]])) {
          df = df %>%
            dplyr::select(!!x_var, !!y_var) %>%
            dplyr::mutate(dplyr::across(tidyselect:::where(is.double),
                                        ~round(.x, 2)))
        } else {
          df = df %>%
            dplyr::filter(!is.na(!!rlang::sym(y_var)))
        }

        DT::datatable(df, rownames = FALSE)
      })

      return(clicked_category)
    })
}
