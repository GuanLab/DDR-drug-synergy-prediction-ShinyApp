#' @title 
#' Wrapper for indexed fst access
#' 
#' @description 
#' This method uses an fst file  (file.fst) and a corresponding index file 
#' (file_index.fst) to optimize the access of fst files. The index file
#' should contain a column with the name of the column that is indexed in the 
#' fst file and to columns (start, end) that list the first and last value of
#' an value in the column. The fst file must be sorted by the column.
#' 
#' @param fst_name A character with the fst file name (optionally the same 
#' location contains a corresponding index file).
#' @param filter A named list with column name, value pairs (used for filtering)
#' @param distinct A logical indicating if only unique values should be returned
#' @param order A character vector, listing columns by that the returned tibble
#' should be sorted
#' @param stringAsFactors A logical, if characters should be cast to facttors
#' @param peek A logical, if only the top 10 rows should be returned
#' @param folder A character with the path to the fst file
#' @param cols A character vector of columns to be loaded
#' 
#' @return A tibble listing entries for the defined parameters
#' 
load_fst = function(fst_name, cols = NULL, filter = NULL, distinct = FALSE, 
                    order = NULL, stringAsFactors = FALSE, peek = FALSE,
                    folder=NULL) {
  
  if (is.null(fst_name)) 
    stop("No fst file selected")
  
  if (!is.character(fst_name)) 
    stop("'fst_name' must be a character value. Please check your parameters.")
  
  if (is.null(folder))
    folder = getOption("fst_folder")
  
  fst_file = paste(folder, fst_name, sep="/")
  
  if (!file.exists(fst_file))
    stop("FST file unavailable. Please check the view name and folder or ",
         "options('fst_folder') variables")
  
  q = fst::fst(fst_file)
  if (is.null(cols)) 
    cols = colnames(q[1, ])
  
  # check if filtering is defined correctly
  filter_columns = names(filter) %in% colnames(q[1, ])
  if (any(!filter_columns)) 
    stop(sprintf("Filter columns not available: %s", 
                 paste(names(filter)[!filter_columns],  collapse = ", ")))
  
  idx = 1:nrow(q)
  col_names = list()
  idx_file = gsub("\\.fst", "_index.fst", fst_file)
  if (file.exists(idx_file)) {
    q_idx = fst::fst(idx_file)
    
    col_names = filter[intersect(names(filter), colnames(q_idx))]
    filter = filter[setdiff(names(filter), colnames(q_idx))]
    index_idx = rep(T, nrow(q_idx))
    if (length(col_names) > 0) {
      for (i in 1:length(col_names)) {
        if (length(col_names[[i]]) == 1) {
          index_idx = index_idx & (q_idx[[names(col_names)[i]]] == 
                                     col_names[[i]])
        } else {
          index_idx = index_idx & (q_idx[[names(col_names)[i]]] %in% 
                                     col_names[[i]])
        }
      }
      if (isTRUE(any(index_idx))) {
        idx_seqs = plyr::alply(q_idx[index_idx, ], 1, 
                               function(x) seq(x[["start"]], x[["end"]]))
        index_idx = idx_seqs %>% unlist %>% unname %>% 
          unique
        idx = index_idx
      } else {
        idx = c()
      }
    }
  }
  
  if (!is.null(filter)) {
    if (length(filter) > 0) {
      if (length(col_names) > 0) {
        consecutives = R.utils::seqToIntervals(idx)
      } else {
        consecutives = matrix(c(1, nrow(q)), nrow = 1, 
                              ncol = 2)
      }
      res_idx = list()
      for (j in 1:nrow(consecutives)) {
        rg = consecutives[j, 1]:consecutives[j, 2]
        for (i in 1:length(filter)) {
          que = q[rg, names(filter)[i]]
          if (is.call(filter[[i]])) {
            stop("TODO")
          } else if (length(filter[[i]]) == 1) {
            test = que == filter[[i]] & !is.na(que)
          } else {
            test = que %in% filter[[i]]
          }
          rg = rg[test]
        }
        res_idx[[j]] = rg
      }
      idx = unlist(res_idx)
    }
  }
  
  if (peek) {
    idx = idx %>% .[1:10]
  }
  
  if (length(idx) == 0) {
    dt <- tibble::tibble(matrix(nrow = 0, ncol = length(cols)))
    dt <- setNames(dt, cols)
  } else {
    consecutives = R.utils::seqToIntervals(idx) %>% as.data.frame() %>% 
      dplyr::mutate(rows = to - from)
    if (any(consecutives$rows > 10) & nrow(consecutives) < 
        20) {
      dt_l = list()
      for (i in 1:nrow(consecutives)) {
        idx_i_start = consecutives[i, 1]
        idx_i_end = consecutives[i, 2]
        idx_i = idx_i_start:idx_i_end
        dt_l[[i]] = q[idx_i, cols, drop = F]
      }
      dt = dt_l %>% dplyr::bind_rows()
    }
    else {
      dt = q[idx, cols, drop = F]
    }
  }
  if (stringAsFactors) {
    dt = dt %>% dplyr::mutate_if(is.character, as.factor)
  }
  else {
    dt = dt %>% dplyr::mutate_if(is.factor, as.character)
  }
  if (is.data.frame(dt)) {
    dt = dt %>% dplyr::as_tibble()
  }
  
  if (distinct) 
    dt = dt %>% dplyr::distinct()
  if (!is.null(order)) 
    if (all(order %in% colnames(dt))) 
      dt = dt %>% dplyr::arrange(!!!rlang::syms(order))
  
  dt
}