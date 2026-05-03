# dichotomize: Transform character and factor objects into dichotomies.
## dichotomize ----
dichotomize <- function(data, variables, sep = "", regex = FALSE, 
                        text = c("none", "lower", "upper"),
                        min = 1, length = 0, sort = TRUE,
                        values = NULL, stopwords = NULL, nas = "None",
                        sparse = FALSE, add = FALSE) {
  
  text <- match.arg(text)
  
  if (!is.data.frame(data)) stop("Se requiere un data.frame")
  if (min > 0 && min < 1) min <- min * nrow(data)
  
  # ====  SEP="C" Rule ====.
  if (!is.null(nas) && sep[1] == "C") {
    oldData <- data
    cols_to_paste <- data[, variables, drop = FALSE]
    sep_char <- "|"
    String <- apply(cols_to_paste, 1, function(row) {
      valid_vals <- na.omit(row)
      if (base::length(valid_vals) == 0) return(nas)
      paste(valid_vals, collapse = sep_char)
    })
    data <- data.frame(String = String, stringsAsFactors = FALSE)
    variables <- "String"
    sep <- sep_char
  }
  
  if (base::length(sep) != base::length(variables)) {
    sep <- rep(sep[1], base::length(variables))
  }
  names(sep) <- variables
  
  # 1. Text extraction and cleaning
  list_of_lists <- list()
  for (col_name in variables) {
    current_sep <- sep[col_name]
    
    raw_values <- as.character(data[[col_name]])
    raw_values[is.na(data[[col_name]])] <- ""
    is_blank <- trimws(raw_values) == ""
    
    if (current_sep != "" && current_sep != "C") {
      parts <- strsplit(raw_values, current_sep, fixed = !regex)
    } else {
      parts <- as.list(raw_values)
    }
    
    if (text %in% c("lower", "upper")) {
      if (!is.null(stopwords)) {
        if (text == "lower") stopwords <- tolower(stopwords)
        else if (text == "upper") stopwords <- toupper(stopwords)
      }
      parts <- lapply(parts, function(x) {
        if (text == "lower") x <- tolower(x)
        else if (text == "upper") x <- toupper(x)
        
        x <- gsub("[^[:alpha:][:space:]_]", "", x)
        x <- trimws(x)
        x <- x[x != ""]
        
        if (!is.null(stopwords)) x <- x[!(x %in% stopwords)]
        return(x)
      })
    }
    
    parts[is_blank] <- list(character(0))
    list_of_lists[[col_name]] <- parts
  }
  
  # ==== POOLING ====.
  pool <- regex && base::length(variables) > 1
  
  if (pool) {
    split_list <- do.call(Map, c(list(c), list_of_lists))
    split_list <- lapply(split_list, unique)
    pool_name <- variables[1] 
    processing_groups <- list()
    processing_groups[[pool_name]] <- split_list
  } else {
    processing_groups <- lapply(list_of_lists, function(sl) lapply(sl, unique))
  }
  
  # Solo aplica el prefijo si hay más de 1 grupo independiente
  apply_category_prefix <- base::length(processing_groups) > 1
  
  result_list <- list()
  
  for (group_name in names(processing_groups)) {
    split_list <- processing_groups[[group_name]]
    all_terms <- unlist(split_list)
    
    if (base::length(all_terms) == 0) next
    
    # 2. Filtering and cuttering
    if (is.null(values)) {
      term_counts <- table(all_terms)
      if (sort) term_counts <- sort(term_counts, decreasing = TRUE)
      if (min > 0) term_counts <- term_counts[term_counts >= min]
      if (length > 0 && base::length(term_counts) > length) {
        term_counts <- term_counts[1:length] 
      }
      valid_terms <- names(term_counts)
    } else {
      valid_terms <- as.character(values)
    }
    
    if (base::length(valid_terms) == 0) next
    
    # 3. Matrix mapping
    row_indices <- rep(seq_along(split_list), lengths(split_list))
    col_match <- match(all_terms, valid_terms)
    keep <- !is.na(col_match)
    
    if (!is.null(values) && (min > 0 || length > 0)) {
      freqs <- table(col_match[keep])
      valid_idx <- as.integer(names(freqs))
      
      if (min > 0) valid_idx <- valid_idx[freqs >= min]
      if (sort) valid_idx <- valid_idx[order(-freqs[match(valid_idx, names(freqs))])]
      if (length > 0 && base::length(valid_idx) > length) {
        valid_idx <- valid_idx[1:length]
      }
      
      keep <- keep & (col_match %in% valid_idx)
      valid_terms <- valid_terms[valid_idx]
      col_match <- match(all_terms, valid_terms)
      keep <- !is.na(col_match)
    }
    
    if (base::length(valid_terms) == 0) next
    
    # ==== PREFIXES and CATEGORIES ====.
    if (apply_category_prefix) {
      col_names_formatted <- paste(group_name, valid_terms, sep = ":")
    } else {
      col_names_formatted <- valid_terms
    }
    
    # 4. Matrix creation
    if (sparse) {
      if (!requireNamespace("Matrix", quietly = TRUE)) stop("Falta 'Matrix'.")
      Q <- Matrix::sparseMatrix(
        i = row_indices[keep], j = col_match[keep], x = 1,
        dims = c(nrow(data), base::length(valid_terms))
      )
      colnames(Q) <- col_names_formatted
      result_list[[group_name]] <- Q
      
    } else {
      Q <- matrix(0L, nrow = nrow(data), ncol = base::length(valid_terms))
      colnames(Q) <- col_names_formatted
      Q[cbind(row_indices[keep], col_match[keep])] <- 1L
      
      W <- as.data.frame(Q, check.names = FALSE)
      
      if (!is.null(nas)) {
        nas_name <- paste(group_name, nas, sep = ":")
        if (!(nas %in% valid_terms) || nas_name != nas) {
          W[[nas_name]] <- as.integer(rowSums(Q) == 0)
        }
      }
      result_list[[group_name]] <- W
    }
  }
  
  if (base::length(result_list) == 0) return(if(exists("oldData")) oldData else data)
  
  # list without name
  names(result_list) <- NULL 
  final_res <- do.call(cbind, result_list)
  
  # 5. Data return 
  if (add && !sparse) {
    base_data <- if (exists("oldData")) oldData else data
    return(cbind(base_data, final_res))
  } else {
    return(final_res)
  }
}
