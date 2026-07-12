# Mon Jan  5 20:04:05 2026 ------------------------------
# Modesto Escobar + copilot + Gemini 

# ============================================================ #
# GW contrasts and Nodes/Links summary (R Base + marginaleffects)
# Unified with: auto-quasi, robust weights, binomial labels,
# Nodes sanitization and scaling to % in binomial (GW + slopes).
# MULTINOMIAL SUPPORT ADDED (nnet::multinom)
# ============================================================ #

# ------------------------ Utilities ------------------------- #

# Is it a binomial family?
.is_binomial_family <- function(fam) {
  if (is.null(fam)) return(FALSE)
  if (is.character(fam)) return(tolower(trimws(fam)) == "binomial")
  if (is.list(fam) && !is.null(fam$family)) return(grepl("binomial", fam$family, fixed = TRUE))
  FALSE
}

# Positive label according to type of response
.positive_label_from_y <- function(y) {
  if (is.factor(y)) {
    if (nlevels(y) == 2L) return(levels(y)[2L]) # event=1 is second level
    return(NA_character_)
  } else if (is.logical(y)) {
    return("TRUE")
  } else if (is.numeric(y)) {
    u <- sort(unique(na.omit(y)))
    if (length(u) == 2L && all(u %in% c(0, 1))) return("1")
    return(NA_character_)
  } else if (is.matrix(y) && ncol(y) == 2L) {
    nm <- colnames(y)
    if (!is.null(nm) && nzchar(nm[1L])) return(nm[1L])
    return("success")
  }
  NA_character_
}

# Concatenate response label + positive
.cat_label <- function(target, positive, sep = ":") {
  if (is.na(positive) || !nzchar(positive)) return(target)
  paste0(target, sep, positive)
}

# rbind with column joining
.rbind_union <- function(dfs) {
  dfs <- Filter(function(z) is.data.frame(z) && nrow(z) > 0, dfs)
  if (length(dfs) == 0) return(data.frame())
  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs2 <- lapply(dfs, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss)) for (mc in miss) df[[mc]] <- NA
    df[all_cols]
  })
  out <- do.call(rbind, dfs2)
  rownames(out) <- NULL
  out
}

# Key Deduction for Nodes
.pick_key <- function(nodes, .key = NULL) {
  if (!is.null(.key)) {
    if (!.key %in% names(nodes)) stop("The `.key` key is not in `Nodes`.")
    return(.key)
  }
  candidates <- c("name","Name","id","ID","Node","node")
  hit <- candidates[candidates %in% names(nodes)]
  if (length(hit) >= 1) return(hit[1L])
  names(nodes)[1L]
}

# Merging/deduplication of Nodes by key
.merge_nodes <- function(existing, incoming, .key = NULL,
                         .mode = c("keep_first","update","sum_numeric"),
                         .sum_cols = NULL) {
  .mode <- match.arg(.mode)
  if (is.null(incoming) || !is.data.frame(incoming) || nrow(incoming) == 0) return(existing)
  if (is.null(existing)) return(incoming)
  key <- .pick_key(incoming, .key)
  all_cols <- unique(c(names(existing), names(incoming)))
  for (mc in setdiff(all_cols, names(existing))) existing[[mc]] <- NA
  for (mc in setdiff(all_cols, names(incoming))) incoming[[mc]] <- NA
  existing <- existing[all_cols]
  incoming <- incoming[all_cols]
  ex_keys <- existing[[key]]
  in_keys <- incoming[[key]]
  comunes <- intersect(ex_keys, in_keys)
  nuevas  <- setdiff(in_keys, ex_keys)
  out <- existing
  if (length(nuevas)) {
    add_rows <- incoming[incoming[[key]] %in% nuevas, , drop = FALSE]
    out <- rbind(out, add_rows)
  }
  if (length(comunes)) {
    for (k in comunes) {
      i_ex <- which(out[[key]] == k)[1L]
      i_in <- which(incoming[[key]] == k)[1L]
      for (col in all_cols) {
        if (col == key) next
        ex_val <- out[i_ex, col]
        in_val <- incoming[i_in, col]
        if (.mode == "keep_first") {
          next
        } else if (.mode == "update") {
          out[i_ex, col] <- if (!is.na(in_val)) in_val else ex_val
        } else { # sum_numeric
          if (is.numeric(ex_val) && is.numeric(in_val) &&
              (is.null(.sum_cols) || col %in% .sum_cols)) {
            out[i_ex, col] <- ifelse(is.na(ex_val), 0, ex_val) + ifelse(is.na(in_val), 0, in_val)
          } else {
            out[i_ex, col] <- if (!is.na(in_val)) in_val else ex_val
          }
        }
      }
    }
  }
  rownames(out) <- NULL
  out
}

# --------------------------- Specifications parser --------------------------- #

# Detectors
.is_model <- function(x) inherits(x, c("lm","glm","mlm","multinom"))
.is_list_of_models <- function(x) {
  is.list(x) && !.is_model(x) && length(x) > 0L && all(vapply(x, .is_model, logical(1L)))
}

# Parser top-level (per line / comma top level)
.split_top_level <- function(s, seps = c("\n", ",")) {
  if (!nzchar(s)) return(list(lhs = "", rhs = NULL, sep = NULL))
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  depth_par <- 0L; depth_brk <- 0L; depth_crl <- 0L
  in_squote <- FALSE; in_dquote <- FALSE
  esc <- FALSE; pos <- NA_integer_; sep_found <- NULL
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (esc) { esc <- FALSE; next }
    if (ch == "\\") { esc <- TRUE; next }
    if (!in_squote && !in_dquote) {
      if (ch == "(") depth_par <- depth_par + 1L
      else if (ch == ")") depth_par <- max(0L, depth_par - 1L)
      else if (ch == "[") depth_brk <- depth_brk + 1L
      else if (ch == "]") depth_brk <- max(0L, depth_brk - 1L)
      else if (ch == "{") depth_crl <- depth_crl + 1L
      else if (ch == "}") depth_crl <- max(0L, depth_crl - 1L)
      else if (depth_par == 0L && depth_brk == 0L && depth_crl == 0L && ch %in% seps) {
        pos <- i; sep_found <- ch; break
      }
    }
    if (!in_dquote && ch == "'") in_squote <- !in_squote
    else if (!in_squote && ch == "\"") in_dquote <- !in_dquote
  }
  if (is.na(pos)) {
    list(lhs = trimws(s), rhs = NULL, sep = NULL)
  } else {
    lhs <- trimws(paste0(chars[seq_len(pos - 1L)], collapse = ""))
    rhs <- trimws(paste0(chars[seq(pos + 1L, length(chars))], collapse = ""))
    list(lhs = lhs, rhs = if (nzchar(rhs)) rhs else NULL, sep = sep_found)
  }
}

.safe_eval_family <- function(fam_str) {
  if (!nzchar(fam_str)) return(NULL)
  fam_str <- sub("^family\\s*=\\s*", "", fam_str)
  
  # 1. Factor out the base BEFORE adding the parentheses
  base <- sub("\\s*\\(.*$", "", fam_str)
  
  # 2. Special support for “multinomial” strings (we're nipping that in the bud here)
  if (base == "multinomial") return("multinomial")
  
  # 3. If it's not a multinomial, we make sure it has parentheses
  if (!grepl("\\(", fam_str)) fam_str <- paste0(fam_str, "()")
  
  allowed <- c("gaussian","binomial","poisson","Gamma","inverse.gaussian",
               "quasi","quasibinomial","quasipoisson")
  
  if (!(base %in% allowed)) {
    stop(sprintf("Family not allowed: '%s'. Allowed: %s, multinomial", base, paste(allowed, collapse = ", ")))
  }
  
  env <- new.env(parent = emptyenv())
  env$gaussian <- stats::gaussian
  env$binomial <- stats::binomial
  env$poisson <- stats::poisson
  env$Gamma <- stats::Gamma
  env$`inverse.gaussian` <- stats::`inverse.gaussian`
  env$quasi <- stats::quasi
  env$quasibinomial <- stats::quasibinomial
  env$quasipoisson <- stats::quasipoisson
  call <- try(str2lang(fam_str), silent = TRUE)
  if (inherits(call, "try-error")) stop(sprintf("The family could not be parsed: %s", fam_str))
  fam <- try(eval(call, envir = env), silent = TRUE)
  if (inherits(fam, "try-error")) stop(sprintf("The family could not be evaluated: %s", fam_str))
  fam
}
.parse_specs <- function(spec) {
  lines <- if (length(spec) == 1L) unlist(strsplit(spec, "\n")) else spec
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!grepl("^\\s*#", lines)]
  out <- vector("list", length(lines))
  names(out) <- lines
  for (i in seq_along(lines)) {
    sp <- .split_top_level(lines[i], seps = c("\n", ",")) # separador '\n' o coma top-level
    fm_text <- sp$lhs
    fam_obj <- NULL
    if (!is.null(sp$rhs)) fam_obj <- .safe_eval_family(sp$rhs)
    fm <- try(stats::as.formula(fm_text), silent = TRUE)
    if (inherits(fm, "try-error")) stop(sprintf("It could not be converted into a formula: '%s'", fm_text))
    out[[i]] <- list(formula = fm, family = fam_obj)
  }
  out
}

# --------------------------- Adjustment from specifications --------------------------- #
.fit_models_from_specs <- function(spec, data, default_family = gaussian(), weights = NULL, auto_quasi = TRUE) {
  entries <- .parse_specs(spec)
  if (is.null(data)) stop("You must pass `data=` to adjust models from specifications.")
  
  # Prepare weight vector (optional)
  w <- NULL
  if (!is.null(weights)) {
    if (is.numeric(weights) && length(weights) == nrow(data)) {
      w <- weights
    } else {
      stop("`weights` must be numeric and of length nrow(data) when adjusting from specifications.")
    }
  }
  
  # Copy of data + column .wg__ if there are vector weights
  if (!is.null(w)) {
    data_w <- data
    .wg__  <- w
    data_w[[".wg__"]] <- .wg__
  } else {
    data_w <- data
  }
  
  # Default family (with auto-quasi if binomial)
  default_family_use <- default_family
  if (isTRUE(auto_quasi) && .is_binomial_family(default_family)) {
    default_family_use <- stats::quasibinomial()
  }
  
  mods <- lapply(entries, function(e) {
    fam_use <- e$family
    f_use <- e$formula
    # THE MAGIC TRICK! We tell the formula to remember this environment 
    # so that marginaleffects doesn't get lost looking for the data later.
    environment(f_use) <- environment()
    
    # If the family is in the spec and is binomial, apply auto-quasi
    if (isTRUE(auto_quasi) && .is_binomial_family(fam_use)) {
      fam_use <- stats::quasibinomial()
    }
    
    # If family comes as a string, map to object
    if (is.character(fam_use)) {
      fam_str <- tolower(trimws(fam_use))
      if (fam_str == "binomial") {
        fam_use <- if (isTRUE(auto_quasi)) stats::quasibinomial() else stats::binomial()
      } else if (fam_str == "multinomial") {
        fam_use <- "multinomial"
      } else if (fam_str == "gaussian") {
        fam_use <- stats::gaussian()
      } else if (fam_str == "poisson") {
        fam_use <- stats::poisson()
      } else if (fam_str == "gamma") {
        fam_use <- stats::Gamma()
      }
    }
    
    # Model adjustment (We use eval(bquote(...)) to destroy the reference to the variable 'e')
    if (identical(fam_use, "multinomial")) {
      if (!requireNamespace("nnet", quietly = TRUE)) stop("Install 'nnet' for multinomial models.")
      if (is.null(w)) eval(bquote(nnet::multinom(formula = .(f_use), data = data_w, trace = FALSE)))
      else            eval(bquote(nnet::multinom(formula = .(f_use), data = data_w, weights = .wg__, trace = FALSE)))
    } else if (is.null(fam_use)) {
      if (identical(default_family_use, gaussian())) {
        if (is.null(w)) eval(bquote(stats::lm(formula = .(f_use), data = data_w)))
        else            eval(bquote(stats::lm(formula = .(f_use), data = data_w, weights = .wg__)))
      } else {
        if (is.null(w)) eval(bquote(stats::glm(formula = .(f_use), data = data_w, family = .(default_family_use))))
        else            eval(bquote(stats::glm(formula = .(f_use), data = data_w, family = .(default_family_use), weights = .wg__)))
      }
    } else {
      if (is.null(w)) eval(bquote(stats::glm(formula = .(f_use), data = data_w, family = .(fam_use))))
      else            eval(bquote(stats::glm(formula = .(f_use), data = data_w, family = .(fam_use), weights = .wg__)))
    }
  })
  
  names(mods) <- names(entries)
  mods
}

# --------------------------- GW contrasts (group) --------------------------- #
contrastes_gw <- function(model, group_var, weights = NULL, vcov = "HC1") {
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("You need the 'marginaleffects' package. Install 'marginaleffects'.")
  }
  if (length(group_var) < 1) stop("'group_var' must contain at least one variable name.")
  mf <- stats::model.frame(model)
  if (!is.data.frame(mf)) stop("Could not extract a model.frame from object 'model'.")
  
  # Validate group variables (existence, factor, >=2 levels)
  for (v in group_var) {
    if (!v %in% names(mf)) stop(sprintf("The variable '%s' is not in model.frame(model).", v))
    if (!is.factor(mf[[v]])) stop(sprintf("The variable '%s' is not a factor.", v))
    if (length(unique(mf[[v]])) < 2) stop(sprintf("The variable '%s' has only one observed level.", v))
  }
  
  # Target
  fchar <- as.character(stats::formula(model))
  target <- if (length(fchar) >= 2) fchar[2] else NA_character_
  
  # To detect binomial y multinomial
  is_binom <- inherits(model, "glm") && !is.null(model$family) && grepl("binomial", model$family$family, fixed = TRUE)
  is_multinom <- inherits(model, "multinom")
  
  # --- Weights (robust) ---
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L && weights %in% names(mf)) {
      w <- mf[[weights]]
    } else {
      w_expr <- substitute(weights)
      w_name <- if (is.symbol(w_expr) || is.name(w_expr)) deparse(w_expr) else NULL
      if (!is.null(w_name) && w_name %in% names(mf)) {
        w <- mf[[w_name]]
      } else if (is.numeric(weights) && length(weights) == nrow(mf)) {
        w <- weights
      } else if ("(weights)" %in% names(mf)) {
        w <- mf[["(weights)"]]
      } else {
        warning("'weights' argument not recognized; using w = 1 (no weighting).")
        w <- rep(1, nrow(mf))
      }
    }
  } else if ("(weights)" %in% names(mf)) {
    w <- mf[["(weights)"]]
  } else {
    w <- rep(1, nrow(mf))
  }
  
  # Average margins per level combination
  margenes <- marginaleffects::avg_predictions(
    model,
    variables = group_var,
    wts = w,
    vcov = vcov
  )
  
  # Weights added per level combination (alignment)
  df_tmp <- mf[, group_var, drop = FALSE]
  df_tmp$w <- w
  rhs <- paste(group_var, collapse = " + ")
  form_agg <- stats::as.formula(paste("w ~", rhs))
  df_pesos <- stats::aggregate(form_agg, data = df_tmp, FUN = sum)
  names(df_pesos)[names(df_pesos) == "w"] <- "peso_abs"
  df_pesos$prop <- df_pesos$peso_abs / sum(df_pesos$peso_abs)
  
  # Alignment by keys
  marg_df <- as.data.frame(margenes)
  cols_ok <- group_var[group_var %in% names(marg_df)]
  if (length(cols_ok) != length(group_var)) {
    stop("The 'margins' columns do not contain all the variables of 'group_var'.")
  }
  key_marg <- apply(marg_df[, cols_ok, drop = FALSE], 1,
                    function(z) paste(as.character(z), collapse = "\u0001"))
  key_peso <- apply(df_pesos[, group_var, drop = FALSE], 1,
                    function(z) paste(as.character(z), collapse = "\u0001"))
  indices <- match(key_marg, key_peso)
  if (any(is.na(indices))) stop("Weights could not be aligned with 'margin' levels")
  vector_pesos_fijo <- df_pesos$prop[indices]
  
  # GW Hypothesis: subtract the global weighted average (adapted for Multinomial groups)
  funcion_gw <- function(x) {
    est <- x$estimate
    res <- numeric(length(est))
    
    if ("group" %in% names(x)) {
      for (g in unique(x$group)) {
        idx <- which(x$group == g)
        w_g <- vector_pesos_fijo[idx]
        w_g <- w_g / sum(w_g) # rescale weights within the group
        gm_g <- sum(est[idx] * w_g)
        res[idx] <- est[idx] - gm_g
      }
    } else {
      w_norm <- vector_pesos_fijo / sum(vector_pesos_fijo)
      gm <- sum(est * w_norm)
      res <- est - gm
    }
    return(res)
  }
  
  res <- marginaleffects::hypotheses(margenes, hypothesis = funcion_gw)
  resultado_final <- as.data.frame(res)
  
  # --- SOLUTION: Restore the category columns that the hypotheses() function deleted ---
  for (col_cat in c("group", "class", "response")) {
    if ((col_cat %in% names(marg_df)) && !(col_cat %in% names(resultado_final))) {
      resultado_final[[col_cat]] <- marg_df[[col_cat]]
    }
  }
  # ----------------------------------------------------------------------------------------.

  # --- Contrasts vs the first (base) category: same margins, delta-method SEs
  base_lev <- vapply(group_var, function(v) as.character(levels(mf[[v]])[1L]), character(1))
  is_base_row <- rep(TRUE, nrow(marg_df))
  for (i in seq_along(group_var))
    is_base_row <- is_base_row & as.character(marg_df[[group_var[i]]]) == base_lev[i]
  funcion_base <- function(x) {
    est <- x$estimate
    res <- numeric(length(est))
    pick <- function(idx) { b <- idx[which(is_base_row[idx])[1L]]
                            if (is.na(b)) idx[1L] else b }
    if ("group" %in% names(x)) {
      for (g in unique(x$group)) {
        idx <- which(x$group == g)
        res[idx] <- est[idx] - est[pick(idx)]
      }
    } else {
      res <- est - est[pick(seq_along(est))]
    }
    return(res)
  }
  res_base <- try(marginaleffects::hypotheses(margenes, hypothesis = funcion_base), silent = TRUE)
  if (!inherits(res_base, "try-error")) {
    rb <- as.data.frame(res_base)
    resultado_final$estimate.base  <- rb$estimate
    resultado_final$std.error.base <- rb$std.error
    resultado_final$p.value.base   <- rb$p.value
    resultado_final$std.error.base[is_base_row] <- NA_real_
    resultado_final$p.value.base[is_base_row]   <- NA_real_
  }

  # --- Scaled to % if the model is binomial or multinomial ---
  if ((is_binom || is_multinom) && nrow(resultado_final) > 0) {
    pct_cols <- intersect(c("estimate","std.error","conf.low","conf.high",
                            "estimate.base","std.error.base"), names(resultado_final))
    for (cc in pct_cols) resultado_final[[cc]] <- resultado_final[[cc]] * 100
  }
  
  # Source/Target Labels
  var_combo <- paste(group_var, collapse = "*")
  niveles_mat <- marg_df[, cols_ok, drop = FALSE]
  niveles_comb <- if (ncol(niveles_mat) == 1) {
    as.character(niveles_mat[[1]])
  } else {
    apply(niveles_mat, 1, function(z) paste(as.character(z), collapse = "*"))
  }
  resultado_final$Source <- paste0(var_combo, ":", niveles_comb)
  
  # Configure Target to support multinomial models (searching for “group,” “class,” or “response”)
  if ("group" %in% names(resultado_final)) {
    resultado_final$Target <- paste0(target, ":", resultado_final$group)
  } else if ("class" %in% names(resultado_final)) {
    resultado_final$Target <- paste0(target, ":", resultado_final$class)
  } else if ("response" %in% names(resultado_final)) {
    resultado_final$Target <- paste0(target, ":", resultado_final$response)
  } else {
    resultado_final$Target <- target
  }
  
  # Cleaning
  rm_cols <- intersect(c("variable","term","contrast","predicted_lo","predicted_hi","predicted","group", "class", "response"), names(resultado_final))
  if (length(rm_cols) > 0) resultado_final <- resultado_final[, setdiff(names(resultado_final), rm_cols), drop = FALSE]
  resultado_final
}

# --------------------------- Nodes/Links summary (unique model) --------------------------- #
contrastes_gw_one <- function(model, weights = NULL, vcov = "HC1",
                              drop_empty = TRUE, stdcov = FALSE,
                              outliers = TRUE, q_low = 0.05, q_high = 0.95) {
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("You need the 'marginaleffects' package. Install 'marginaleffects'.")
  }
  mf <- stats::model.frame(model)
  if (!is.data.frame(mf)) stop("Could not extract a model.frame from object 'model'.")
  N_model <- nrow(mf)
  
  # Quantile validations
  if (!is.numeric(q_low) || !is.numeric(q_high) ||
      length(q_low) != 1 || length(q_high) != 1 ||
      is.na(q_low) || is.na(q_high) ||
      q_low < 0 || q_high > 1 || q_low >= q_high) {
    stop("Invalid arguments: requires 0 <= q_low < q_high <= 1.")
  }
  
  # Target (response) and binomial detection
  fchar  <- as.character(stats::formula(model))
  target <- if (length(fchar) >= 2) fchar[2] else NA_character_
  is_binom <- inherits(model, "glm") && !is.null(model$family) && grepl("binomial", model$family$family, fixed = TRUE)
  is_multinom <- inherits(model, "multinom")
  y <- if (!is.na(target) && target %in% names(mf)) mf[[target]] else NULL
  
  # --- Weights (robust) ---
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L && weights %in% names(mf)) {
      w <- mf[[weights]]
    } else {
      w_expr <- substitute(weights)
      w_name <- if (is.symbol(w_expr) || is.name(w_expr)) deparse(w_expr) else NULL
      if (!is.null(w_name) && w_name %in% names(mf)) {
        w <- mf[[w_name]]
      } else if (is.numeric(weights) && length(weights) == nrow(mf)) {
        w <- weights
      } else if ("(weights)" %in% names(mf)) {
        w <- mf[["(weights)"]]
      } else {
        warning("'weights' argument not recognized; using w = 1 (no weighting).")
        w <- rep(1, nrow(mf))
      }
    }
  } else if ("(weights)" %in% names(mf)) {
    w <- mf[["(weights)"]]
  } else {
    w <- rep(1, nrow(mf))
  }
  w_tot <- sum(w, na.rm = TRUE)
  
  # Weighted helpers
  wmean <- function(x, w) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
  wsd   <- function(x, w) {
    m <- wmean(x, w)
    sqrt(sum(w * (x - m)^2, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE))
  }
  
  # Extremes for normalization
  range_fun <- function(x, outliers, ql, qh) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(c(NA_real_, NA_real_))
    if (isTRUE(outliers)) {
      c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    } else {
      q <- stats::quantile(x, probs = c(ql, qh), na.rm = TRUE, names = FALSE, type = 7)
      c(q[1], q[2])
    }
  }
  
  # Determine positive category in binomial (forma robusta)
  pred_cat <- NA_character_
  if (is_binom && !is.null(y)) {
    if (is.factor(y)) {
      pred_cat <- levels(y)[2] # En glm binomial se predice el segundo nivel
    } else if (is.logical(y)) {
      pred_cat <- "TRUE"
    } else {
      # Si es 0/1, toma el segundo valor ordenado (el 1)
      pred_cat <- as.character(sort(unique(na.omit(y)))[2]) 
    }
  }
  target_label <- if (is_binom && !is.na(pred_cat)) paste0(target, ":", pred_cat) else target
  
  # Terms in the order of the model
  trm         <- stats::terms(model)
  term_labels <- attr(trm, "term.labels")
  
  # ---------------------------- #
  # 1) LINKS (GW + slopes)
  # ---------------------------- #
  res_list <- list()
  for (lab in term_labels) {
    if (grepl(":", lab, fixed = TRUE)) {
      parts <- strsplit(lab, ":", fixed = TRUE)[[1]]
      ok <- TRUE
      for (p in parts) {
        if (!(p %in% names(mf) && is.factor(mf[[p]]))) { ok <- FALSE; break }
        if (drop_empty && length(unique(mf[[p]])) < 2) { ok <- FALSE; break }
      }
      if (!ok) next
      tmp <- try(contrastes_gw(model, group_var = parts, weights = weights, vcov = vcov), silent = TRUE)
      if (!inherits(tmp, "try-error")) {
        if (is_binom && !is.na(target_label) && "Target" %in% names(tmp)) tmp$Target <- target_label
        res_list[[length(res_list) + 1]] <- tmp
      } else {
        warning(sprintf("GW interaction failed '%s': %s", lab, as.character(tmp)))
      }
    } else {
      if (!(lab %in% names(mf))) next
      if (is.factor(mf[[lab]])) {
        if (drop_empty && length(unique(mf[[lab]])) < 2) next
        tmp <- try(contrastes_gw(model, group_var = lab, weights = weights, vcov = vcov), silent = TRUE)
        if (!inherits(tmp, "try-error")) {
          if (is_binom && !is.na(target_label) && "Target" %in% names(tmp)) tmp$Target <- target_label
          res_list[[length(res_list) + 1]] <- tmp
        } else {
          warning(sprintf("GW factor failed '%s': %s", lab, as.character(tmp)))
        }
      } else if (is.numeric(mf[[lab]])) {
        cont_res <- try(marginaleffects::avg_slopes(model, variables = lab, wts = w, vcov = vcov), silent = TRUE)
        if (!inherits(cont_res, "try-error")) {
          cont_df <- as.data.frame(cont_res)
          # Optional standardization
          if (isTRUE(stdcov)) {
            s <- wsd(mf[[lab]], w)
            if ("estimate" %in% names(cont_df)) cont_df$estimate <- cont_df$estimate * s
            if ("std.error" %in% names(cont_df)) cont_df$std.error <- cont_df$std.error * s
            if ("conf.low" %in% names(cont_df)) cont_df$conf.low <- cont_df$conf.low * s
            if ("conf.high" %in% names(cont_df)) cont_df$conf.high <- cont_df$conf.high * s
          }
          # --- Scaled to % if the model is binomial or multinomial ---
          if ((is_binom || is_multinom) && nrow(cont_df) > 0) {
            pct_cols <- intersect(c("estimate","std.error","conf.low","conf.high"), names(cont_df))
            for (cc in pct_cols) cont_df[[cc]] <- cont_df[[cc]] * 100
          }

          # slopes do not depend on the contrast scheme; reuse them as base contrasts
          cont_df$estimate.base  <- cont_df$estimate
          cont_df$std.error.base <- cont_df$std.error
          cont_df$p.value.base   <- cont_df$p.value

          # Labels and cleaning (Soporte robusto para Multinomial Target)
          cont_df$Source <- if (isTRUE(stdcov)) paste0(lab, " (std)") else lab
          
          if ("group" %in% names(cont_df)) {
            cont_df$Target <- paste0(target, ":", cont_df$group)
          } else if ("class" %in% names(cont_df)) {
            cont_df$Target <- paste0(target, ":", cont_df$class)
          } else if ("response" %in% names(cont_df)) {
            cont_df$Target <- paste0(target, ":", cont_df$response)
          } else {
            cont_df$Target <- target_label
          }
          
          for (nm in c("variable","term", "group", "class", "response")) if (nm %in% names(cont_df)) cont_df[[nm]] <- NULL
          res_list[[length(res_list) + 1]] <- cont_df
        } else {
          warning(sprintf("avg_slopes failed '%s': %s", lab, as.character(cont_res)))
        }
      }
    }
  }
  
  # Apilado de Links
  if (length(res_list) == 0) {
    Links <- data.frame(Source=character(0), Target=character(0), stringsAsFactors = FALSE)
  } else {
    Links <- res_list[[1]]
    if (length(res_list) > 1) {
      for (j in 2:length(res_list)) {
        mc <- union(names(Links), names(res_list[[j]]))
        for (nm in setdiff(mc, names(Links)))             Links[[nm]] <- NA
        for (nm in setdiff(mc, names(res_list[[j]]))) res_list[[j]][[nm]] <- NA
        Links <- rbind(Links[, mc, drop = FALSE], res_list[[j]][, mc, drop = FALSE])
      }
    }
    rownames(Links) <- NULL
  }
  
  rm_cols <- intersect(c("variable","term","contrast","predicted_lo","predicted_hi","predicted", "group", "class", "response"), names(Links))
  if (length(rm_cols) > 0) Links <- Links[, setdiff(names(Links), rm_cols), drop = FALSE]
  if(nrow(Links) > 0) Links <- Links[, c("Source","Target", setdiff(names(Links), c("Source","Target"))), drop = FALSE]
  
  # ---------------------------- #
  # 2) NODES (frequencies/% + means)
  # ---------------------------- #
  sum_rows <- list()
  add_node <- function(df) { sum_rows[[length(sum_rows) + 1]] <<- df }
  
  for (lab in term_labels) {
    if (grepl(":", lab, fixed = TRUE)) {
      parts <- strsplit(lab, ":", fixed = TRUE)[[1]]
      ok <- TRUE
      for (p in parts) {
        if (!(p %in% names(mf) && is.factor(mf[[p]]))) { ok <- FALSE; break }
        if (drop_empty && length(unique(mf[[p]])) < 2) { ok <- FALSE; break }
      }
      if (!ok) next
      df_tmp <- mf[, parts, drop = FALSE]; df_tmp$w <- w
      rhs <- paste(parts, collapse = " + ")
      form_agg <- stats::as.formula(paste("w ~", rhs))
      agg <- stats::aggregate(form_agg, data = df_tmp, FUN = function(z) sum(z, na.rm = TRUE))
      names(agg)[names(agg) == "w"] <- "Frequency"
      cat_vec <- apply(agg[, parts, drop = FALSE], 1, function(z) paste(as.character(z), collapse = "*"))
      node_df <- data.frame(
        name = paste0(paste(parts, collapse = "*"), ":", cat_vec),
        variable = paste(parts, collapse = "*"),
        Frequency = agg$Frequency,
        check.names = FALSE, stringsAsFactors = FALSE
      )
      node_df[["%"]] <- (node_df$Frequency / w_tot) * 100
      node_df$Type <- "factor_interaction"
      node_df$Mean <- NA_real_; node_df$Min <- NA_real_; node_df$Max <- NA_real_
      for (i in seq_along(parts)) node_df[[parts[i]]] <- as.character(agg[[parts[i]]])
      add_node(node_df)
    } else {
      if (!(lab %in% names(mf))) next
      if (is.factor(mf[[lab]])) {
        df_tmp <- data.frame(grp = mf[[lab]], w = w)
        agg <- stats::aggregate(w ~ grp, data = df_tmp, FUN = function(z) sum(z, na.rm = TRUE))
        names(agg) <- c("grp","Frequency")
        node_df <- data.frame(
          name = paste0(lab, ":", as.character(agg$grp)),
          variable = lab,
          Frequency = agg$Frequency,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- (node_df$Frequency / w_tot) * 100
        node_df$Type <- "factor"
        node_df$Mean <- NA_real_; node_df$Min <- NA_real_; node_df$Max <- NA_real_
        node_df[[lab]] <- as.character(agg$grp)
        add_node(node_df)
      } else if (is.numeric(mf[[lab]])) {
        x <- mf[[lab]]
        m <- wmean(x, w)
        xr <- range_fun(x, outliers = outliers, ql = q_low, qh = q_high)
        xmin <- xr[1]; xmax <- xr[2]
        norm <- if (!is.finite(xmax - xmin) || (xmax <= xmin)) NA_real_ else (m - xmin) / (xmax - xmin) * 100
        node_df <- data.frame(
          name = lab,
          variable = lab,
          Frequency = w_tot,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- norm
        node_df$Type <- "numeric"
        node_df$Mean <- m; node_df$Min <- xmin; node_df$Max <- xmax
        add_node(node_df)
      }
    }
  }
  
  # Response (binomial or multinomial positive category)
  if (!is.na(target) && target %in% names(mf)) {
    if (is.factor(mf[[target]])) {
      df_tmp <- data.frame(grp = mf[[target]], w = w)
      agg <- stats::aggregate(w ~ grp, data = df_tmp, FUN = function(z) sum(z, na.rm = TRUE))
      names(agg) <- c("grp","Frequency")
      # Filtramos por pred_cat SOLO si es binomial explícito. Si es multinomial, crea Nodos para todos.
      if (is_binom && !is.na(pred_cat)) {
        keep <- which(as.character(agg$grp) == pred_cat)
        agg  <- if (length(keep) > 0) agg[keep, , drop = FALSE] else agg[0, , drop = FALSE]
      }
      if (nrow(agg) > 0) {
        node_df <- data.frame(
          name = paste0(target, ":", as.character(agg$grp)),
          variable = target,
          Frequency = agg$Frequency,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- (node_df$Frequency / w_tot) * 100
        node_df$Type <- "response_factor"
        node_df$Mean <- NA_real_; node_df$Min <- NA_real_; node_df$Max <- NA_real_
        node_df[[target]] <- as.character(agg$grp)
        sum_rows[[length(sum_rows) + 1]] <- node_df
      }
    } else if (is.numeric(mf[[target]])) {
      x <- mf[[target]]
      u <- sort(unique(na.omit(x)))
      is01 <- length(u) == 2L && all(u %in% c(0, 1))
      if (is_binom && is01) {
        freq1 <- sum(w[!is.na(x) & x == 1], na.rm = TRUE)
        node_df <- data.frame(
          name = paste0(target, ":", "1"),
          variable = target,
          Frequency = freq1,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- if (is.finite(w_tot) && w_tot > 0) (freq1 / w_tot) * 100 else NA_real_
        node_df$Type <- "response_factor"
        node_df$Mean <- NA_real_; node_df$Min <- NA_real_; node_df$Max <- NA_real_
        node_df[[target]] <- "1"
        sum_rows[[length(sum_rows) + 1]] <- node_df
      } else {
        m <- wmean(x, w)
        xr <- range_fun(x, outliers = outliers, ql = q_low, qh = q_high)
        xmin <- xr[1]; xmax <- xr[2]
        norm <- if (!is.finite(xmax - xmin) || (xmax <= xmin)) NA_real_ else (m - xmin) / (xmax - xmin) * 100
        node_df <- data.frame(
          name = target,
          variable = target,
          Frequency = w_tot,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- norm
        node_df$Type <- "response_numeric"
        node_df$Mean <- m; node_df$Min <- xmin; node_df$Max <- xmax
        sum_rows[[length(sum_rows) + 1]] <- node_df
      }
    } else if (is.logical(mf[[target]]) && is_binom) {
      df_tmp <- data.frame(grp = mf[[target]], w = w)
      agg <- stats::aggregate(w ~ grp, data = df_tmp, FUN = function(z) sum(z, na.rm = TRUE))
      names(agg) <- c("grp","Frequency")
      keep <- which(as.character(agg$grp) == "TRUE")
      agg  <- if (length(keep) > 0) agg[keep, , drop = FALSE] else agg[0, , drop = FALSE]
      if (nrow(agg) > 0) {
        node_df <- data.frame(
          name = paste0(target, ":", as.character(agg$grp)),
          variable = target,
          Frequency = agg$Frequency,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        node_df[["%"]] <- (node_df$Frequency / w_tot) * 100
        node_df$Type <- "response_factor"
        node_df$Mean <- NA_real_; node_df$Min <- NA_real_; node_df$Max <- NA_real_
        node_df[[target]] <- as.character(agg$grp)
        sum_rows[[length(sum_rows) + 1]] <- node_df
      }
    }
  }
  
  # Stack NODES (column junction)
  if (length(sum_rows) == 0) {
    Nodes <- data.frame()
  } else {
    Nodes <- sum_rows[[1]]
    if (length(sum_rows) > 1) {
      for (k in 2:length(sum_rows)) {
        mc <- union(names(Nodes), names(sum_rows[[k]]))
        for (nm in setdiff(mc, names(Nodes)))         Nodes[[nm]]         <- NA
        for (nm in setdiff(mc, names(sum_rows[[k]]))) sum_rows[[k]][[nm]] <- NA
        Nodes <- rbind(Nodes[, mc, drop = FALSE], sum_rows[[k]][, mc, drop = FALSE])
      }
    }
    rownames(Nodes) <- NULL
  }
  
  # --- Nodes sanitization ---
  if (nrow(Nodes) > 0) {
    keep_cols <- c("name","variable","Frequency","%","Type","Mean","Min","Max")
    existing  <- intersect(keep_cols, names(Nodes))
    Nodes     <- Nodes[, existing, drop = FALSE]
  }
  if (nrow(Nodes) > 0) {
    mm_cols <- intersect(c("Mean","Min","Max"), names(Nodes))
    if (length(mm_cols) > 0) {
      all_na_mm <- all(vapply(mm_cols, function(cc) all(is.na(Nodes[[cc]])), logical(1)))
      if (isTRUE(all_na_mm)) {
        Nodes <- Nodes[, setdiff(names(Nodes), mm_cols), drop = FALSE]
      }
    }
  }
  
  list(Nodes = Nodes, Links = Links)
}

# --------------------------- Multi-model public function --------------------------- #
contr.gw <- function(model,
                     weights = NULL, vcov = "HC1",
                     drop_empty = TRUE, stdcov = FALSE,
                     outliers = TRUE, q_low = 0.05, q_high = 0.95,
                     data = NULL, family = gaussian(),
                     .names = NULL,
                     .key = "name",
                     .nodes_merge = c("keep_first","update","sum_numeric"),
                     .sum_cols = NULL,
                     .links_dedup = FALSE,
                     .keep_models = FALSE,
                     .auto_quasi = TRUE) {
  
  .nodes_merge <- match.arg(.nodes_merge)
  
  # Solve vector of weights with respect to `data`
  w_vector <- NULL
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L) {
      if (!is.null(data) && weights %in% names(data)) w_vector <- data[[weights]]
    } else {
      w_expr <- substitute(weights)
      w_name <- if (is.symbol(w_expr) || is.name(w_expr)) deparse(w_expr) else NULL
      if (!is.null(w_name) && !is.null(data) && w_name %in% names(data)) {
        w_vector <- data[[w_name]]
      } else if (is.numeric(weights)) {
        w_vector <- weights
      }
    }
    if (!is.null(w_vector) && !is.null(data) && length(w_vector) != nrow(data)) {
      stop("The length of 'weights' does not match 'nrow(data)'.")
    }
  }
  
  # Build model list
  if (is.character(model)) {
    modelos <- .fit_models_from_specs(spec = model, data = data,
                                      default_family = family, weights = w_vector,
                                      auto_quasi = .auto_quasi)
  } else if (.is_list_of_models(model)) {
    modelos <- model
    if (is.null(names(modelos)) || any(!nzchar(names(modelos)))) {
      if (!is.null(.names)) {
        if (length(.names) != length(modelos)) {
          stop("`.names` should be the same length as the model list.")
        }
        names(modelos) <- .names
      } else {
        names(modelos) <- paste0("Model_", seq_along(modelos))
      }
    }
  } else if (.is_model(model)) {
    modelos <- list(Model_1 = model)
    if (!is.null(.names) && length(.names) == 1) names(modelos) <- .names
  } else {
    stop("`model` must be a model object, a list of models, or a vector of string-type specifications.")
  }
  
  # Process each model and collect Nodes/Links
  all_nodes <- list()
  all_links <- list()
  
  for (i in seq_along(modelos)) {
    m_name <- names(modelos)[i]
    res <- contrastes_gw_one(modelos[[i]], weights = w_vector, vcov = vcov,
                             drop_empty = drop_empty, stdcov = stdcov,
                             outliers = outliers, q_low = q_low, q_high = q_high)
    if (nrow(res$Links) > 0) {
      res$Links$Model <- m_name
      all_links[[length(all_links) + 1]] <- res$Links
    }
    if (nrow(res$Nodes) > 0) {
      all_nodes[[length(all_nodes) + 1]] <- res$Nodes
    }
  }
  
  # Assemble unified matrix
  Links <- .rbind_union(all_links)
  Nodes <- NULL
  if (length(all_nodes) > 0) {
    Nodes <- all_nodes[[1]]
    if (length(all_nodes) > 1) {
      for (j in 2:length(all_nodes)) {
        Nodes <- .merge_nodes(Nodes, all_nodes[[j]], .key = .key, .mode = .nodes_merge, .sum_cols = .sum_cols)
      }
    }
  }
  
  if (isTRUE(.links_dedup) && nrow(Links) > 0) {
    Links <- unique(Links)
  }
  
  out <- list(Nodes = Nodes, Links = Links)
  if (isTRUE(.keep_models)) out$models <- modelos
  
  return(out)
}