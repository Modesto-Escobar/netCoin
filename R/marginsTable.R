## ===========================================================================
## marginsTable(G): marginal-effects table built ONLY from G, the object
## returned by glmCoin() of the netCoin package. Uses base R only.
##
##   Reads G$nodes (frequency and % of each category) and
##         G$links (estimate = average marginal effect in percentage points,
##                  std.error and p.value).
##
##   Marg.(category) = overall % of the dependent + marginal effect.
##
##   by = "category" : one row per category (n = frequency of the category).
##   by = "variable" : one row per variable, using its representative category
##                     (smallest p; or largest effect if criterion = "effect").
##
##   contrast = "gw"   : effects relative to the weighted overall mean (GW).
##   contrast = "base" : effects relative to the first (base) category of each
##                       qualitative variable (the base row shows "ref.").
##                       Needs the estimate.base columns stored in G$links.
##
##   columns : which columns to display; any subset of "n" (global column)
##             and "marg", "effect", "se", "sig" (per-dependent block).
##
##   sig  : significance thresholds for the stars (default .05, .01, .001;
##          sorted internally, so sig = c(.1, .05, .01, .001) adds a 4th star).
##   note : if TRUE, a legend of the stars is printed below the table.
##
##   Several dependent variables are supported: each Target produces its own
##   block of columns (Marg./M.Ef./St.Er./sig.) headed by the dependent's name.
##
##   Export (base R only):
##     file / format = "text" (.txt) | "latex" (.tex) | "csv" (.csv, opens
##     directly in Excel). If 'file' has an extension, 'format' is inferred.
##     True .xlsx, PDF and Word need extra packages and are not provided here.
##
##   Footer rows: statistics of each model taken from G$equations, one value
##   per model even with multinomial responses. 'footer' selects the columns
##   (default n and pseudo-R2; e.g. footer = c("n", "r2.nagelkerke", "aic")).
##
##   Usage:
##     source("marginsTable.R")
##     G <- glmCoin(eq, data)
##     marginsTable(G)                                  # full, by category
##     marginsTable(G, by = "variable")                 # summary, by variable
##     marginsTable(G, file = "table.tex")              # export to LaTeX
##     marginsTable(G, file = "table.csv")              # export to Excel (CSV)
## ===========================================================================
marginsTable <- function(G,
                         by        = c("category", "variable"),
                         criterion = c("p", "effect"),
                         contrast  = c("gw", "base"),
                         columns   = c("n", "marg", "effect", "se", "sig"),
                         sig       = c(0.05, 0.01, 0.001),
                         note      = FALSE,
                         digits    = 1,
                         footer    = c("n", "r2"),
                         file      = NULL,
                         format    = NULL) {
  by        <- match.arg(by)
  criterion <- match.arg(criterion)
  contrast  <- match.arg(contrast)
  columns   <- match.arg(columns, several.ok = TRUE)
  showN     <- "n" %in% columns
  blockSel  <- c(marg = "Marg.", effect = "M.Ef.", se = "St.Er.", sig = "sig.")
  blockSel  <- unname(blockSel[intersect(names(blockSel), columns)])
  if (!length(blockSel))
    stop("'columns' must include at least one of 'marg', 'effect', 'se', 'sig'.")
  k <- length(blockSel)
  nodes <- G$nodes; links <- G$links
  if (contrast == "base" &&
      !all(c("estimate.base", "std.error.base", "p.value.base") %in% names(links)))
    stop("G$links does not include base-category contrasts; ",
         "rebuild G with the current version of glmCoin().")

  varOf  <- setNames(nodes$variable,  nodes$name)   # node name -> variable
  freqOf <- setNames(nodes$Frequency, nodes$name)   # node name -> frequency
  pctOf  <- setNames(nodes[["%"]],    nodes$name)   # node name -> percentage
  typeOf <- setNames(nodes$Type,      nodes$name)   # node name -> type

  if (!is.numeric(sig) || !length(sig) || anyNA(sig) || any(sig <= 0 | sig >= 1))
    stop("'sig' must be numeric thresholds strictly between 0 and 1.")
  sig <- sort(unique(sig), decreasing = TRUE)
  stars <- function(p) vapply(p, function(pp)
    if (is.na(pp)) "" else strrep("*", sum(pp < sig)), "")
  noteTxt <- if (isTRUE(note))
    paste(vapply(seq_along(sig), function(i)
      paste0(strrep("*", i), " p<", format(sig[i], scientific = FALSE)), ""),
      collapse = ", ")

  targets <- unique(links$Target)                   # dependent variables
  deps    <- sub(":.*$", "", targets)
  # multinomial: several targets share a dependent; label blocks by category
  blockLab <- ifelse(duplicated(deps) | duplicated(deps, fromLast = TRUE),
                     targets, deps)
  totalN  <- function(t) if (grepl("factor", typeOf[t]))
                           round(freqOf[t] / (pctOf[t] / 100)) else round(freqOf[t])
  varOrder <- unique(varOf[links$Source[links$Target == targets[1]]])

  # long table: one row per (category, dependent)
  rows <- do.call(rbind, lapply(targets, function(t) {
    li <- links[links$Target == t, ]
    data.frame(variable = varOf[li$Source],
               category = sub("^[^:]*:", "", li$Source),
               n        = freqOf[li$Source],
               dep      = sub(":.*$", "", t),
               target   = t,
               marg     = pctOf[t] + li$estimate,    # marginal = total% + GW effect
               effect   = if (contrast == "base") li$estimate.base  else li$estimate,
               se       = if (contrast == "base") li$std.error.base else li$std.error,
               p        = if (contrast == "base") li$p.value.base   else li$p.value,
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  # collapse to one representative category per variable (shared across deps)
  if (by == "variable") {
    repCat <- tapply(seq_len(nrow(rows)), rows$variable, function(ix) {
      d <- rows[ix, ]
      rows$category[ix][ if (criterion == "p") which.min(d$p)
                         else which.max(abs(d$effect)) ]
    })
    rows <- rows[rows$category == repCat[rows$variable], ]
  }

  fmtCell <- function(d) {
    num <- function(fmt, v) ifelse(is.na(v), "", sprintf(fmt, v))
    isRef <- !is.na(d$effect) & d$effect == 0 & is.na(d$se)  # base category
    cbind(
      "Marg."  = num(paste0("%.", digits, "f"), d$marg),
      "M.Ef."  = ifelse(isRef, "ref.", num(paste0("%+.", digits, "f"), d$effect)),
      "St.Er." = num("(%.2f)", d$se),
      "sig."   = stars(d$p))
  }

  rows$key <- paste(rows$variable, rows$category, sep = "\r")
  base <- unique(rows[, c("variable", "category", "n", "key")])
  base <- base[order(match(base$variable, varOrder), base$category), ]

  label <- if (by == "variable") paste0(base$variable, " (", base$category, ")")
           else paste0(base$variable, ": ", base$category)

  # cell matrix (as text) + group header (dependent) + sub-headers
  M     <- cbind(Variable = label)
  group <- ""; sub <- "Variable"
  if (showN) {
    M     <- cbind(M, "(n)" = as.character(base$n))
    group <- c(group, ""); sub <- c(sub, "(n)")
  }
  for (i in seq_along(targets)) {
    di    <- rows[rows$target == targets[i], ]
    block <- fmtCell(di[match(base$key, di$key), ])[, blockSel, drop = FALSE]
    M     <- cbind(M, block)
    group <- c(group, rep(blockLab[i], k)); sub <- c(sub, colnames(block))
  }

  # Total row on top
  totalRow <- c("Total", if (showN) as.character(totalN(targets[1])))
  for (i in seq_along(targets)) {
    cell <- rep("", k)
    if ("Marg." %in% blockSel)
      cell[match("Marg.", blockSel)] <- sprintf(paste0("%.", digits, "f"), pctOf[targets[i]])
    totalRow <- c(totalRow, cell)
  }
  M <- rbind(totalRow, M)

  # footer rows from G$equations: one value per model, printed only under the
  # first block of each dependent (multinomial blocks share a single value)
  eqs <- G$equations
  if (!is.null(eqs) && !is.null(eqs$model) && length(footer)) {
    unknown <- setdiff(footer, names(eqs))
    if (length(unknown))
      warning("column(s) not found in G$equations: ",
              paste(unknown, collapse = ", "))
    lhs  <- gsub("`", "", trimws(sub("~.*$", "", eqs$model)))
    labs <- c(n = "n", n_obs = "n (obs.)", r2 = "R2",
              r2.mcfadden = "R2 McFadden", r2.nagelkerke = "R2 Nagelkerke",
              logLik = "logLik", aic = "AIC", bic = "BIC",
              deviance = "Deviance", null.deviance = "Null dev.",
              df.residual = "df resid.", family = "Family")
    fmtVal <- function(col, v) {
      if (is.na(v)) return("")
      if (!is.numeric(v)) return(as.character(v))
      fmt <- if (col %in% c("n", "n_obs", "df.residual")) "%.0f"
             else if (grepl("^r2", col)) "%.3f" else "%.2f"
      sprintf(fmt, v)
    }
    firstBlock <- !duplicated(deps)          # first target of each dependent
    for (col in intersect(footer, names(eqs))) {
      vals <- eqs[[col]][match(deps, lhs)]
      row  <- c(if (col %in% names(labs)) labs[[col]] else col, if (showN) "")
      for (i in seq_along(targets))
        row <- c(row, if (firstBlock[i]) fmtVal(col, vals[i]) else "", rep("", k - 1))
      M <- rbind(M, row)
    }
  }

  colnames(M) <- sub
  out <- structure(M, group = group, note = noteTxt,
                   class = c("marginsTable", "matrix"))

  # optional export ----------------------------------------------------------
  if (!is.null(file) || !is.null(format)) {
    if (is.null(format)) format <- tolower(sub(".*\\.", "", file))
    format <- switch(format,
                     text = , txt = "text",
                     latex = , tex = "latex",
                     csv = , excel = , xlsx = , xls = "csv",
                     stop("Unsupported format: ", format,
                          ". Use 'text', 'latex' or 'csv' (Excel)."))
    if (is.null(file)) file <- paste0("marginsTable.",
                                      c(text = "txt", latex = "tex", csv = "csv")[format])
    text <- switch(format,
                   text  = .marginsText(out),
                   latex = .marginsLatex(out),
                   csv   = .marginsCsv(out))
    con <- file(file, open = "w", encoding = "UTF-8"); writeLines(text, con); close(con)
    message("Table exported to '", file, "' (", format, ").")
    return(invisible(out))
  }
  out
}

## --- shared fixed-width renderer (used by print and text export) -----------
.marginsRender <- function(x) {
  group <- attr(x, "group"); sub <- colnames(x)
  M <- unclass(x); M[is.na(M)] <- ""; sep <- "  "; left <- 1
  w <- vapply(seq_len(ncol(M)), function(j) max(nchar(sub[j]), nchar(M[, j])), integer(1))
  # widen blocks whose group label is wider than the columns it spans
  j <- 1
  while (j <= ncol(M)) {
    if (group[j] == "") { j <- j + 1; next }
    k <- which(group == group[j]); k <- k[k >= j]
    need <- nchar(group[j]) - (sum(w[k]) + nchar(sep) * (length(k) - 1))
    if (need > 0) {
      w[k] <- w[k] + need %/% length(k)
      extra <- need %% length(k)
      if (extra) w[k[seq_len(extra)]] <- w[k[seq_len(extra)]] + 1L
    }
    j <- max(k) + 1
  }
  padC <- function(s, k) { t <- max(0L, k - nchar(s)); l <- t %/% 2
                           paste0(strrep(" ", l), s, strrep(" ", t - l)) }
  al <- function(s, j) if (j == left) formatC(s, width = w[j], flag = "-")
                       else formatC(s, width = w[j])
  banner <- character(0); j <- 1
  while (j <= ncol(M)) {
    if (group[j] == "") { banner <- c(banner, strrep(" ", w[j])); j <- j + 1 }
    else { k <- which(group == group[j]); k <- k[k >= j]
           banner <- c(banner, padC(group[j], sum(w[k]) + nchar(sep) * (length(k) - 1)))
           j <- max(k) + 1 }
  }
  c(paste(banner, collapse = sep),
    paste(vapply(seq_len(ncol(M)), function(j) al(sub[j], j), ""), collapse = sep),
    vapply(seq_len(nrow(M)),
           function(r) paste(vapply(seq_len(ncol(M)), function(j) al(M[r, j], j), ""),
                             collapse = sep), ""),
    if (!is.null(attr(x, "note"))) attr(x, "note"))
}

print.marginsTable <- function(x, ...) { cat(.marginsRender(x), sep = "\n"); invisible(x) }
.marginsText <- function(x) .marginsRender(x)

## --- LaTeX (booktabs; \multicolumn for the dependent-variable headers) -----
.marginsLatex <- function(x) {
  group <- attr(x, "group"); sub <- colnames(x)
  M <- unclass(x); M[is.na(M)] <- ""
  esc <- function(s) gsub("([&%#_$])", "\\\\\\1", s)
  nc  <- ncol(M)

  banner <- character(0); mids <- character(0); j <- 1
  while (j <= nc) {
    if (group[j] == "") { banner <- c(banner, ""); j <- j + 1 }
    else { k <- which(group == group[j]); k <- k[k >= j]
           banner <- c(banner, sprintf("\\multicolumn{%d}{c}{%s}", length(k), esc(group[j])))
           mids <- c(mids, sprintf("\\cmidrule(lr){%d-%d}", min(k), max(k)))
           j <- max(k) + 1 }
  }
  row <- function(v) paste(vapply(v, esc, ""), collapse = " & ")
  c(paste0("\\begin{tabular}{l", strrep("r", nc - 1), "}"),
    "\\toprule",
    paste0(paste(banner, collapse = " & "), " \\\\"),
    paste(mids, collapse = " "),
    paste0(row(sub), " \\\\"),
    "\\midrule",
    vapply(seq_len(nrow(M)), function(r) paste0(row(M[r, ]), " \\\\"), ""),
    "\\bottomrule",
    if (!is.null(attr(x, "note")))
      sprintf("\\multicolumn{%d}{l}{\\footnotesize %s} \\\\", nc,
              gsub("<", "$<$", esc(attr(x, "note")), fixed = TRUE)),
    "\\end{tabular}")
}

## --- CSV (two header rows; opens directly in Excel) ------------------------
.marginsCsv <- function(x) {
  group <- attr(x, "group"); sub <- colnames(x)
  M <- unclass(x); M[is.na(M)] <- ""
  top <- character(ncol(M)); j <- 1
  while (j <= ncol(M)) {
    if (group[j] == "") { top[j] <- ""; j <- j + 1 }
    else { k <- which(group == group[j]); k <- k[k >= j]
           top[k] <- ""; top[min(k)] <- group[j]; j <- max(k) + 1 }
  }
  full <- rbind(top, sub, M)
  if (!is.null(attr(x, "note")))
    full <- rbind(full, c(attr(x, "note"), rep("", ncol(M) - 1)))
  q <- function(s) paste0("\"", gsub("\"", "\"\"", s), "\"")
  apply(full, 1, function(r) paste(q(r), collapse = ","))
}
