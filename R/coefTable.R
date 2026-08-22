## coefTable() -- tabla comparativa de coeficientes de uno o varios glmCoin.
## Companera de marginsTable(): misma arquitectura (matriz de caracteres con
## atributos 'group' y 'note'), mismos renderers, sin dependencias nuevas.
##
## Lee G$coeffs (un data.frame por ecuacion, guardado por glmCoin) y NO G$links,
## que contiene efectos marginales, no coeficientes. Como $coeffs y $equations
## comparten orden, las columnas se indexan por posicion: no hay que emparejar
## cadenas de modelo (links$Model conserva el sufijo de familia y $equations no).
##
## Unidad de columna = (objeto, ecuacion, Target). Eso cubre de un golpe:
##   - jerarquicas: 1 VD, varios conjuntos de VI  -> list(r1, ..., r5)
##   - varias VD  : glmCoin(c("y1~x", "y2~x"), d) -> un objeto, varias ecuaciones
##   - politomicas: Target = "dep:categoria"      -> varios Target por ecuacion
##   - cualquier combinacion de las anteriores
##
## Estandarizacion (columna "sta"), en dos ejes independientes:
##   numerador   = sd(x) del termino; los terminos factoriales llevan 1 (es
##                 decir, no se estandarizan) si stfactors = "nostandard".
##   denominador = sd(y) si la respuesta esta en su propia escala (gaussiana);
##                 sd(y*) de la latente si stfactors = "long" y el enlace es
##                 logit o probit; 1 en cualquier otro caso.
## Asi "nostandard" deja las dummies en su metrica original, "x" las estandariza
## por la sd de la dummy, y "long" ademas divide por la sd de la latente.

coefTable <- function(G,
                      labels    = NULL,
                      columns   = c("coef", "sig", "se"),
                      stfactors = c("nostandard", "x", "long"),
                      blocks    = c("entry", "variable", "none"),
                      intercept = FALSE,
                      sig       = c(0.05, 0.01, 0.001),
                      digits    = 3,
                      digits.se = digits,
                      footer    = c("n_obs", "r2", "dr2"),
                      note      = FALSE,
                      file      = NULL, format = NULL) {

  columns   <- match.arg(columns, c("coef", "se", "sig", "sta", "exp", "ci", "p"),
                         several.ok = TRUE)
  stfactors <- match.arg(stfactors)
  if (!is.list(blocks)) blocks <- match.arg(blocks)
  if (!is.numeric(digits) || !is.numeric(digits.se))
    stop("'digits' and 'digits.se' must be numeric.")

  ## ---- 1. Normalizar la entrada a lista de objetos --------------------------
  if (inherits(G, "netCoin")) G <- list(G)
  if (!is.list(G) || !length(G)) stop("'G' must be a glmCoin object or a list of them.")
  ok <- vapply(G, function(x) is.list(x) && all(c("coeffs", "equations") %in% names(x)), NA)
  if (!all(ok))
    stop("element(s) ", paste(which(!ok), collapse = ", "),
         " of 'G' have no $coeffs: rebuild them with the current version of glmCoin().")

  ## ---- 2. Indice de columnas: (objeto, ecuacion, target) --------------------
  idx <- do.call(rbind, lapply(seq_along(G), function(g) {
    cf   <- G[[g]]$coeffs
    eqs  <- G[[g]]$equations
    keep <- which(!vapply(cf, is.null, NA))
    if (!length(keep)) return(NULL)
    do.call(rbind, lapply(keep, function(e)
      data.frame(g = g, eq = e,
                 model  = if (!is.null(eqs) && e <= nrow(eqs)) eqs$model[e] else names(cf)[e],
                 target = unique(cf[[e]]$Target),
                 stringsAsFactors = FALSE)))
  }))
  if (is.null(idx) || !nrow(idx)) stop("no equation in 'G' has coefficients to show.")
  nb <- nrow(idx)

  cfOf  <- function(j) G[[idx$g[j]]]$coeffs[[idx$eq[j]]]
  eqval <- function(f) vapply(seq_len(nb), function(j) {
    eqs <- G[[idx$g[j]]]$equations
    if (is.null(eqs) || is.null(eqs[[f]])) NA_real_ else as.numeric(eqs[[f]][idx$eq[j]])
  }, 0)

  ## ---- 3. Etiquetas de columna ---------------------------------------------
  ## Solo entran en la etiqueta los componentes que efectivamente varian.
  dep  <- sub(":.*$", "", idx$target)          # VD; tras ':' va la categoria
  cat_ <- ifelse(grepl(":", idx$target), sub("^[^:]*:", "", idx$target), "")
  key  <- paste(idx$g, idx$eq)
  varyMod <- length(unique(key)) > 1
  varyDep <- length(unique(dep)) > 1
  ## La categoria solo desambigua si esa VD aporta mas de una (multinomial):
  ## en una binomial, Target es "low:1" y el ":1" sobra en el rotulo.
  useCat <- nzchar(cat_) &
    vapply(dep, function(d) length(unique(cat_[dep == d])), 1L) > 1
  if (is.null(labels)) {
    nm <- if (!is.null(names(G)) && all(nzchar(names(G)))) names(G)[idx$g] else
      paste("Model", match(key, unique(key)))
    labels <- vapply(seq_len(nb), function(j) {
      p <- c(if (varyMod) nm[j], if (varyDep) dep[j], if (useCat[j]) cat_[j])
      if (length(p)) paste(p, collapse = ": ") else "Model"
    }, "")
  }
  if (length(labels) != nb)
    stop("'labels' must have ", nb, " element(s), one per model/target column.")

  ## ---- 4. Filas: union de predictores en orden de aparicion -----------------
  allRows <- lapply(seq_len(nb), function(j) {
    d <- cfOf(j)
    d[d$Target == idx$target[j], , drop = FALSE]
  })
  srcOf  <- lapply(allRows, function(d) d$Source[d$term != "(Intercept)"])
  vars   <- unique(unlist(srcOf))
  hasInt <- any(vapply(allRows, function(d) any(d$term == "(Intercept)"), NA))

  ## Los Source politomicos vienen como "variable:categoria"; los dummies, planos.
  varOf <- unlist(lapply(G, function(x) setNames(x$nodes$variable, x$nodes$name)))
  varOf <- varOf[!duplicated(names(varOf))]
  vvar  <- ifelse(is.na(varOf[vars]), sub(":.*$", "", vars), varOf[vars])
  vcat  <- ifelse(grepl(":", vars), sub("^[^:]*:", "", vars), "")
  vlab  <- ifelse(nzchar(vcat), paste0(vvar, ": ", vcat), vars)

  ## ---- 5. Celdas ------------------------------------------------------------
  sig   <- sort(unique(sig), decreasing = TRUE)
  stars <- function(p) vapply(p, function(pp)
    if (is.na(pp)) "" else strrep("*", sum(pp < sig)), "")
  num <- function(fmt, v) ifelse(is.na(v), "", sprintf(fmt, v))
  fb  <- paste0("%.", digits, "f")            # coeficientes
  fs  <- paste0("%.", digits.se, "f")         # errores tipicos e intervalos

  sub_ <- c(coef = "b", se = "(SE)", sig = "sig.", sta = "beta",
            exp = "exp(b)", ci = "95% CI", p = "p")[columns]
  k <- length(sub_)

  ## Estandarizacion: ver la cabecera del fichero.
  staOf <- function(d, j) {
    cf  <- cfOf(j)
    sdy <- attr(cf, "sd.y")
    lat <- unname(attr(cf, "latent.sd")[idx$target[j]])
    den <- if (!is.null(sdy) && !is.na(sdy)) sdy
           else if (identical(stfactors, "long") && length(lat) && !is.na(lat)) lat
           else 1
    nu  <- if (identical(stfactors, "nostandard"))
             ifelse(is.na(d$factor), NA_real_, ifelse(d$factor, 1, d$sd.x)) else d$sd.x
    d$estimate * nu / den
  }
  ## exp(b) solo tiene lectura en enlaces log (RR) y logit (odds ratio): se
  ## decide por el enlace, no por la familia, para no exponenciar un probit.
  expLink <- vapply(seq_len(nb), function(j) {
    lk <- attr(cfOf(j), "link"); !is.null(lk) && lk %in% c("log", "logit") }, NA)
  if ("exp" %in% columns && !any(expLink))
    warning("no equation uses a log or logit link; column 'exp' is left empty.",
            call. = FALSE)
  expOf <- function(d, j) if (expLink[j]) exp(d$estimate) else rep(NA_real_, nrow(d))

  body <- function(vv) do.call(cbind, lapply(seq_len(nb), function(j) {
    d <- allRows[[j]][match(vv, allRows[[j]]$Source), , drop = FALSE]
    out <- cbind(
      b        = num(fb, d$estimate),
      `(SE)`   = num(paste0("(", fs, ")"), d$std.error),
      sig.     = stars(d$p.value),
      beta     = num(fb, staOf(d, j)),
      `exp(b)` = num(fb, expOf(d, j)),
      `95% CI` = ifelse(is.na(d$estimate), "",
                        sprintf(paste0("[", fs, ", ", fs, "]"),
                                d$conf.low, d$conf.high)),
      p        = num("%.3f", d$p.value))
    out[, sub_, drop = FALSE]
  }))

  ## ---- 6. Agrupacion de filas ----------------------------------------------
  ## "entry": bloque = ecuacion en la que el predictor entra por primera vez.
  ##          Reproduce los bloques teoricos de una regresion jerarquica.
  grp <- if (identical(blocks, "none")) NULL
  else if (is.list(blocks)) {
    g <- rep(NA_character_, length(vars))
    for (b in names(blocks)) g[vars %in% blocks[[b]]] <- b
    g
  } else if (identical(blocks, "variable")) vvar
  else {
    first <- vapply(vars, function(v) which(vapply(srcOf, function(s) v %in% s, NA))[1], 1L)
    paste("Block", match(first, sort(unique(first))))
  }
  ## Si la variable ya encabeza el bloque, la fila solo lleva la categoria.
  if (identical(blocks, "variable")) vlab <- ifelse(nzchar(vcat), vcat, vars)

  ## ---- 7. Cuerpo ------------------------------------------------------------
  M <- cbind(Variable = vlab, body(vars))
  group <- c("", rep(labels, each = k))
  colnames(M) <- c("Variable", rep(unname(sub_), nb))

  if (!is.null(grp)) {                       # insertar filas de bloque
    o <- order(match(grp, unique(grp)))
    M <- M[o, , drop = FALSE]; grp <- grp[o]
    for (r in rev(which(!duplicated(grp))))
      M <- rbind(M[seq_len(r - 1), , drop = FALSE],
                 c(grp[r], rep("", ncol(M) - 1)),
                 M[seq(r, nrow(M)), , drop = FALSE])
  }
  ## La constante va al final, fuera de los bloques de predictores.
  if (isTRUE(intercept) && hasInt)
    M <- rbind(M, c("(Intercept)", body("(Intercept)")))

  ## ---- 8. Pie de ajuste -----------------------------------------------------
  labs <- c(n = "n", n_obs = "n (obs.)", k = "Predictors", r2 = "R2",
            dr2 = "Delta R2", ref = "Reference model", r2.mcfadden = "R2 McFadden",
            r2.nagelkerke = "R2 Nagelkerke", logLik = "logLik", aic = "AIC",
            bic = "BIC", deviance = "Deviance", null.deviance = "Null dev.",
            df.residual = "df resid.")
  known <- c(names(labs), names(G[[1]]$equations))
  if (length(setdiff(footer, known)))
    warning("column(s) not found in $equations: ",
            paste(setdiff(footer, known), collapse = ", "), call. = FALSE)

  r2 <- eqval("r2")
  ## Modelo de referencia: el mayor estrictamente anidado, dentro del mismo Target
  ref <- vapply(seq_len(nb), function(j) {
    c_ <- which(idx$target == idx$target[j] & seq_len(nb) != j &
                vapply(srcOf, function(s) all(s %in% srcOf[[j]]) &&
                         length(s) < length(srcOf[[j]]), NA))
    if (!length(c_)) NA_integer_ else c_[which.max(vapply(srcOf[c_], length, 1L))]
  }, 1L)

  fmtVal <- function(f, v) if (is.na(v)) "" else
    sprintf(if (grepl("^r2", f)) "%.3f"
            else if (f %in% c("n", "n_obs", "k", "df.residual")) "%.0f"
            else "%.2f", v)
  ## Varios Target de una misma ecuacion (multinomial) comparten ajuste:
  ## se muestra una sola vez, bajo la primera de sus columnas.
  firstEq <- !duplicated(paste(idx$g, idx$eq))
  for (f in footer) {
    row <- if (identical(f, "ref")) ifelse(is.na(ref), "—", labels[ref])
    else if (identical(f, "dr2")) ifelse(is.na(ref), "—", sprintf("%+.3f", r2 - r2[ref]))
    else {
      v <- if (identical(f, "k")) vapply(srcOf, length, 1L) else eqval(f)
      vapply(seq_len(nb), function(j) fmtVal(f, v[j]), "")
    }
    M <- rbind(M, c(if (f %in% names(labs)) labs[[f]] else f,
                    as.vector(rbind(ifelse(firstEq, row, ""), matrix("", k - 1, nb)))))
  }

  ## ---- 9. Aviso: muestras distintas hacen incomparables los R2 --------------
  n <- eqval("n_obs")
  difN <- length(unique(n[!is.na(n)])) > 1
  if (difN)
    warning("models are fitted on different samples (n = ",
            paste(range(n, na.rm = TRUE), collapse = "-"),
            "); R2 and Delta R2 are not directly comparable across columns.",
            call. = FALSE)

  staTxt <- if ("sta" %in% columns) switch(stfactors,
    nostandard = "beta: b x sd(x); qualitative terms not standardized",
    x          = "beta: b x sd(x)",
    long       = "beta: b x sd(x) / sd(latent y*)")
  noteTxt <- if (isTRUE(note))
    paste0(paste(vapply(seq_along(sig), function(i)
      paste0(strrep("*", i), " p<", format(sig[i], scientific = FALSE)), ""),
      collapse = ", "),
      if (!is.null(staTxt)) paste0(". ", staTxt),
      if (difN) ". Warning: models fitted on different samples.")

  rownames(M) <- NULL
  out <- structure(M, group = group, note = noteTxt,
                   class = c("coefTable", "matrix"))

  ## ---- 10. Exportacion: reutiliza los renderers de marginsTable -------------
  if (!is.null(file) || !is.null(format)) {
    if (is.null(format)) format <- tolower(sub(".*[.]", "", file))
    format <- switch(format, text = , txt = "text", latex = , tex = "latex",
                     csv = , excel = , xlsx = , xls = "csv",
                     stop("Unsupported format: ", format,
                          ". Use 'text', 'latex' or 'csv' (Excel)."))
    if (is.null(file))
      file <- paste0("coefTable.", c(text = "txt", latex = "tex", csv = "csv")[format])
    text <- switch(format, text = .marginsText(out), latex = .marginsLatex(out),
                   csv = .marginsCsv(out))
    con <- file(file, open = "w", encoding = "UTF-8")
    writeLines(text, con); close(con)
    message("Table exported to '", file, "' (", format, ").")
    return(invisible(out))
  }
  out
}

print.coefTable <- function(x, ...) { cat(.marginsRender(x), sep = "\n"); invisible(x) }
