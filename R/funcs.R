setseq <- function(seq, levs, text, equal.length = TRUE) {
  name <- deparse(substitute(seq))
  if (missing(text))
      text <- paste0("Argument '", name, "'")
  ##
  if (length(levs) != length(seq) & equal.length)
    stop("Length of argument '", name,
         "' different from number of treatments.", call. = FALSE)
  ##
  if (length(unique(seq)) != length(seq))
    stop("Values for argument '", name,
         "' must all be disparate.", call. = FALSE)
  ##
  if (is.numeric(seq)) {
    if (anyNA(seq))
      stop("Missing values not allowed in argument '",
           name, "'.", call. = FALSE)
    if (any(!(seq %in% seq_len(length(levs)))))
      stop(paste("Argument '", name,
                 "' must be a permutation of the integers from 1 to ",
                 length(levs), ".", sep = ""), call. = FALSE)
    res <- levs[seq]
  }
  else if (is.character(seq)) {
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(seq), tolower(levs), nomatch = NA)
    else
      idx <- charmatch(seq, levs, nomatch = NA)
    ##
    if (equal.length && (anyNA(idx) || any(idx == 0)))
      stop(paste(text,
                 " must be a permutation of the following values:\n  ",
                 paste(paste("'", levs, "'", sep = ""),
                       collapse = " - "), sep = ""), call. = FALSE)
    res <- levs[idx]
    if (!equal.length)
      res <- res[!is.na(res)]
  }
  else
    stop("Argument '", name, "' must be either a numeric or character vector.",
         call. = FALSE)
  
  res
}

chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(NULL)
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      stop("Argument '", name, "' must be negative.",
           call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  invisible(NULL)
}

chklogical <- function(x, name = NULL) {
  ##
  ## Check whether argument is logical
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.numeric(x))
    x <- as.logical(x)
  ##
  if (length(x) !=  1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical.", call. = FALSE)
  ##
  invisible(NULL)
}

chkclass <- function(x, class, name = NULL) {
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  n.class <- length(class)
  if (n.class == 1)
    text.class <- paste0('"', class, '"')
  else if (n.class == 2)
    text.class <- paste0('"', class, '"', collapse = " or ")
  else
    text.class <- paste0(paste0('"', class[-n.class], '"', collapse = ", "),
                    ', or ', '"', class[n.class], '"')
  ##
  if (!inherits(x, class))
    stop("Argument '", name,
         "' must be an object of class \"",
         text.class, "\".", call. = FALSE)
  ##
  invisible(NULL)
}

catch <- function(argname, matchcall, data, encl) {
  ##
  ## Catch value for argument
  ##
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
}


setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "") {
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  ##
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    ##
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  ##
  if (anyNA(idx) || any(idx == 0)) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    ##
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      ##
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, addtext, ".", call. = FALSE)
      else
        return(NULL)
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else
        return(NULL)
    }
  }
  ##
  val[idx]
}


isCol <- function(data, varname)
  !is.null(data) & varname %in% names(data)


setref <- function(reference.group, levs, length = 1,
                   varname = "reference.group", error.text) {

  if (missing(error.text)) {
    text.start <- paste0("Argument '", varname, "'")
    text.within <- paste0("argument '", varname, "'")
  }
  else {
    text.start <- paste0(toupper(substring(error.text, 1, 1)),
                         substring(error.text, 2))
    text.within <- error.text
  }
  
  
  if (length && length(reference.group) != length)
    stop(text.start,
         if (length == 1)
           " must be a numeric or a character string"
         else
           paste(" must be a numeric of character vector of length", length),
         ".",
         call. = FALSE)
  ##
  if (is.numeric(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    if (!all(reference.group %in% seq_len(length(levs))))
      stop(paste(text.start, " must ",
                 if (length == 1) "be any of the " else "contain ",
                 "integers from 1 to ",
                 length(levs), ".", sep = ""),
           call. = FALSE)
    res <- levs[reference.group]
  }
  else if (is.character(reference.group)) {
    if (any(is.na(reference.group)))
      stop("Missing value not allowed in ", text.within, ".",
           call. = FALSE)
    ##
    if (length(unique(levs)) == length(unique(tolower(levs))))
      idx <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
    else {
      idx1 <- charmatch(reference.group, levs, nomatch = NA)
      idx2 <- charmatch(tolower(reference.group), tolower(levs), nomatch = NA)
      if (anyNA(idx1) & !anyNA(idx2))
        idx <- idx2
      else
        idx <- idx1
    }
    ##
    if (anyNA(idx) || any(idx == 0))
      stop("Admissible values for ", text.within, ":\n  ",
           paste(paste("'", levs, "'", sep = ""), collapse = " - "),
           "\n  (unmatched value", if (sum(is.na(idx)) > 1) "s",
           ": ",
           paste(paste("'", reference.group[is.na(idx)], "'", sep = ""),
                 collapse = " - "),
           ")",
           call. = FALSE)
    res <- levs[idx]
  }
  
  res
}

replaceNULL <- function(x, replace = NA) {
  if (is.null(x))
    return(replace)
  x
}
