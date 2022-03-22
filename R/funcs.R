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
