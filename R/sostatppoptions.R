.onLoad <- function(libname, pkgname) {
    assign(".sostatppStuff", new.env(), envir=parent.env(environment()))
    Kstyle = list(
      theo = simplist(col = "black", lty = "dashed"),
      iso = simplist(col = "red", lty = "solid"),
      trans = simplist(col = "blue", lty = "solid"),
      border = simplist(col = "green", lty = "solid")
      )
    assign("opt", simplist(Kstyles = Kstyle), .sostatppStuff)
}

#'@title Options for package sostatpp
#'@description Get or set options for package sostatpp
#'@param ... option names given as character string, or given as \code{name = value}.
#'Alternatively, options can be passed as \code{\link{simplist}}s.
#'@return For \code{ploptions()}, a \code{simplist} of set options for package \code{plutils};
#'if invoked with names of options, a list of all options contained in the arguments. 
#'If called with \code{name = value}, no visible result is returned.
#'@details 
#'It is possible to set any options - this might be useful for global options of
#'other packages. This can however also be done with the base function \code{options}, the 
#'only difference is that \code{sostatppoptions} returns a \code{\link{simplist}}, and resolves
#'\code{\link{simplist}} arguments, but not \code{\link{list}} arguments.
#'
#'Currently, only one options is used by package \code{sostatpp}:
#'\tabular{ll}{
#'\code{Kstyles}\tab list with elements \code{theo}, \code{iso}, \code{trans} and
#'\code{border}, that control plotting of $K$-functions obtained by \code{\link{estK}},
#'or other \code{funsamples} containing elements \code{theo}, \ldots, \code{border}.
#'By default, each element of \code{Kstyles} is a \code{simplist} with entries \code{col} and
#'\code{lty}.
#'}
#'@seealso The base function \code{\link{options}} is similar.
#'@export
#@examples
#sostatppoptions() 
## add a new option
#ploptions(schnurz = "piep")
#ploptions()
## save options
#oldopt <- ploptions()
## change old option, and add a new one
#ploptions(schnurz = pi, a = 5)
#ploptions()
## reset old options, and change the new one from the step before
#ploptions(oldopt, a = 3)
#ploptions()
sostatpp.options <- function(...)
{
  arglist <- simplist(...)
  opt <- get("opt", envir = .sostatppStuff)
  if (!length(arglist)) return (opt)
  nams <- names(arglist) 
  if (is.null(nams))
    nams <- rep("", length(arglist))
  reportonly <- sapply(nams, function(nm) identical(nm, ""))
  invisbl <- (sum(!reportonly) > 0)
   
  updates <- arglist[!reportonly]
  firstclass(updates) <- "simplist" 
  # should perhaps change behaviour of simplist when subsetting
  
  newopt <- simplist(opt, updates)
    
 # now check which elements to output
  for (i in seq_along(nams)) {
    if (identical(nams[i], "") && is.character(arglist[[i]]))
      nams[i] <- arglist[[i]]
  }
  result <- opt[nams[nams %in% names(opt)]]
  if (identical(length(arglist), 1L) && length(result)) 
    result <- result[[1]]
  else 
    firstclass(result) <- "simplist"
  assign("opt", newopt, envir = .sostatppStuff)
  if (invisbl) invisible(result) else return(result)
}


#'@title Set plot style for second order functions
#'@description Set default color, line type, etc for plotting second-order functions like
#'the \eqn{K}-function, according to border correction, or fot the theoretical
#'@param ... named arguments for plot styles corresponding to edge corrections.
#'Possible names are \code{theo}, \code{iso}, \code{trans} and code{border}
#'@export
#'@examples
#'setplotstyle(theo = simplist(col = "purple", lty = "dotted"))
#'K <- estK(rpoispp(100), correction = "trans")
#'plot(K)
setplotstyle <- function(...)
{
  Kstyle <- sostatpp.options("Kstyles")
 # class(Kstyle) <- "list"
  newstyle <- list(...)
  for (nam in names(newstyle)) {
    if (nam %in% names(Kstyle))
      Kstyle[[nam]] <- updateList(Kstyle[[nam]], newstyle[[nam]])
  }
  sostatpp.options(Kstyles = Kstyle)
  invisible(NULL)
}