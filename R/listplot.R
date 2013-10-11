# plot objects given in a list

#@title Reassemble a list of arguments according to a list of names
#@description Reorder and reassemble lists of arguments, according to a given vector
# of name tags.
#@param arglist list, the list that is to be reassembled,
#@param tagnames character, the names to be sortet after.
#@return a list of lists with names c(sortnames, "..."). 
#@details The input \code{arglist} is a list of lists or simple elements. Any lists contained in 
#arglist are searched for names contained in \code{tagnames}. If elements with appropriate
# names are found, they constitute the new lists returned by the function. 
# Elements in arglist that are not lists are pooled in the "..." component of the
# returned list.
#
# Note: to avoid taking list arguments into parts, enclose them in a list, see the Examples
# This function was written for processing plot arguments.
# @author Ute Hahn,  \email{ute@@imf.au.dk}
#@export
#@examples
#arglist <- list(x = "a", y = list(a = 1, b = 2), 
#                z = list(a = 3, b = 4, c = 5))
#nametags <- c("a", "b", "d", "...")
#newlist <- reorderList(arglist, nametags)
#str(newlist$a)
##unassigned elements
#str(newlist$...)
## protecting list argument y
#arglist.p <- list(x= "a", y = list(list(a = 1, b = 2)), 
#                  z = list(a = 3, b = 4, c = 5))
#newlist.p <- reorderList(arglist.p, nametags)
#str(newlist.p$...)
#' @rdname sostatpp-internal
#' @keywords internal


retagList <- function (arglist, tagnames)
{
  argnames <- names(arglist)
  result <- as.list(rep(list(list()), length(tagnames)))
  names(result) <- tagnames
  # provide space for the unnamed ones
  result$... <- list() 
  if (length(arglist) > 0) for(i in 1 : length(arglist))
  {
    arg <- arglist[[i]]
    
    if (is.list(arg))
    {
      namesarg <- names(arg)
      if ((length(arg) == 1) && is.null(namesarg)) result[["..."]][[argnames[i]]] <- arg[[1]]
      else{  
      # sort its elements and distribute them to the rights slots in the result
      # use only known names
      slots <- match(names(arg), tagnames)
      for (j in (1 : length(arg)))
      {
        if (!is.na(slots[j])) result[[slots[j]]][[argnames[i]]]  <- arg[[j]]
      }
    }}
    else 
      result[["..."]][[argnames[i]]] <- arg
  }
  return(result)
}

#'@title Plot objects in a list
#'@description Plot all objects given in a list, with parameters that also can be given as a list.
#'@param objects named list of \code{R}-objects.
#'@param allinone logical, if TRUE, all objects are plotted in one window.
#'@param ... parameters and parameterlists passed to plot method of the objects. 
#'@details 
#'Only objects that belong to classes with a plot method are plotted.
#'The plot parameters may be given as lists with the same name as the objects list.
#' Parameters in named lists are assigned to the object with same name in
#'list \code{objects}. Plot parameters not given in a list apply to all objects.
#'
#'Any \code{add}-parameters only affect the first plotted object. Whether or not
#'the plots of the remaining objects are added to the first plot, is controlled 
#'by parameter \code{allinone.}
#'
#'@export
#'@examples
#'# a list of plottable objects: functions
#'curves <- list(b = cos, a = sin)
#'lplot(curves, col = list(a = "green", b = "red"), ylim = c(-1, 1), to = pi)
#'# start a new plot for every object
#'lplot(curves, allinone = FALSE, col = list(a = "green", b = "red"), 
#'      ylim = c(-1, 1), to = pi)
#' @author Ute Hahn,  \email{ute@@imf.au.dk}


lplot <- function (objects=NULL, allinone = TRUE, ...)
{
  stopifnot(is.list(objects))
  
  # check whether objects have plot methods
  plotmethods <- methods(plot)
  plottingclasses <- sapply(plotmethods, function(s) substr(s, 6, 200))
  canplot <- function(obj) any(class(obj) %in% plottingclasses)
  objects <- objects[sapply(objects, canplot)]
 
  nobjects <- length(objects)
  if (nobjects == 0) stop("nothing to plot here")
  unnamed <- is.null(obnames <- names(objects))
  
  dotargs <- list(...)
  
  addfirst <- any(dotargs$add)
  #override add arguments: plot either all in one, or all separately
  dotargs$add <- NULL # plot all in one plot
   
  orderedArgs <- retagList(dotargs, obnames)
  do.call(plot, c(list(objects[[1]]), orderedArgs[[1]], 
                  orderedArgs[["..."]], add = addfirst))
  
  if (nobjects > 1) for (i in 2 : nobjects)
     do.call(plot, c(list(objects[[i]]), orderedArgs[[i]], 
       orderedArgs[["..."]], add = allinone))
}
