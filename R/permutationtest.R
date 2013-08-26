#' Studentized Permutation Test for Two Samples of Functional Data
#' 
#' Perform a studentized permutation test of equal distribution (mean) for two 
#' samples of functional data.
#' @param foos1,foos2 arrays of real-valued data, of dimension \eqn{n x m_1} and 
#'   \eqn{n x m_2}, where \eqn{m_1}, \eqn{m_2} are the group size, and \eqn{n} 
#'   is the length of the function vectors.
#' @param use.tbar logical, defaults to \code{FALSE}. Whether to apply 
#'   studentization after integration, see `Details'.
#' @param nperm \code{NULL} or an integer giving the number of random 
#'   permutations. If \code{NULL}, all permutations are used. Only feasible if 
#'   group sizes \eqn{m_1}, \eqn{m_2} do not exceed 10.
#' @details Test statistics are integral means of studentized square distances 
#'   between the group means. The test statistics are closely related, but not 
#'   equal to Hotelling's two-sample T-squared statistic. It is assumed that 
#'   the functions all have the same, equidistant arguments. Depending on the 
#'   value of \code{use.tbar}, the test statistic is either
#'   \itemize{
#'   \item \eqn{T = mean [ (mu_1(x)-mu_2(x))^2 / (s_1^2(x)/m_1 + s_2^2(x)/m_2) ]}    or
#'   \item \eqn{Tbar = mean [ (mu_1(x)-mu_2(x))^2 ] / mean[ s_1^2(x)/m_1 + s_2^2(x)/m_2 ]}
#'   }
#'   where \eqn{mu_1(x), mu_2(x)} and \eqn{s_1^2(x), s_2^2(x)} are within group 
#'   means and variances at a fixed argument \eqn{x}, and the overall
#'   the mean is taken over all \eqn{n} arguments \eqn{x}.
#'   
#'   If \code{nperm == NULL}, the exact test with all permutations is used
#'   (combinations, for symmetry reasons). This may cause memory or computing 
#'   time issues.
#'   If \code{nperm} is given as an integer, the permutations are \code{sample}d randomly, 
#'   unless \code{nperm} is larger than the number of disjoint combinations. In that
#'   case, the exact version is used. 
#' @return 
#' A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic,}
#' \item{p.value}{the p-value of the test,}
#' \item{alternative}{a character string describing the alternative hypothesis,}
#' \item{method}{a character string indicating what type of test was performed,}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'    
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @source Hahn(2012), with slight modification (using the mean instead of the 
#'    integral, in order to avoid having to pass the arguments of the functions)
#' @references Hahn, U. (2012) A Studentized Permutation Test for the Comparison 
#' of Spatial Point Patterns. \emph{Journal of the American Statistical 
#' Association},  \bold{107} (498), 754--764.
#' @keywords htest
#' @keywords robust
#' @keywords nonparametric
#' @keywords ts    

studpermut.test <- function (foos1, foos2, use.tbar=FALSE, nperm = NULL)
{
  ##### preparations ----------------
  n <- dim(foos1)[1]
  m1 <- dim(foos1)[2]
  stopifnot (m1 > 1) # need at least two per group
  stopifnot(dim(foos2)[1] == n) # dimensions
  m2 <- dim(foos2)[2]
  stopifnot (m2 > 1) # need at least two per group
  m <- m1+m2
  foos <- cbind(foos1, foos2)
  # get the permutations. 
  # If m1 == m2, break the symmetry and save half time and memory!
  
  allcomb <- is.null(nperm)
  # if nperm is larger than the actual number of combinations, also use all of them
  if (!allcomb)
  {
    ncomb <- if (m1 == m2) choose(m-1, m1-1) else choose(m, m1)
    if (ncomb < (nperm + 1)) allcomb <- TRUE
  }
  if (allcomb) 
    {
      if (m1 == m2) index1 <- rbind(1, combn(m - 1, m1 - 1) + 1)
      else index1 <- combn(m, m1)
    } else {
      if (m1 == m2) index1 <- rbind(1, replicate(nperm, sample(m - 1, m1 - 1) + 1)) 
      else index1 <- replicate(nperm, sample(m, m1)) 
      index1 <- cbind((1 : m1), index1) # the first is the original
    }
   
  # do the calculations the good old fashioned way with sums and sqs, to save time
  
  SX. <- apply (foos, 1, sum)
  SXX. <- apply (foos^2, 1, sum)
  
  Tstatistic <- function (ind) # could be further optimized in symmetric case 
  {
    SX1 <- apply(foos[, ind], 1, sum)
    SXX1 <- apply(foos[, ind]^2, 1, sum)
    SX2 <- SX. - SX1
    SXX2 <- SXX. - SXX1
    mu1 <- SX1 / m1 
    mu2 <- SX2 / m2
    ss1 <- (SXX1 - (SX1^2 / m1)) / ((m1-1) / m1)
    ss2 <- (SXX2 - (SX2^2 / m2)) / ((m2-1) / m2)
    
    if (use.tbar) return (sum((mu1 -mu2)^2) / sum((ss1 + ss2))) else 
                  return (mean((mu1 -mu2)^2 / (ss1 + ss2), na.rm=T))
  }
  
  Tvals <- apply(index1, 2, Tstatistic)
  
  pval <- mean(Tvals >= Tvals[1])           
  stat <- Tvals[1]
  names(stat) <- if(use.tbar) "Tbar" else "T"
  datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
  method <- paste("Studentized two sample permutation test for fda, using T",
                  ifelse(use.tbar, "bar", ""), sep="")
  alternative <- "samples not exchangeable"
  ptt <- list(statistic = stat, 
              p.value = pval, 
              alternative = alternative, 
              method = method, 
              data.name = datname)
  class(ptt) <- "htest"
  return(ptt)
}
  