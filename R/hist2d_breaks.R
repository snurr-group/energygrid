# derived from base R hist2d function, but with defined breaks fed to it, not number of bins
# this new Function doesn't have the same_scale option, since both x and y bounds are fed
hist2d_break <- function(x, y, x_breaks, y_breaks, na.rm = TRUE, 
                         show = TRUE, col = c("black", heat.colors(12)), FUN = base::length, 
                         xlab, ylab, ...)
{
  if (is.null(y)) {
    if (ncol(x) != 2) 
      stop("If y is ommitted, x must be a 2 column matirx")
    y <- x[, 2]
    x <- x[, 1]
  }
  nas <- is.na(x) | is.na(y)
  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
  }
  else stop("missing values not permitted if na.rm=FALSE")
  x.cuts <- x_breaks
  y.cuts <- y_breaks
  index.x <- cut(x, x.cuts, include.lowest = TRUE)
  index.y <- cut(y, y.cuts, include.lowest = TRUE)
  m <- tapply(x, list(index.x, index.y), FUN)
  if (identical(FUN, base::length)) 
    m[is.na(m)] <- 0
  if (missing(xlab)) 
    xlab <- deparse(substitute(xlab))
  if (missing(ylab)) 
    ylab <- deparse(substitute(ylab))
  if (show) 
    image(x.cuts, y.cuts, m, col = col, xlab = xlab, ylab = ylab, 
          ...)
  midpoints <- function(x) (x[-1] + x[-length(x)])/2
  retval <- list()
  retval$counts <- m
  retval$x.breaks = x.cuts
  retval$y.breaks = y.cuts
  retval$x = midpoints(x.cuts)
  retval$y = midpoints(y.cuts)
  retval$nobs = length(x)
  retval$call <- match.call()
  class(retval) <- "hist2d"
  retval
}  
