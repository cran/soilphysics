\name{maxcurv}
\alias{maxcurv}

\title{ Maximum Curvature Point }

\description{
Function to determine the maximum curvature point of an univariate nonlinear
function of x.
}

\usage{ maxcurv(x.range, fun, graph = TRUE, ...) }

\arguments{
  \item{x.range}{ a numeric vector of length two, the range of x. }
  \item{fun}{ a function of x; it must be a one-line-written function, with no curly braces '\{\}'. }
  \item{graph}{ logical; if TRUE (default) a curve of \code{fun} is plotted. }
  \item{\dots}{ further graphical arguments. }
}


\value{
A list of
  \item{fun}{ the function of x. }
  \item{x0}{ the x critical value. }
  \item{y0}{ the y critical value. }
}

\author{ Anderson Rodrigo da Silva <anderson.agro@hotmail.com> }

\section{Warning }{ 
Numerical values into \code{fun} with three or more digits may cause miscalculation.
}

\seealso{
\code{\link{function}}, \code{\link{curve}}
}

\examples{
# an exponential model
f <- function(x) exp(-x)
maxcurv(x.range = c(-2, 5), fun = f)

# Gompertz Growth Model
Asym <- 8.5
b2 <- 2.3
b3 <- 0.6
g <- function(x) Asym * exp(-b2 * b3 ^ x)
maxcurv(x.range = c(-5, 20), fun = g)

# Lessman & Atkins (1963) model for optimum plot size
a = 40.1
b = 0.72
cv <- function(x) a * x^-b
maxcurv(x.range = c(1, 50), fun = cv)

# End (not run)
}
