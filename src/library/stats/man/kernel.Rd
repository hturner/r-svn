% File src/library/stats/man/kernel.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2024 R Core Team
% Distributed under GPL 2 or later

\name{kernel}
\alias{kernel}
\alias{bandwidth.kernel}
\alias{df.kernel}
\alias{is.tskernel}
\alias{plot.tskernel}
\title{Smoothing Kernel Objects}
\description{
  The \code{"tskernel"} class is designed to represent discrete
  symmetric normalized smoothing kernels.  These kernels can be used to
  smooth vectors, matrices, or time series objects.

  There are \code{\link{print}}, \code{\link{plot}} and \code{\link{[}}
  methods for these kernel objects.
}
\usage{
kernel(coef, m = 2, r, name)

df.kernel(k)
bandwidth.kernel(k)
is.tskernel(k)

\method{plot}{tskernel}(x, type = "h", xlab = "k", ylab = "W[k]",
     main = attr(x,"name"), \dots)
}
\arguments{
  \item{coef}{the upper half of the smoothing kernel coefficients
    (including coefficient zero) \emph{or} the name of a kernel
    (currently \code{"daniell"}, \code{"dirichlet"}, \code{"fejer"} or
    \code{"modified.daniell"}).}
  \item{m}{the kernel dimension(s) if \code{coef} is a name.  When \code{m}
    has length larger than one, it means the convolution of
    kernels of dimension \code{m[j]}, for \code{j in 1:length(m)}.
    Currently this is supported only for the named \code{"*daniell"} kernels.}
  \item{name}{the name the kernel will be called.}
  \item{r}{the kernel order for a \I{Dirichlet} or \I{Fejer} kernel.}
  \item{k, x}{a \code{"tskernel"} object.}
  \item{type, xlab, ylab, main, \dots}{arguments passed to
    \code{\link{plot.default}}.}
}
\details{
  \code{kernel} is used to construct a general kernel or named specific
  kernels.  The modified \I{Daniell} kernel halves the end coefficients.

  The \code{\link{[}} method allows natural indexing of kernel objects
  with indices in \code{(-m) : m}.  The normalization is such that for
  \code{k <- kernel(*)}, \code{sum(k[ -k$m : k$m ])} is one.

  \code{df.kernel} returns the \sQuote{equivalent degrees of freedom} of
  a smoothing kernel as defined in
  \bibcite{Brockwell and Davis (1991), page 362},
  and \code{bandwidth.kernel} returns the equivalent bandwidth as
  defined in \bibcite{Bloomfield (1976), p.\sspace{}201},
  with a continuity correction.
}
\value{
  \code{kernel()} returns an object of class \code{"tskernel"} which is
  basically a list with the two components \code{coef} and the kernel
  dimension \code{m}.  An additional attribute is \code{"name"}.
}
\author{A. Trapletti; modifications by B.D. Ripley}
\seealso{
    \code{\link{kernapply}}
}
\references{
  Bloomfield, P. (1976)
  \emph{Fourier Analysis of Time Series: An Introduction.}
  Wiley.

  Brockwell, P.J. and Davis, R.A. (1991)
  \emph{Time Series: Theory and Methods.}
  Second edition. Springer, pp.\sspace{}350--365.
}
\examples{
require(graphics)

## Demonstrate a simple trading strategy for the
## financial time series German stock index DAX.
x <- EuStockMarkets[,1]
k1 <- kernel("daniell", 50)  # a long moving average
k2 <- kernel("daniell", 10)  # and a short one
plot(k1)
plot(k2)
x1 <- kernapply(x, k1)
x2 <- kernapply(x, k2)
plot(x)
lines(x1, col = "red")    # go long if the short crosses the long upwards
lines(x2, col = "green")  # and go short otherwise

## More interesting kernels
kd <- kernel("daniell", c(3, 3))
kd # note the unusual indexing
kd[-2:2]
plot(kernel("fejer", 100, r = 6))
plot(kernel("modified.daniell", c(7,5,3)))

# Reproduce example 10.4.3 from Brockwell and Davis (1991)
spectrum(sunspot.year, kernel = kernel("daniell", c(11,7,3)), log = "no")
}
\keyword{ts}
