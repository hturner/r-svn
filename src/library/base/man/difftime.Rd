% File src/library/base/man/difftime.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2025 R Core Team
% Distributed under GPL 2 or later

\name{difftime}
\title{Time Intervals / Differences}
\alias{difftime}
\alias{as.difftime}
\alias{*.difftime}
\alias{/.difftime}
\alias{Math.difftime}
\alias{Ops.difftime}
\alias{Summary.difftime}
\alias{[.difftime}
\alias{[<-.difftime}
\alias{as.double.difftime}
\alias{as.list.difftime}
\alias{c.difftime}
\alias{diff.difftime}
\alias{format.difftime}
\alias{is.numeric.difftime}
\alias{length<-.difftime}
\alias{mean.difftime}
\alias{print.difftime}
\alias{rep.difftime}
\alias{units.difftime}
\alias{units<-.difftime}
\alias{time interval} % used by mean.Rd
\alias{summary.difftime}
\description{
  Time intervals creation, printing, and some arithmetic.  The
  \code{\link{print}()} method calls these \dQuote{time differences}.
}
\usage{
\special{time1 - time2}

difftime(time1, time2, tz,
         units = c("auto", "secs", "mins", "hours",
                   "days", "weeks"))

as.difftime(tim, format = "\%X", units = "auto", tz = "UTC")

\method{format}{difftime}(x, ..., with.units = TRUE)
\method{units}{difftime}(x)
\method{units}{difftime}(x) <- value
\method{as.double}{difftime}(x, units = "auto", ...)

## Group methods, notably for round(), signif(), floor(),
## ceiling(), trunc(), abs(); called directly, *not* as Math():
\method{Math}{difftime}(x, \dots)
}
\arguments{
  \item{time1, time2}{\link{date-time} or \link[=Dates]{date} objects.}
  \item{tz}{an optional \link{time zone} specification to be used for the
    conversion, mainly for \code{"POSIXlt"} objects.}
  \item{units}{character string.  Units in which the results are
    desired.  Can be abbreviated.}
  \item{value}{character string.  Like \code{units}, except that
    abbreviations are not allowed.}
  \item{tim}{character string or numeric value specifying a time interval.}
  \item{format}{character specifying the format of \code{tim}: see
    \code{\link{strptime}}.  The default is a locale-specific time format.}
  \item{x}{an object inheriting from class \code{"difftime"}.}
  \item{\dots}{arguments to be passed to or from other methods.}
  \item{with.units}{(for the  \code{format()} method:) logical indicating
    the units should be part, e.g., "3.5 hours"; if false, the units are
    suppressed.}
}
\details{
  Function \code{difftime} calculates a difference of two date/time
  objects, \code{time1} (end) and \code{time2} (beginning), and
  returns an object of class \code{"difftime"} with an
  attribute indicating the units.  The
  \code{\link[=S3groupGeneric]{Math}} group method provides
  \code{\link{round}}, \code{\link{signif}}, \code{\link{floor}},
  \code{\link{ceiling}}, \code{\link{trunc}}, \code{\link{abs}}, and
  \code{\link{sign}} methods for objects of this class, and there are
  methods for the group-generic (see
  \code{\link[=S3groupGeneric]{Ops}}) logical and arithmetic
  operations.

  If \code{units = "auto"}, a suitable set of units is chosen, the largest
  possible (excluding \code{"weeks"}) in which all the absolute
  differences are greater than one.

  Subtraction of date-time objects gives an object of this class,
  by calling \code{difftime} with \code{units = "auto"}.  Alternatively,
  \code{as.difftime()} works on character-coded or numeric time
  intervals; in the latter case, units must be specified, and
  \code{format} has no effect.

  Limited arithmetic is available on \code{"difftime"} objects: they can
  be added or subtracted, and multiplied or divided by a numeric vector.
  In addition, adding or subtracting a numeric vector by a
  \code{"difftime"} object implicitly converts the numeric vector to a
  \code{"difftime"} object with the same units as the \code{"difftime"}
  object.  There are methods for \code{\link{mean}} and
  \code{\link{sum}} (via the \code{\link[=S3groupGeneric]{Summary}}
  group generic), and \code{\link{diff}} via \code{\link{diff.default}}
  building on the \code{"difftime"} method for arithmetic, notably
  \code{-}.

  The units of a \code{"difftime"} object can be extracted by the
  \code{units} function, which also has a replacement form.  If the
  units are changed, the numerical value is scaled accordingly.  The
  replacement version keeps attributes such as names and dimensions.

  Note that \code{units = "days"} means a period of 24 hours, hence
  takes no account of Daylight Savings Time.  Differences in objects
  of class \code{"\link{Date}"} are computed as if in the UTC time zone.

  The \code{as.double} method returns the numeric value expressed in
  the specified units.  Using  \code{units = "auto"} means the units of the
  object.

  The \code{format} method simply formats the numeric value and appends
  the units as a text string.
}
\section{Warning}{
  Because \R follows POSIX (and almost all computer clocks) in ignoring
  leap seconds, so do time differences.  So in a UTC time zone
\preformatted{
    z <- as.POSIXct(c("2016-12-31 23:59:59", "2017-01-01 00:00:01"))
    z[2] - z[1]
}
  reports \samp{Time difference of 2 secs} but 3 seconds elapsed while
  the computer clock advanced by 2 seconds.

  If you want the elapsed time interval, you need to add in any
  leap seconds for yourself.
}
\note{
  Units such as \code{"months"} are not possible as they are not of
  constant length.  To create intervals of months, quarters or years
  use \code{\link{seq.Date}} or \code{\link{seq.POSIXt}}.
}
\seealso{
  \code{\link{DateTimeClasses}}.
}
\examples{\donttest{
(z <- Sys.time() - 3600)
Sys.time() - z                # just over 3600 seconds.

## time interval between release days of R 1.2.2 and 1.2.3.
ISOdate(2001, 4, 26) - ISOdate(2001, 2, 26)

as.difftime(c("0:3:20", "11:23:15"))
as.difftime(c("3:20", "23:15", "2:"), format = "\%H:\%M") # 3rd gives NA
(z <- as.difftime(c(0,30,60), units = "mins"))
as.numeric(z, units = "secs")
as.numeric(z, units = "hours")
format(z)
}}
\keyword{utilities}
\keyword{chron}
