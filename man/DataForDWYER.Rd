% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataForDWYER.R
\docType{data}
\name{DWYERcat}
\alias{DWYERcat}
\title{Data for function dwyerCA}
\format{
\code{DWYERcat} is a data.frame with columns \code{t} (strength), \code{v} (number of levels),
\code{k} (number of CA columns), code{N} (number of CA rows) and \code{Source}
(a brief descriptive source entry by Aaron Dwyer).
}
\source{
The source is Dwyer (2024), who cites personal communications with Colbourn as the
source of many of the arrays; there is some overlap between this repository and
\code{\link{CKRScat}}.\cr
The creation of the object is documented by an R file in the \code{extdata}
folder of this package, which can be located using
\code{system.file("extdata", package="CAs")}; the many
raw arrays are not provided there, as they are available on the web.
}
\usage{
DWYERcat
}
\description{
The data.frame DWYERcat supports the identification and download of CAs from the
DWYER GitHub repository via function dwyerCA.
}
\section{Details}{

The repository is large and has not been included into the package.
Function \code{\link{dwyerCA}} calculates the requested array
(in case of an OA) or downloads it on the fly.
}

\examples{
head(DWYERcat) ## will be constructed with function SCA_Busht
tail(DWYERcat) ## will be downloaded on the fly, except the last one

}
\references{
Dwyer (2024)
}
\keyword{datasets}
