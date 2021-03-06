% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bootTimeInference}
\alias{bootTimeInference}
\title{calculate bootstrap inference}
\usage{
bootTimeInference(ret, b, M, Delta_null = 0)
}
\arguments{
\item{ret}{refers to a returns matrix (columns = 2)}

\item{b}{blocksize to use. Choose the optimal blocksize by
\code{\link{blockSizeCalibrate}}}

\item{M}{number of bootstrap resamples (see Eq. 9 in vignette)}

\item{Delta_null}{refers to null hypothesis, where Delta means the difference
in sharpe ratio}
}
\value{
a list with estimated difference in sharpe ratio and p-value of the
test with \code{Delat_null} as null hypothesis
}
\description{
calculate bootstrap inference
}
\examples{
\dontrun{
DATA <- bootTimeInference:::ret.hedge
opt <- blockSizeCalibrate(DATA, b_vec = c(1,2,4,6,8,10), K = 5000)
bootTimeInference(DATA, b = opt$b.optimal, M = 4999)
}
DATA <- bootTimeInference:::ret.hedge
# small example, please choose an appropiate K (see vignette)
opt <- blockSizeCalibrate(DATA, b_vec = c(1,2,4,6), K = 10)
bootTimeInference(DATA, b = opt$b.optimal, M = 4999)
}

