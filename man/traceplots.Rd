\name{traceplots}
\alias{traceplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create MCMC trace plots
}
\description{
Create MCMC trace plots for the transmission parameters
}
\usage{
traceplots(mcmcoutput, tracecol="lightgrey", labels=c("p","z","beta","gamma",
"gamma_G","genpar"), burn=1, filter=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcoutput}{
MCMC output from the \code{\link{transmission_analysis}} function.
}
  \item{tracecol}{
Trace color.
}
  \item{labels}{
Character vector for parameter names.
}
  \item{burn}{
Burn-in length. This is not plotted.
}
  \item{filter}{
Retain every \code{n}th iteration.
}
}
\details{
Red lines indicate the mean (solid) and 95\% credible intervals (dashed) for each parameter.
}
\examples{
  data(mcmcoutput)
  traceplots(mcmcoutput)
}
