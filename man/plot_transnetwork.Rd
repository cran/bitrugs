\name{plot_transnetwork}
\alias{plot_transnetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot transmission networks
}
\description{
Plot inferred, naive and true transmission trees in a circular figure.
}
\usage{
plot_transnetwork(mcmcoutput, epidata=NULL, type=c(1,2,3), 
plotthresh=0.05, labels="ID", text.cex=1, adj=0.25, ID=NULL, n=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcoutput}{
MCMC output returned by \code{\link{transmission_analysis}} function.
}
  \item{epidata}{
Epidemiological data in the format provided by \code{\link{simulate_data}}. Required if true or naive networks are to be plotted.  
}
  \item{type}{
Type of tranmission network to plot. 1=Inferred transmission network from MCMC output; 2=naive network, in which equal weighting is provided to all potential sources of infection; 3=True network, as provided by the infection source details in the \code{epidata} object.
}
  \item{plotthresh}{
Probability above which edges  in the inferred network are to be marked.
}
  \item{labels}{
Text labels to be marked beside each node. Options are "ID", which prints patient ID, "import", which provides the probability of importation (probability that patient was infected outside of the ward prior to admission), and "secondary" which provides the estimated number of secondary infections.
}
  \item{text.cex}{
Text size.
}
  \item{adj}{
Factor governing the space between nodes and labels.
}
  \item{ID}{
Vector of patient IDs for plotting order. By default this is a random reordering of patient IDs to position nodes at random around the circle. Must contain at least one positive host ID.
}
  \item{n}{
Number of patients. This is automatically calculated from the dimensions of \code{mcmcoutput} if not specified, however, if this is nonstandard, \code{n} should be specified here.
}
}
\details{
Nodes are randomly ordered around the circle. Red corresponds to acquisitions, while black nodes are importations. For the inferred and naive networks, the color ranges from red to black depending on the probability of importation. Edges are colored and weighted according to the probability of each potential transmission route (blue & bold high probability, green & light, low probability). To generate an identical ordering in successive plots, specify an order with \code{ID}.
}
\examples{
  data(mcmcoutput)
  data(hospitaldata)
  par(mar=c(0,0,2,0))
  w <- sample(1:100,100)
  plot_transnetwork(mcmcoutput, epidata=hospitaldata$epi, 
    type=c(3,1), plotthresh=0.05, labels="ID", text.cex=1.2, 
    adj=0.5, ID=w, n=100)
}
