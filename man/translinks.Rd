\name{translinks}
\alias{translinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate tree distance
}
\description{
Compute the number of transmission links separating two nodes in a transmission tree.
}
\usage{
translinks(p1, p2, patID, psource)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
Patient ID 1.
}
  \item{p2}{
Patient ID 2.
}
  \item{patID}{
Vector of patient IDs.
}
  \item{psource}{
Vector of transmission sources corresponding to \code{patID}.
}
}
\details{
Returns the number of transmission links separating two patients in the transmission tree. If patients belong to unconnected transmission trees, -1 is returned. Directions of transmission links are disregarded.
}
\examples{
data(hospitaldata)
# Patients 48 and 58 are 3 links apart
translinks(48,58,hospitaldata$epi[,1],hospitaldata$epi[,5])
# Patients 89 and 58 are in separate transmission chains
translinks(89,58,hospitaldata$epi[,1],hospitaldata$epi[,5])

}
