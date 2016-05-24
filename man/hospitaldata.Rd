\name{hospitaldata}
\alias{hospitaldata}
\docType{data}
\title{
Example hospital data
}
\description{
An example dataset, comprising 100 hospital admissions over 50 days, including swab results and genetic distance data.
}
\usage{data(hospitaldata)}
\format{
  A list comprising
  \describe{
\item{\code{epi}}{Matrix of epidemiological data, consisting of columns patient ID, day of admission, day of discharge, time of colonization, source of infection, and infection group.}
\item{\code{resmat}}{Matrix of test results. Each row represents the corresponding patient in \code{epi}. Each column represents a day, spanning the interval over which data were collected. 0=Negative, 1=Positive, -1=Missing}
\item{\code{distmat}}{Pairwise genetic distance data. Each row and column corresponds to patient ID in \code{patientseqIDs}.}
\item{\code{patientseqIDs}}{Vector of patient IDs corresponding to the rows and columns of \code{distmat}.}
  }
}
\details{
In the \code{epi} matrix, the source of infection is the patient ID of the true source of infection. A value of -1 indicates that the patient was infected prior to ward admission. Infection group denotes the genotype group to which this patient belongs. Multiple importations can belong to the same genotype group. Within-group genetic diversity is typically lower than between-group diversity. The group is assigned the ID of the first patient of the genotype group who enters the ward.
}
\examples{
data(hospitaldata)
head(hospitaldata$epi)
}
\keyword{datasets}
