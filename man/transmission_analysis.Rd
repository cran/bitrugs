\name{transmission_analysis}
\alias{transmission_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate transmission dynamics and parameters
}
\description{
MCMC algorithm to sample transmission parameters, infection times and transmission routes.
}
\usage{
transmission_analysis(epidata, distmat, seqIDs, resmat, iterations = 1e+05, 
augmoves = 10, feedb = 1, tag = 1, noaug = 1, model = 2, path = NULL, 
sensprior = c(1, 1), impprior = c(1, 1), betaprior = 1e+06, 
gammaprior = c(1, 1), gammaGprior = c(1, 1), tparprior = 1e+06, 
sigma = c(0.004, 0.03, 0.005, 0.25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epidata}{
Epidemiological data, in the form of an integer matrix consisting of columns: patient ID, day of admission, day of discharge.
}
  \item{distmat}{
Genetic distance matrix. Entry [i,j] provides the pairwise SNP distance between patients \code{seqIDs[i]} and \code{seqIDs[j]}.
}
  \item{seqIDs}{
Vector of patient IDs corresponding to the rows and columns of \code{distmat}.
}
  \item{resmat}{
Matrix of test results, each row corresponding to the patient ID in \code{epidata}. Entry 0=negative, 1=positive, -1=missing.
}
  \item{iterations}{
Number of iterations for the MCMC algorithm to run.
}
  \item{augmoves}{
Number of data augmentation moves to make per MCMC iteration.
}
  \item{feedb}{
Frequency of console feedback; provided every 10^\code{feeback} iterations, with parameter snapshot every 10^\code{feeback+1} iterations.
}
  \item{tag}{
Integer tag to attach to output file.
}
  \item{noaug}{
Level of data augmentation. 0=none, 1=sample infection times and routes for patients with positive swabs only, 2=sample infection times and routes for all patients.
}
  \item{model}{
Genetic diversity model to use. 1=importation clustering model, 2=transmission chain diversity model.
}
  \item{path}{
Location to store output files.
}
  \item{sensprior}{
Prior Beta distribution parameters for test sensitivity (\code{z}).
}
  \item{impprior}{
Prior Beta distribution parameters for importation probability (\code{p}).
}
  \item{betaprior}{
Prior exponential distribution mean for transmission rate (\code{beta}).
}
  \item{gammaprior}{
Prior Beta distribution parameters for within host/group genetic diversity (\code{gamma}).
}
  \item{gammaGprior}{
Prior Beta distribution parameters for between host/group genetic diversity (\code{gamma_gl}).
}
  \item{tparprior}{
Prior exponential distribution mean for \code{genpar}.
}
  \item{sigma}{
Vector of initial variances for Normal proposal distributions for \code{beta}, \code{gamma}, \code{gamma_gl} and \code{genpar}. MCMC algorithm automatically updates variances to reach acceptance rate of 20-40\%.
}
}
\details{
MCMC algorithm runs in C, and writes output file to specified \code{path}.
}
\value{
Returns a matrix in which each row corresponds to MCMC iteration. Columns are as follows: p, z, beta, gamma, gamma_gl, genpar, number of importations, number of acquisitions, number of groups, likelihood, infection source [cols 11:(n+10)], infection group [cols (n+10):(2n+10)].
}
\examples{
  \dontrun{
  data(hospitaldata)
  # Short example run
  mcmcoutput <- transmission_analysis(epidata=hospitaldata$epi, distmat=hospitaldata$distmat, 
                           seqIDs=hospitaldata$patientseqIDs, resmat=hospitaldata$resmat, 
                           path=getwd(), iterations=10000, augmoves=5)
  traceplots(mcmcoutput)
  }
}
