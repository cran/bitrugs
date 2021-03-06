\name{simulate_data}
\alias{simulate_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nosocomial transmission data simulation
}
\description{
Simulation of hospital admission, discharge and disease transmission. Swab test results and genetic information are generated.
}
\usage{
simulate_data(n, D, LOS = 7, p = 0.05, z = 0.8, b = 0.005, gamma = 0.3, gamma_gl = 0.03, 
genpar = 0.8, testdays = 3, model = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of patient admissions to be generated.
}
  \item{D}{
Number of days over which patient admissions may be randomly distributed.
}
  \item{LOS}{
Mean length of hospital stay; lengths are drawn from a Poisson distribution.
}
  \item{p}{
Probability of carriage on admission.
}
  \item{z}{
Swab test sensitivity.
}
  \item{b}{
Transmission rate.
}
  \item{gamma}{
Within host/group genetic diversity (geometric distribution parameter).
}
  \item{gamma_gl}{
Between host/group genetic diversity (geometric distribution parameter).
}
  \item{genpar}{
Group clustering parameter (model 1), or transmission chain parameter (model 2). See details below.
}
  \item{testdays}{
Number of days between each carriage test/isolate collection.
}
  \item{model}{
Genetic diversity model (1: Importation clustering model, 2: Transmission chain model). See details below.
}
}
\details{
We implement two different models of genetic diversification. Model 1 is the importation clustering model, in which each infected host belongs to a group. The probability that a new importation belongs to a new group is \code{genpar}. All hosts within the same transmission chain belong to the same group. Pairwise genetic distances are sampled from geometric distributions with parameters \code{gamma} and \code{gamma_gl} for within-group and between-group pairs respectively. Model 2 is the transmission chain model, in which genetic distances increase as hosts are further separated in the transmission tree. Genetic distances between hosts separated by \code{k} transmission links are drawn independently from the geometric distribution with parameter \code{gamma} x \code{genpar}^\code{k}. Genetic distances between hosts in independent transmission chains are drawn from a geometric distribution with parameter \code{gamma_gl}.
}
\value{
Returns a list of simulated data;
\item{epi}{Matrix of epidemiological data, consisting of columns patient ID, day of admission, day of discharge, time of colonization, source of infection, and infection group.}
\item{resmat}{Matrix of test results. Each row represents the corresponding patient in \code{epi}. Each column represents a day, spanning the interval over which data were collected. 0=Negative, 1=Positive, -1=Missing}
\item{distmat}{Pairwise genetic distance data. Each row and column corresponds to patient ID in \code{patientseqIDs}.}
\item{patientseqIDs}{Vector of patient IDs corresponding to the rows and columns of \code{distmat}.}
}
\seealso{
\code{\link{simulate_data_dates}} to pre-specify admission and discharge dates.
}
\examples{
  Tsim <- simulate_data(n=100,D=50,p=0.1,b=0.01)
}
