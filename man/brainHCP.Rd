\name{brainHCP}
\alias{brainHCP}
\docType{data}
\title{Time series obtained by an fMRI experiment on the brain}
\description{
  Time series for each region of interest in the brain. These series are obtained by SPM preprocessing.
}
\usage{data(brainHCP)}
\format{
  A data frame with 1200 observations on the following 89 variables.
}
\source{
contact S. Achard (sophie.achard@gipsa-lab.fr)
}
\references{
M. Termenon, A. Jaillard, C. Delon-Martin, S. Achard (2016)
Reliability of graph analysis of resting state fMRI using test-retest dataset from the Human Connectome Project, \emph{Neuroimage}, Vol 142, pages 172-187. 
 }
\examples{
data(brainHCP)
## maybe str(brainHCP) ; plot(brainHCP) ...
}
\keyword{datasets}
