\docType{data}
\name{RandomVA1}
\alias{RandomVA1}
\title{100 records of Sample Input}
\format{100 arbitrary input records.}
\description{
This is a  dataset consisting of 100 arbitrary sample input 
deaths in the acceptable format of InSilicoVA. It is in the
format of a list. The first element of the list is the matrix
data input file in InterVA format, while the second is a vector
of length 100 specifying subpopulation for the 100 deaths.
}
\examples{
data(RandomVA1)
dim(RandomVA1$data)
dim(RandomVA1$subpop)
}
\keyword{datasets}
