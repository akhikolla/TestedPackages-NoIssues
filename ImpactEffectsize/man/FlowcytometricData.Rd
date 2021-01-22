\name{FlowcytometricData}
\alias{FlowcytometricData}

\docType{data}
\title{Example data of hematologic marker expression.}

\description{Data set of 8 flow cytometry-based lymphoma makers from 1,494 cells from healthy subjects (class 1) and 1,302 cells from lymphoma patients (class 2).}

\usage{data("FlowcytometricData")}

\details{
Size 2796 x 9 , Dimensions 1, stored in \code{FlowcytometricData$$[CD3,CD4,CD8,CD11,CD19,CD103,CD200,IgM]}

Classes 2, stored in \code{FlowcytometricData$Classes}
}

\examples{
data(FlowcytometricData)
str(FlowcytometricData)
}
