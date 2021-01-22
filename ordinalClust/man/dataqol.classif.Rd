\name{dataqol.classif}
\docType{data}
\alias{dataqol.classif}
\title{Questionnaire Responses Of Patients Affected By Breast Cancer}
\description{
 This dataset contains the responses of 40 patients to 30 questions
 about their quality of life. Furthermore, a variable indicates if the patient
 survived the disease.
}
\usage{dataqol.classif}

\format{ A dataframe with 40 lines and 32 columns. A line represents
 a patient and a column contains information about the patients.
	\describe{
	    \item{Id}{patient Id}
	    \item{q1-q28}{responses to 28 questions with the number of levels equal to 4}
	    \item{q29-q30}{responses to 22 questions with the number of levels equal to 7}
		\item{death}{if the patient survived (1) or not (2)}
	}
}
\source{The table was determined based on data associated with the package
 available at: \url{https://cran.r-project.org/package=QoLR}}
%%\references{
%% McNeil, D. R. (1977) \emph{Interactive Data Analysis}.
%% New York: Wiley.
%%}
\keyword{datsaets}