\name{match.split}
\alias{match.split}
\title{
  Find Parallel Matches for Values in Groups (wordspace)
}
\description{

  Given a set of values and a grouped vector \code{x}, find parallel matches of each value in the different groups and return their positions in the original vector \code{x}.  If there are multiple matches of the same value in a group, only the position of the first match is returned.

}
\usage{

match.split(x, f, values=NULL, groups=NULL, nomatch=NA_integer_)

}
\arguments{

  \item{x}{vector to be divided into groups and matched against}

  \item{f}{a factor that defines the grouping (or a vector that can be converted to a factor)}

  \item{values}{values to be matched in \code{x}.
    Defaults to values that occur in all groups of \code{x} as determined by \code{f} and \code{groups}
  }

  \item{groups}{
    a character vector listing the set of groups to be formed.  Defaults to the levels of \code{f} and should be a subset of these levels if given explicitly
  }

  \item{nomatch}{
    the value to be returned in cases where no match is found (coerced to an integer)
  }
  
}
%% \details{
%% }
\value{

  An integer matrix with one row for each value (in \code{values}) and one column for each group (in \code{groups}), specifying the index in \code{x} of the first match of a value within the respective group.  If not match is found for a given combination of value and group, \code{nomatch} is inserted (defaults to \code{NA}).

}
%% \seealso{
%% }
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
