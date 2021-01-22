\name{convert.lemma}
\alias{convert.lemma}
\title{

  Transform CWB/Penn-Style Lemmas into Other Notation Formats (wordspace)

}
\description{

  Transform POS-disambiguated lemma strings in CWB/Penn format (see Details)
  into several other notation formats.

}
\usage{

convert.lemma(lemma, format=c("CWB", "BNC", "DM", "HW", "HWLC"), hw.tolower=FALSE)

}
\arguments{

  \item{lemma}{a character vector specifying one or more POS-disambiguated lemmas in CWB/Penn notation}

  \item{format}{the notation format to be generated (see Details)}

  \item{hw.tolower}{convert headword part to lowercase, regardless of output format}
  
}
\details{
  
  Input strings must be POS-disambiguated lemmas in CWB/Penn notation, i.e. in the form
  \preformatted{
    <headword>_<P> }
  where \code{<headword>} is a dictionary headword (which may be case-sensitive) and \code{<P>} is
  a one-letter code specifying the simple part of speech.  Standard POS codes are
  \preformatted{
    N ... nouns
    Z ... proper nouns
    V ... lexical and auxiliary verbs
    J ... adjectives
    R ... adverbs
    I ... prepositions (including all uses of "to")
    D ... determiners
    . ... punctuation
  }
  For other parts of speech, the first character of the corresponding Penn tag may be used.
  Note that these codes are not standardised and are only useful for distinguishing between content
  words and function words.
  
  The following output formats are supported:
  \describe{

    \item{\code{CWB}}{
      returns input strings without modifications, but validates that they are in CWB/Penn format
    }

    \item{\code{BNC}}{
      BNC-style POS-disambiguated lemmas based on the simplified CLAWS tagset.
      The headword part of the lemma is unconditionally converted to lowercase.
      The standard POS codes listed above are translated into
      \code{SUBST} (nouns and proper nouns), \code{VERB} (verbs), \code{ADJ} (adjectives), \code{ADV} (adverbs),
      \code{ART} (determiners), \code{PREP} (prepositions), and \code{STOP} (punctuation).
      Other POS codes have no direct CLAWS equivalents and are mapped to \code{UNC} (unclassified),
      so the transformation should only be used for the categories listed above.
    }

    \item{\code{DM}}{
      POS-disambiguated lemmas in the format used by Distributional Memory (Baroni & Lenci 2010),
      viz. \code{<headword>-<p>} with POS code in lowercase and headword in its original capitalisation.
      For example, \code{light_N} will be mapped to \code{light-n}.
    }

    \item{\code{HW}}{
      just the undisambiguated headword
    }

    \item{\code{HWLC}}{
      undisambiguated headword mapped to lowercase (same as \code{HW} with \code{hw.tolower=TRUE})
    }
  }

}
\value{
  
  A character vector of the same length as \code{lemma}, containing the transformed lemmas.
  See Details above for the different output formats.

}
\references{
  
  Baroni, Marco and Lenci, Alessandro (2010).
  Distributional Memory: A general framework for corpus-based semantics.
  \emph{Computational Linguistics}, \bold{36}(4), 673--712.
  
}
\examples{

convert.lemma(RG65$word1, "CWB") # original format
convert.lemma(RG65$word1, "BNC") # BNC-style (simple CLAWS tags)
convert.lemma(RG65$word1, "DM")  # as in Distributional Memory
convert.lemma(RG65$word1, "HW")  # just the headword

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
