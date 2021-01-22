convert.lemma <- 
function (lemma, format=c("CWB", "BNC", "DM", "HW", "HWLC"), hw.tolower=FALSE) {
  format <- match.arg(format)
  lemma <- as.character(lemma)

  errors <- grep("_[A-Z.]$", lemma, invert=TRUE, perl=TRUE, value=TRUE)
  if (length(errors) > 0) stop("invalid lemma string(s): ", paste(errors, collapse=", "))

  len <- nchar(lemma)
  hw <- substr(lemma, 1, len - 2) # headword part
  pos <- substr(lemma, len, len)  # POS code
  
  if (format == "BNC") {
    penn2bnc <- c(N="SUBST", Z="SUBST", V="VERB", J="ADJ", R="ADV", D="ART", .="STOP", I="PREP", T="PREP")
    bncpos <- penn2bnc[pos]
    bncpos[is.na(bncpos)] <- "UNC"
    #### old code (should no longer be needed)
    ## penn <- c("N",     "Z",     "V",    "J",   "R",   "D",   ".",    "I",    "T")
    ## bnc  <- c("SUBST", "SUBST", "VERB", "ADJ", "ADV", "ART", "STOP", "PREP", "PREP", "UNC")
    ## pos.idx <- match(pos, penn, nomatch=length(bnc))
    ## bncpos <- bnc[pos.idx]
  }
  if (hw.tolower) {
    hw <- tolower(hw)
  }
  
  switch(format,
    CWB = paste(hw, pos, sep="_"),
    BNC = paste(tolower(hw), bncpos, sep="_"),
    DM  = paste(hw, tolower(pos), sep="-"),
    HW  = hw,
    HWLC = tolower(hw)
  )
}
