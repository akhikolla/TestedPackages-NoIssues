# function to compute p-value
PvalueFam <- function(U) {
  stat <- colSums(U) * colSums(U)/colSums(U^2)
  p.value <- pchisq(stat, df = 1, lower.tail = FALSE)
  if(length(p.value)>10){
		FDR <- qvalue(p.value)$qvalues
		FWER <- p.adjust(p.value, "bonferroni")
		res <- data.frame(stat = stat, pvalue = p.value, FDR=FDR, FWER=FWER)}
	else{
		cat("pi0 is can not be evaluated due to small number of p-values\nFDR and FWER is not calculated.\n")
		res <- data.frame(stat = stat, pvalue = p.value)}
	res
}
