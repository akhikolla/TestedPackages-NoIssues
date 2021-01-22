mycorMatrix <-
function(mycorcoefs,myDATA){
    myDATA <- data.frame(myDATA)
    index <- myDATA$index
    string <- names(mycorcoefs)
    substring <- substring(string,1,3)
    ar.ord <- length(substring[substring=="Phi"])
    ma.ord <- length(substring[substring=="The"])
    csARMA <- corARMA(mycorcoefs, form = ~index, p = ar.ord,q=ma.ord);
    csARMA <- Initialize(csARMA, data = myDATA);

    corMatrix(csARMA);
}
