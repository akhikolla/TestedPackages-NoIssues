\name{Mc}
\alias{Mc}
\docType{data}
\title{McDonald's data}
\description{
Nutrition analysis of McDonald's menu items.
}
\usage{data(Mc)}
\format{
 A data.frame with 81 rows and 16 columns.
}
\details{
Data are from McDonald's USA Nutrition Facts for Popular Menu Items. A subset of menu items is reported. Beverages are excluded.
In case of duplications, regular size or medium size information is reported. The variable Type is a factor the levels of which specify the kind of the menu items.
Although some menu items could be well described by more than one level, only one level of the variable Type specifies each menu item.
Percent Daily Values (\%DV) are based on a 2,000 calorie diet. Some menu items are registered trademarks.
}
\author{Paolo Giordani, Maria Brigida Ferraro, Alessio Serafini}
\seealso{\code{\link{Fclust}}, \code{\link{FKM}}, \code{\link{FKM.ent}}, \code{\link{FKM.med}}}
\examples{
## McDonald's data
data(Mc)
names(Mc)
## data normalization by dividing the nutrition facts by the Serving Size (column 1)
for (j in 2:(ncol(Mc)-1))
Mc[,j]=Mc[,j]/Mc[,1]
## removing the column Serving Size
Mc=Mc[,-1]
p=(ncol(Mc)-1)
## fuzzy k-means (excluded the factor column Type (last column))
clust.FKM=FKM(Mc[,1:p],k=6,m=1.5,stand=1)
## new factor column Cluster.FKM containing the cluster assignment information
## using fuzzy k-means
Mc[,ncol(Mc)+1]=factor(clust.FKM$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM")
levels(Mc$Cluster.FKM)=paste("Clus FKM",1:clust.FKM$k,sep=" ")
## contingency table (Cluster.FKM vs Type)
## to assess whether clusters can be interpreted in terms of the levels of Type
table(Mc$Type,Mc$Cluster.FKM)
## prototypes using the original units of measurement
clust.FKM$Hraw=Hraw(clust.FKM$X,clust.FKM$H)
clust.FKM$Hraw
## fuzzy k-means with entropy regularization
## (excluded the factor column Type (last column))
\dontrun{
## It may take more than a few seconds
clust.FKM.ent=FKM.ent(Mc[,1:p],k=6,ent=3,RS=10,stand=1)
## new factor column Cluster.FKM.ent containing the cluster assignment information
## using fuzzy k-medoids with entropy regularization
Mc[,ncol(Mc)+1]=factor(clust.FKM.ent$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM.ent")
levels(Mc$Cluster.FKM.ent)=paste("Clus FKM.ent",1:clust.FKM.ent$k,sep=" ")
## contingency table (Cluster.FKM.ent vs Type)
## to assess whether clusters can be interpreted in terms of the levels of Type
table(Mc$Type,Mc$Cluster.FKM.ent)
## prototypes using the original units of measurement
clust.FKM.ent$Hraw=Hraw(clust.FKM.ent$X,clust.FKM.ent$H)
clust.FKM.ent$Hraw}
## fuzzy k-medoids
## (excluded the factor column Type (last column))
clust.FKM.med=FKM.med(Mc[,1:p],k=6,m=1.1,RS=10,stand=1)
## new factor column Cluster.FKM.med containing the cluster assignment information
## using fuzzy k-medoids with entropy regularization
Mc[,ncol(Mc)+1]=factor(clust.FKM.med$clus[,1])
colnames(Mc)[ncol(Mc)]=("Cluster.FKM.med")
levels(Mc$Cluster.FKM.med)=paste("Clus FKM.med",1:clust.FKM.med$k,sep=" ")
## contingency table (Cluster.FKM.med vs Type)
## to assess whether clusters can be interpreted in terms of the levels of Type
table(Mc$Type,Mc$Cluster.FKM.med)
## prototypes using the original units of measurement
clust.FKM.med$Hraw=Hraw(clust.FKM.med$X,clust.FKM.med$H)
clust.FKM.med$Hraw
## or, equivalently,
Mc[clust.FKM.med$medoid,1:p]}

\keyword{data}

\keyword{multivariate}
