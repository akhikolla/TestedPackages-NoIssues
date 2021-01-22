#'@import ggplot2 mclust scales tools
NULL

#'EMMIXgene: 
#'
#'Selects genes using the EMMIXgene algorithm, following the methodology of 
#'G. J. McLachlan, R. W. Bean, D. Peel; A mixture model-based approach to the 
#'clustering of microarray expression data , Bioinformatics, Volume 18, 
#'Issue 3, 1 March 2002, Pages 413–422,
#'https://doi.org/10.1093/bioinformatics/18.3.413
#'
#'@section Functions:
#'\code{\link{select_genes}}: Selects the most differentially expressed genes.
#'
#'\code{\link{cluster_genes}}: Clusters the genes using a mixture model
#'approach.
#'
#'\code{\link{cluster_tissues}}: Clusters the tissues based on the differences 
#'between the tissue samples among the gene groups.
#'
#'See \code{vignette('The-EMMIXgene-Workflow')} for more details.
#' @docType package
#' @name EMMIXgene
NULL


#'@title Selects genes using the EMMIXgene algorithm.
#'
#'@description Follows the gene selection methodology of 
#'G. J. McLachlan, R. W. Bean, D. Peel; A mixture model-based approach to the 
#'clustering of microarray expression data , Bioinformatics, Volume 18, 
#'Issue 3, 1 March 2002, Pages 413–422,
#'https://doi.org/10.1093/bioinformatics/18.3.413
#'
#'@param dat A matrix or dataframe containing gene expression data.
#'Rows are genes and columns are samples. Must supply one of filename and dat.
#'@param filename Name of file containing gene data. Can be either .csv 
#'or space separated .dat. Rows are genes and columns are samples. 
#'Must supply one of filename and dat.
#'@param random_starts The number of random initializations used per gene when
#'fitting mixtures of t-distributions. Initialization uses k-means by default.
#'@param max_it The maximum number of iterations per mixture fit. 
#'Default value is 100.
#'@param ll_thresh The difference in -2 log lambda used as a threshold 
#'for selecting between g=1 and g=2 for each gene. Default value is 8,
#'which was chosen arbitrarily in the original paper.
#'@param min_clust_size The minimum number of observations per cluster
#'used when fitting mixtures of t-distributions for each gene.
#'Default value is 8. 
#'@param tol Tolerance value used for detecting convergence of EMMIX fits. 
#'@param start_method Default value is "both". 
#'Can also choose "random" for purely random starts.
#'@param three Also test g=2 vs g=3 where appropriate. Defaults to FALSE.
#'@return An EMMIXgene object containing:
#'\item{stat}{The difference in log-likelihood for g=1 
#'and g=2 for each gene (or for g=2 and g=3 where relevant).}
#'\item{g}{The selected number of components for each gene.}
#'\item{it}{The number of iterations for each genes selected fit.} 
#'\item{selected}{An indicator for each genes selected status}
#'\item{ranks}{selected gene ids ranked by stat}
#'\item{genes}{A dataframe of selected genes.}
#'\item{all_genes}{Returns dat or contents of filename.}
#'
#'@examples
#'#only run on first 100 genes for speed
#'alon_sel <- select_genes(alon_data[seq_len(100), ]) 
#'
#'@export
select_genes<-function(dat, filename, random_starts=4, max_it = 100, 
                ll_thresh = 8, min_clust_size = 8, tol = 0.0001,
                start_method ="both",  three=FALSE){
    #housekeeping   
    
    if(missing(filename)&missing(dat)){
        stop("Must supply one of filename and dat.")
    }
    
    if(!missing(filename)&!missing(dat)){
        stop("Must supply ONLY one of filename and dat.")
    }  
    
    if(!missing(filename)&missing(dat)){
        if(tools::file_ext(filename)=="dat"){
            data<-as.matrix(utils::read.delim(filename, sep=" ", header=FALSE))}
        if(tools::file_ext(filename)=="csv"){
            data<-as.matrix(utils::read.csv(filename, header=FALSE))}
    }
    
    if(missing(filename)&!missing(dat)){
        data<-as.matrix(dat)
    }
    
    
    if(any(!stats::complete.cases((data)))| 
        any(!stats::complete.cases(t(data)))){
        warning("Incomplete cases removed removed.")
    }
    
    if(any(apply(data,1,stats::var)==0)){
        warning("Some rows are constant across samples and have been removed.")
    }
    
    data<-data[!apply(data,1,stats::var)==0,]
    
    
    #remove missing data
    data<-data[stats::complete.cases((data)),]
    data<-data[,stats::complete.cases(t(data))]
    
    
    #actual method
    
    #do emmix_gene C++ routine
    a<-emmix_gene(data,random_starts, max_it, ll_thresh,
            min_clust_size,tol,start_method, three)
    a$selected<-a$g>1

    #add selected genes to result
    a$genes <- data[a$g>1,]
    a$all_genes <- data
    
    #ranks of selected 
    a$stat2<-a$stat
    a$stat2[!a$selected]<-NA
    a$ranks<-order(a$stat2, decreasing = TRUE, na.last = TRUE)
    a$stat2 <- NULL
    
    result<-structure(a, class="EMMIXgene")
    #return selected genes
    return(result)
}

#'Clusters genes using mixtures of normal distributions
#'
#'Sorts genes into clusters using mixtures of normal distributions with
#'covariance matrices restricted to be multiples of the identity matrix.
#'
#'@param gen an EMMIXgene object produced by select_genes().  
#'@param g The desired number of gene clusters. If not specified will be 
#'selected automatically on the basis of BIC.
#'@return An array containing the clustering. 
#'@examples
#'
#'#only run on first 100 genes for speed
#'alon_sel <- select_genes(alon_data[seq_len(100), ]) 
#'alon_clust<- cluster_genes(alon_sel , 2)
#'@export
cluster_genes<-function(gen, g=NULL){
    
    clust_genes_k<-tkmeans(as.matrix(gen$genes), k=g, 0.00001, 
            rep(1, ncol(gen$genes)), 10, 100, 0.001, FALSE)
    clusters<-nearest_cluster(as.matrix(gen$genes), clust_genes_k)
    #clust_genes<-mclust::Mclust((gen$genes), G=g, modelNames = "VII")
    #g<-clust_genes$G
    
    ll_rank_stat<-array(0,g)
    
    for(i in seq_len(g)){
        ll_rank_stat[i]<-
            each_gene(colMeans(gen$genes[clusters==i, ,drop=FALSE]))$Ratio
    }
    
    #reorder clusters as ranked by mean of ll statistic from EMMIXgene fit
    tmp<-factor(clusters)
    levels(tmp) <- as.character(order(unlist(ll_rank_stat), decreasing = TRUE))
    clust_genes<-as.numeric(as.character(tmp))
    
    return(clust_genes)
}


#'Clusters tissues using all group means
#'
#'@param gen EMMIXgene object
#'@param clusters mclust object
#'@param q number of factors if using mfa
#'@param G number of components if using mfa
#'@return a clustering for each sample (columns) by each group(rows)
#'@examples
#'
#'example <- plot_single_gene(alon_data,1) 
#'#only run on first 100 genes for speed
#'alon_sel <- select_genes(alon_data[seq_len(100), ]) 
#'alon_clust<- cluster_genes(alon_sel , 2)
#'\donttest{alon_tissue_all<-all_cluster_tissues(alon_sel, alon_clust, q=1, G=2)}
#'@export
all_cluster_tissues<-function(gen, clusters, q=6, G=2){
    g<- length(table(clusters))
    
    p<-ncol(gen$genes)
    group_means<-array(0,c(g,p))
    
    for(i in seq_len(g)){
        group_means[i,] <- colMeans(gen$genes[clusters==i,,drop=FALSE])
    }
    mfa_fit<-mcfa(t(group_means), G, q, itmax=100, nkmeans=5, nrandom=5)
    clustering<- as.numeric(predict.mcfa(mfa_fit, t(group_means))-1) 
    
    return(clustering)
}

#'Clusters tissues
#'
#'@param gen EMMIXgene object
#'@param clusters mclust object
#'@param method Method for separating tissue classes. Can be either 't' for a 
#'univariate mixture of t-distributions on gene cluster means, or 'mfa'
#'for a mixture of factor analyzers. 
#'@param q number of factors if using mfa
#'@param G number of components if using mfa
#'@return a clustering for each sample (columns) by each group(rows)
#'@examples
#'#only run on first 100 genes for speed
#'alon_sel <- select_genes(alon_data[seq_len(100), ]) 
#'alon_clust<- cluster_genes(alon_sel,2)
#'alon_tissue_t<-
#'    cluster_tissues(alon_sel,alon_clust,method='t')
#'alon_tissue_mfa<-
#'    cluster_tissues(alon_sel, alon_clust,method='mfa',q=2,G=2) 
#'@export
cluster_tissues<-function(gen, clusters, method='t', q=6, G=2){
    g<-length(table(clusters))
    p<-ncol(gen$genes)
    
    clustering<-array(0,c(g,p))
    clustering2<-array(0,c(g,p))
    
    if(method=='t'){
        for(i in seq_len(g)){
            group_means <- colMeans(gen$genes[clusters==i,,drop=FALSE])
            t_fit<-emmix_t(group_means, G)
            clustering[i,]<-t_fit$Clusters
        }
        
    }
    
    if(method=='mfa'){
        for(i in seq_len(g)){
            
            group <- as.matrix((gen$genes[clusters==i,,drop=FALSE]))
            #actually mixture of common factor analysers. consider fixing.
            if(dim(group)[1]>20){
                q1<-min(q, sum(clusters==i)-1 )
                mfa_fit<-mcfa(t(group), G, q1, itmax=100,
                nkmeans=50, nrandom=50)
                clustering[i,]<- as.numeric(predict.mcfa(mfa_fit, t(group))-1) 
            }
            
            
        }
        
        
    }
    
    return(clustering)
}




#'Cluster tissues
#'
#'@param gen An EMMIXgene object produced by select_genes().
#'@param n_top number of top genes (as ranked by likelihood) to be selected
#'@param method Method for separating tissue classes. Can be either 't' for a 
#'univariate mixture of t-distributions on gene cluster means, 
#'or 'mfa' for a mixture of factor analysers. 
#'@param q number of factors if using mfa
#'@param g number of components if using mfa
#'@return An EMMIXgene object containing:
#'\item{stat}{A matrix containing clustering (0 or 1) 
#'for each sample (columns) by each group(rows).}
#'\item{top_gene}{The row numbers of the top genes.}
#'\item{fit}{The fit object used to determine the clustering.}
#'@examples
#'\donttest{alon_sel <- select_genes(alon_data[seq_len(100), ])}
#'\donttest{alon_top_10<-top_genes_cluster_tissues(alon_sel, 10, method='mfa', q=3, g=2)}
#'@export
top_genes_cluster_tissues<-function(gen, n_top=100, method='mfa', q=2, g=2){
    
    
    p<-ncol(gen$genes)
    clustering<-array(0,p)
    
    top_genes<-gen$ranks[seq_len(n_top)]
    
    if(method=='t'){
        
        group_means <- colMeans(gen$all_genes[top_genes,])
        fit<-emmix_t(group_means, g)
        clustering<-fit$Clusters
        
        
    }
    
    if(method=='mfa'){
        
        group <- as.matrix((gen$all_genes[top_genes,]))
        fit<-mcfa(t(group), g, q)
        clustering<- predict.mcfa(fit, t(group))-1
        
        
        
    }
    
    
    
    
    return(list(clustering=clustering, top_genes=top_genes, fit=fit))
}










#'Heat maps
#'
#'Plot heat maps of gene expression data. Optionally sort the x-axis 
#'according to a predetermined clustering.
#'
#'@param dat matrix of gene expression data.
#'@param clustering a vector of sample classifications.
#'Must be same length as the number of columns in dat.
#'@param y_lab optional label for y-axis.
#'@return A ggplot2 heat map.
#'@examples
#' 
#'example <- heat_maps(alon_data[seq_len(100), ])
#'
#'
#'@export
heat_maps<-function(dat, clustering=NULL, y_lab=NULL){
    colnames(dat) <- NULL
    if(!is.null(clustering)){
        dat<-dat[,order(clustering)]
    }
    
    
    df_heatmap<-reshape::melt(dat)
    names(df_heatmap)<-c("genes", "samples",  "expression_level")
    df_heatmap$genes<-factor(df_heatmap$genes)
    df_heatmap$samples<-factor(df_heatmap$samples)
    
    plot<-ggplot(df_heatmap, aes(df_heatmap$samples,df_heatmap$genes ))+
        geom_tile(aes(fill = df_heatmap$expression_level),  color = "white")+
        scale_fill_distiller(palette = "Spectral")+  
        xlab("Samples")  + guides(fill=guide_legend(title="Expression Level"))+
        theme(axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.ticks.y=element_blank()
        )
    
    if(!is.null(clustering)){
        plot<- plot + ylab(y_lab)
    }else{
        plot<- plot + ylab("Genes")
    }
    
    if(is.null(clustering)){
        plot<- plot +  
        scale_x_discrete(breaks = seq(0, length(levels(df_heatmap$genes)), 10) )
    }
    
    if(!is.null(clustering)){
        plot<- plot  +  scale_x_discrete(labels =order(clustering)) 
        #+ theme(axis.text.x = element_blank())
    }
    
    return(plot)
}

#'Plot a single gene expression histogram with best 
#'fitted mixture of t-distributions.
#'
#'Plot a single gene expression histogram with best 
#'fitted mixture of t-distributions according to the EMMIX-gene algorithm.
#'
#'@param dat matrix of gene expression data.
#'@param gene_id row number of gene to be plotted.
#'@param g force number of components, default = NULL
#'@param random_starts The number of random initializations used 
#'per gene when fitting mixtures of t-distributions. 
#'Initialization uses k-means by default.
#'@param max_it The maximum number of iterations per mixture fit.
#'Default value is 100.
#'@param ll_thresh The difference in -2 log lambda used as a 
#'threshold for selecting between g=1 and g=2 for each gene. 
#'Default value is 8, which was chosen arbitrarily in the original paper.
#'@param min_clust_size The minimum number of observations per cluster 
#'used when fitting mixtures of t-distributions for each gene.
#'Default value is 8. 
#'@param tol Tolerance value used for detecting convergence of EMMIX fits. 
#'@param start_method Default value is "both". 
#'Can also choose "random" for purely random starts.
#'@param three Also test g=2 vs g=3 where appropriate. Defaults to TRUE.
#'@param min,max Minimum and maximum x-axis values for the plot window.
#'@return A ggplot2 histogram with fitted t-distributions overlayed. 
#'@examples
#'
#'example <- plot_single_gene(alon_data,1) 
#'#plot(example)
#'
#'@export
plot_single_gene<-function(dat, gene_id, g=NULL, 
                           random_starts=8, max_it = 100,
                           ll_thresh = 8, min_clust_size = 8,
                           tol = 0.0001, start_method = "both",
                           three=TRUE, min = -4, max = 2){ 
    
    df<-data.frame(x=dat[gene_id,])
    n<-length(df$x)/4
    breaks<-seq(min, max, length.out=n)
    p<-ggplot(df, aes(x=df$x)) + geom_histogram(breaks=breaks, alpha=.5)
    p2<-ggplot_build(p+theme_bw())
    
    df_dens<-data.frame(x=p2$data[[1]]$x, y=p2$data[[1]]$density) 
    
    plot<-ggplot(df_dens)
    plot<-plot+geom_bar(aes(x=df_dens$x, y=df_dens$y), alpha=.5,
        stat="identity", position="identity")
    plot<-plot+theme_bw()+xlab("Gene Expression Value")+ylab("Density")
    
    df2<-data.frame(x=seq(min, max, length.out = 1000))
    
    if(!is.null(g)){
        res<-emmix_t(dat[gene_id,], g, random_starts, max_it, tol,start_method)
    }else{
        res<-each_gene(dat[gene_id,],random_starts, max_it, ll_thresh, 
            min_clust_size, tol,start_method,three)
    }
    
    
    for(i in seq_len(res$components)){
        for(j in seq_len(nrow(df2))){
            df2[[paste0('y',i)]][j]<-res$pi[i]*t_dist(df2$x[j], res$mu[i], 
                res$sigma[i], res$nu[i])
        }
    }
    
    
    plot<-plot+geom_line(data=df2, aes(x=df2$x, y=df2$y1))
    
    if(res$components>1){
        plot<-plot+geom_line(data=df2, aes(x=df2$x, y=df2$y2))
    }
    
    if(res$components>2){
        plot<-plot+geom_line(data=df2, aes(x=df2$x, y=df2$y3))
    }
    
    return(plot)
}

#'@title Normalized gene expression values from Alon et al. (1999).
#'
#'@description A dataset containing centred and normalized values of the 
#'logged expression values of a subset of 2000 genes taken from 
#'Alon, Uri, et al. "Broad patterns of gene expression revealed by clustering
#'analysis of tumor and normal colon tissues probed by oligonucleotide arrays."
#'Proceedings of the National Academy of Sciences 96.12 (1999): 6745-6750.
#'The method of subset selection was described in G. J. McLachlan, R. W. Bean, 
#'D. Peel; A mixture model-based approach to the clustering of microarray 
#'expression data ,
#'Bioinformatics, Volume 18, Issue 3, 1 March 2002, Pages 413–422.
#'
#'@docType data
#'@keywords datasets
#'@name alon_data
#'@usage data(alon_data)
#'@format A data frame with 2000 rows (genes) and 62 variables (samples).
#'@examples
#'dim(alon_data)
NULL

#'@title Normalized gene expression values from Golub et al. (1999).
#'
#'@description A dataset containing the centred and normalized values of the
#'logged expression values of a subset of 3731 genes taken from Golub,
#'Todd R., et al. "Molecular classification of cancer: class discovery
#'and class prediction by gene expression monitoring." 
#'Science 286.5439 (1999): 531-537.
#'The method of subset selection was described in G. J. McLachlan, R. W. Bean,
#'D. Peel; A mixture model-based approach to the clustering of microarray
#'expression data ,
#'Bioinformatics, Volume 18, Issue 3, 1 March 2002, Pages 413–422.
#'
#'@docType data
#'@keywords datasets
#'@name golub_data
#'@usage data(golub_data)
#'@format A data frame with 3731 rows (genes) and 72 variables (samples).
#'#'@examples
#'dim(golub_data)
NULL
