

## ClusterR 1.2.2

* I've added the vectorized version of clusters to the output of the Affinity Propagation algorithm
* I've added the *threads* parameter to the *predict_KMeans()* function to return the k-means clusters in parallel  (useful especially for high dimensional data, see: https://stackoverflow.com/q/61551071/8302386)
* I've added a check-duplicated *CENTROIDS* if-condition in the *predict_KMeans()* function similar to the base kmeans function (see: https://stackoverflow.com/q/61551071/8302386). Due to the fact that the *CENTROIDS* output matrix is of class *"k-means clustering"* the base R function *duplicated()* performs a check column-wise rather than row-wise. Therefore before checking for duplicates I have to set the class to NULL.


## ClusterR 1.2.1

* I added a dockerfile in the root of the package directory and instructions in the README.md file on how to build and run the docker image  (https://github.com/mlampros/ClusterR/issues/17)
* I fixed a documentation and Vignette mistake regarding the *KMeans_rcpp* function (https://github.com/mlampros/ClusterR/issues/19)
* I fixed the *"failure: the condition has length > 1"* CRAN error which appeared mainly due to the misuse of the base *class()* function in multiple code snippets in the package (for more info on this matter see: https://developer.r-project.org/Blog/public/2019/11/09/when-you-think-class.-think-again/index.html)


## ClusterR 1.2.0

* I added the 'cosine' distance to the following functions: 'Cluster_Medoids', 'Clara_Medoids', 'predict_Medoids', 'Optimal_Clusters_Medoids' and 'distance_matrix'.
* I fixed an error case in the .pdf manual of the package (https://github.com/mlampros/ClusterR/issues/16)


## ClusterR 1.1.9

* I added parallelization for the *exact* method of the *AP_preferenceRange* function which is more computationally intensive as the *bound* method
* I modified the *Optimal_Clusters_KMeans*, *Optimal_Clusters_GMM* and *Optimal_Clusters_Medoids* to accept also a contiguous or non-contiguous vector besides single values as a *max_clusters* parameter. However, the limitation currently is that the user won't be in place to plot the clusters but only to receive the ouput data ( this can be changed in the future however the plotting function for the contiguous and non-contiguous vectors must be a separate plotting function outside of the existing one).  Moreover, the *distortion_fK* criterion can't be computed in the *Optimal_Clusters_KMeans* function if the *max_clusters* parameter is a contiguous or non-continguous vector ( the *distortion_fK* criterion requires consecutive clusters ). The same applies also to the *Adjusted_Rsquared* criterion which returns incorrect output. For this feature request see the following [Github issue](https://github.com/mlampros/ClusterR/issues/15).


## ClusterR 1.1.8

* I moved the *OpenImageR* dependency in the DESCRIPTION file from 'Imports' to 'Suggests', as it appears only in the Vignette file.


## ClusterR 1.1.7

* I fixed the *clang-UBSAN* errors


## ClusterR 1.1.6

* I updated the README.md file (I removed unnecessary calls of ClusterR in DESCRIPTION and NAMESPACE files)
* I renamed the *export_inst_header.cpp* file in the src folder to *export_inst_folder_headers.cpp*
* I modified the *Predict_mini_batch_kmeans()* function to accept an armadillo matrix rather than an Rcpp Numeric matrix. The function appers both in *ClusterRHeader.h* file ( 'inst' folder ) and in *export_inst_folder_headers.cpp* file ( 'src' folder )
* I added the *mini_batch_params* parameter to the *Optimal_Clusters_KMeans* function. Now, the optimal number of clusters can be found also based on the min-batch-kmeans algorithm (except for the *variance_explained* criterion)
* I changed the license from MIT to GPL-3
* I added the *affinity propagation algorithm* (<span></span>www.psi.toronto.edu/index.php?q=affinity%20propagation). Especially, I converted the matlab files *apcluster.m* and *referenceRange.m*.
* I modified the minimum version of RcppArmadillo in the DESCRIPTION file to 0.9.1 because the Affinity Propagation algorithm requires the *.is_symmetric()* function, which was included in version 0.9.1


## ClusterR 1.1.5

As of version 1.1.5 the ClusterR functions can take [tibble](https://tibble.tidyverse.org/) objects as input too.


## ClusterR 1.1.4

I modified the ClusterR package to a cpp-header-only package to allow linking of cpp code between Rcpp packages. See the update of the README.md file (16-08-2018) for more information.


## ClusterR 1.1.3

I updated the example section of the documentation by replacing the *optimal_init* with the *kmeans++* initializer


## ClusterR 1.1.2

* I fixed an [Issue](https://github.com/mlampros/ClusterR/issues/8) related to *NAs produced by integer overflow* of the *external_validation* function. See, the commented line of the *Clustering_functions.R* file (line 1830).


## ClusterR 1.1.1

* I added a *tryCatch* in *Optimal_Clusters_Medoids()* function to account for the error described in [Error in Optimal_Clusters_Medoids function#5](https://github.com/mlampros/ClusterR/issues/5) issue


## ClusterR 1.1.0

* I added the *DARMA_64BIT_WORD* flag in the Makevars file to allow the package processing big datasets
* I modified the *kmeans_miniBatchKmeans_GMM_Medoids.cpp* file and especially all *Rcpp::List::create()* objects to addrress the clang-ASAN errors.


## ClusterR 1.0.9

* I modified the *Optimal_Clusters_KMeans* function to return a vector with the *distortion_fK* values if criterion is *distortion_fK* (instead of the *WCSSE* values).
* I added the 'Moore-Penrose pseudo-inverse' for the case of the 'mahalanobis' distance calculation.


## ClusterR 1.0.8

* I modified the *OpenMP* clauses of the .cpp files to address the ASAN errors.
* I removed the *threads* parameter from the *KMeans_rcpp* function, to address the ASAN errors ( negligible performance difference between threaded and non-threaded version especially if the *num_init* parameter is less than 10 ). The *threads* parameter was removed also from the *Optimal_Clusters_KMeans* function as it utilizes the *KMeans_rcpp* function to find the optimal clusters for the various methods.


## ClusterR 1.0.7

I modified the *kmeans_miniBatchKmeans_GMM_Medoids.cpp* file in the following lines in order to fix the clang-ASAN errors (without loss in performance):

* lines 1156-1160 : I commented the second OpenMp parallel-loop and I replaced the *k* variable with the *i* variable in the second for-loop [in the *dissim_mat()* function]
* lines 1739-1741 : I commented the second OpenMp parallel-loop [in the *silhouette_matrix()* function]
* I replaced (all) the *silhouette_matrix* (arma::mat) variable names with *Silhouette_matrix*, because the name overlapped with the name of the Rcpp function [in the *silhouette_matrix* function]
* I replaced all *sorted_medoids.n_elem* with the variable *unsigned int sorted_medoids_elem* [in the *silhouette_matrix* function]

I modified the following *functions* in the *clustering_functions.R* file:

* *KMeans_rcpp()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.
* *Optimal_Clusters_KMeans()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.
* *MiniBatchKmeans()* : I added an *experimental* note in the details for the *optimal_init* and *quantile_init* initializers.


## ClusterR 1.0.6

The *normalized variation of information* was added in the *external_validation* function (https://github.com/mlampros/ClusterR/pull/1)


## ClusterR 1.0.5

I fixed the valgrind memory errors


## ClusterR 1.0.4

I removed the warnings, which occured during compilation.
I corrected the UBSAN memory errors which occured due to a mistake in the *check_medoids()* function of the *utils_rcpp.cpp* file.
I also modified the *quantile_init_rcpp()* function of the *utils_rcpp.cpp* file to print a warning if duplicates are present in the initial centroid matrix.


## ClusterR 1.0.3

* I updated the dissimilarity functions to accept data with missing values.
* I added an error exception in the predict_GMM() function in case that the determinant is equal to zero. The latter is possible if the data includes highly correlated variables or variables with low variance.
* I replaced all unsigned int's in the rcpp files with int data types


## ClusterR 1.0.2

I modified the RcppArmadillo functions so that ClusterR passes the Windows and OSX OS package check results


## ClusterR 1.0.1

I modified the RcppArmadillo functions so that ClusterR passes the Windows and OSX OS package check results


## ClusterR 1.0.0




