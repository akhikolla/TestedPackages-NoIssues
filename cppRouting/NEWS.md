
cppRouting v2.0
===============

Major changes

-   implementation of contraction hierarchies algorithm
-   thread-safe implementation of all parallel algorithms with `RcppParallel` package instead of `parallel`
-   implementation of one-to-one query algorithm on contracted graph 
-   implementation of many-to-many query algorithm on contracted graph 
-   implementation of PHAST algorithm on contracted graph 

Minor changes

-   new options `long` and `keep` for `get_path_pair`, `get_isochrone` and `get_multi_paths` functions    
-   remove `allcores` option for `get_detour` function  
-   optimization of `cpp_simplify` function  
-   optimization of `makegraph` function

cppRouting v1.2
===============

Major changes

-   new functions `cpp_simplify`, `get_detour` and `to_df`

Minor changes

-   bug fix in `get_distance_pair` and `get_path_pair` : verify that origin / destination nodes are present in the graph before running c++ function
-   modification of `get_distance_pair` and `get_path_pair` Rd. files : clearer explanation about the choice between algorithms
-   optimization of `get_distance_matrix` function : if length(to) &lt; length(from) then Dijkstra algorithm is ran from destination nodes on reversed graph

cppRouting v1.1
===============

First CRAN release
