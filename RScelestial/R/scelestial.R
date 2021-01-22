#' Read mutation table
#' 
#' A simple read of a sequencing file.
#'
#' @param file.name Name of the file to be loaded
#' @return A table representing the content of the file. 
#'   First column of the file represents the row names.
#' 
#' @examples
#' # An example input without header could be like following:
#' # 1 C/C A/A A/A A/A 
#' # 2 ./. A/A C/C C/C 
#' # 3 C/C A/A C/C ./. 
#' # 4 A/A ./. ./. ./. 
#' # 5 ./. A/A A/A A/A 
#' # 
#' # For this file you can run
#' read.sequence.table(system.file("extdata/sample1.txt", package="RScelestial"))
#' 
read.sequence.table <- function(file.name) {
	utils::read.table(file.name, row.names = 1, header = TRUE, stringsAsFactors = TRUE);
}

#' Infer the single-cell phylogenetic tree
#' 
#' Performs the Scelestial algorithm and calculates the phylogenetic tree reconstruction 
#' based on an approximation algorithm for Steiner tree problem.
#' 
#' @param seq The sequence matrix. Rows represent loci and columns represent samples.
#'   Elements of the matrix represent 10-state genome sequencing results, or missing values.
#'   I.e each element is in the format "X/Y" where X and Y are from the set {A, T, C, G}. 
#'   There is a special case "./." that represents the missing value.
#' @param mink The minimum k used in the calculation of k-restricted Steiner trees. 
#'   It is supposed to be 3.
#' @param maxk The maximum k used in the calculation of k-restricted Steiner trees.  
#'   When maxk=3, the approximation algorithm produces an 11/6-approximation result.
#'   Increasing k increases the running time as well as the approximation ratio of the algorithm.
#'   maxk should be not less than mink.
#' @param return.graph If TRUE, the actual graph through igraph library is generated and produced.
#' @param root.assign.method,root root.assign.method is the method for choosing the root. \itemize{
#'   \item "none" for undirected tree,
#'   \item "fix" for a tree with \code{root} as its root. 
#'   \item "balance" to let the root to be chosen to produce the most balanced tree.
#'   }
#' @return Returns a list containing following elements: 
#'   \itemize{
#'     \item \code{tree}: A data frame representing edges of the tree. \code{tree$src} is the source of the edge,
#'        \code{tree$dest} represents the destination of the edge, and \code{tree$len} represents its weight (evolutionary distance).
#'     \item \code{input}: input sequences. 
#'     \item \code{sequence}: inferred or imputed sequences for the tree nodes. If the node is already
#'       in the input, sequence represents its missing value imputation, in the case of presence of missing values, 
#'       and if the node is not an input node, the sequence represents inferred sequence for the tree node.
#'     \item \code{graph}: graph. If the return.graph is TRUE, there is an element G that represents the graph 
#'       from the igraph library.
#'   }
#' 
#' @examples
#' ## simulates tumor evolution
#' S = synthesis(10, 10, 2, seed=7)
#' ## convert to 10-state matrix
#' seq = as.ten.state.matrix(S$seqeunce)
#' ## runs the scelestial to generate 4-restricted Steiner trees. It represents the tree and graph
#' SP = scelestial(seq, mink=3, maxk=4, return.graph = TRUE)
#' SP
#' ## Expected output: 
#' # $input
#' #    node   sequence
#' # 1     0 AAXACAAXXA
#' # 2     1 AXXXAXAAXA
#' # 3     2 AXAXCAXXAX
#' # 4     3 AXCCCAXAAX
#' # 5     4 AXCXAXXCAX
#' # 6     5 XXCAXXXXXX
#' # 7     6 XACXACAAAC
#' # 8     7 AXAXXAXAXA
#' # 9     8 AXAAXXAXXX
#' # 10    9 AAXXXXCXCX
#' #
#' # $sequence
#' #    node   sequence
#' # 1     0 AAAACAAACA
#' # 2     1 AACAAAAAAA
#' # 3     2 AAAACAAAAA
#' # 4     3 AACCCAAAAA
#' # 5     4 AACAACACAC
#' # 6     5 AACAACAAAC
#' # 7     6 AACAACAAAC
#' # 8     7 AAAACAAACA
#' # 9     8 AAAACAAACA
#' # 10    9 AAAACACACA
#' # 11   10 AAAACAAACA
#' # 12   16 AACAAAAAAA
#' # 13   18 AACACAAAAA
#' #
#' # $tree
#' #    src dest     len
#' # 1    9   10 4.00006
#' # 2    8   10 3.00006
#' # 3    7   10 2.50005
#' # 4    0   10 1.50003
#' # 5    6   16 3.00002
#' # 6    1   16 2.50005
#' # 7    3   18 2.50003
#' # 8    0   18 1.50003
#' # 9   16   18 1.00000
#' # 10   0    2 3.50008
#' # 11   4    6 4.00007
#' # 12   5    6 4.50010
#' #
#' # $graph
#' # IGRAPH 6ba60f3 DNW- 13 12 --
#' # + attr: name (v/c), weight (e/n)
#' # + edges from 6ba60f3 (vertex names):
#' #  [1] 9 ->10 8 ->10 7 ->10 0 ->10 6 ->16 1 ->16 3 ->18 0 ->18 16->18 0 ->2
#' # [11] 4 ->6  5 ->6
#' #
scelestial <- function(seq, mink=3, maxk=3, root.assign.method = c("none", "balance", "fix"), root = NULL, 
					return.graph = FALSE) {
	if (!is.data.frame(seq)) {
		stop("Input should be a dataframe.")
	}
	nuc = c("A", "C", "T", "G");
	states = c(as.vector(outer(nuc, nuc, function(a,b) {paste(a,b,sep="/")})), "./.");
	if (!all(apply(seq, 1:2, function(e, states) {e %in% states}, states))) {
		stop("Elements of data frame should be of the format 'X/Y' where X, Y are nucleic acids ('A', 'C', 'G', 'T'), or './.'.")
	}
	root.assign.method = match.arg(root.assign.method);
	if (root.assign.method == "fix") {
		if (is.null(root))
			stop("When you use 'fix' root.assign.method, you should provide the root as root argument.");
		# if (nrow(seq) != length(root))
		# 	stop("length of root should be equal to the number of columns of sequence matrix");
		# if (!all(sapply(root, function(e, states) {e %in% states}, states)))
		# 	stop("elements of root should be of the form 'X/Y' where X, Y are nucleic acids ('A', 'C', 'G', 'T') or './.'.")
		if (!root %in% colnames(seq))
			stop("root should be name of a column of seq");
	}
	
	Y = .scelestial(seq, mink, maxk);
	tree = data.frame("src" = Y$Esrc, "dest" = Y$Edst, "len" = Y$Ew, stringsAsFactors = TRUE);
	# S = data.frame("node" = Y$sequence$node, "sequence" = Y$sequence$seq, stringsAsFactors = TRUE)
	S = as.ten.state.matrix.from.node.seq(data.frame("node" = Y$sequence$node, "sequence" = Y$sequence$seq, stringsAsFactors = TRUE))
	# I = data.frame("node" = Y$input$node, "sequence" = Y$input$seq, stringsAsFactors = TRUE)
	I = as.ten.state.matrix.from.node.seq(data.frame("node" = Y$input$node, "sequence" = Y$input$seq, stringsAsFactors = TRUE))
	R = list("input" = I, "sequence" = S, "tree" = tree)
	
	if (return.graph) {
		# library(igraph);
		# dir = TRUE;
		# if (root.assign.method == "none")
		# 	dir = FALSE;
		
		if (root.assign.method == "none") {
			G = igraph::graph(edges = as.character(as.vector(t(tree[c("src", "dest")]))), directed = FALSE);
		} else {
			G.dfs = my.dfs(tree, root = root);
			if (root.assign.method == "balance") {
				root = names(which.min(G.dfs$balance.depth));
				if (length(root) != 1) {
					stop("There are more than one vertex with minimum depth")
				}
				G.dfs = my.dfs(tree, root = root);
			}
			E.new = do.call(rbind, apply(tree, 1, function(r, G.dfs) {
						  	# print(r);
						  	# print(r["src"] %in% names(G.dfs$father));
						  	if (!r["src"] %in% names(G.dfs$father) || G.dfs$father[[r["src"]]] != r["dest"])
						  		data.frame(src=r["src"], dest=r["dest"], len=r["len"], stringsAsFactors = TRUE)
							else
						  		data.frame(src=r["dest"], dest=r["src"], len=r["len"], stringsAsFactors = TRUE)
						  }, G.dfs));
			G = igraph::graph(edges = as.character(as.vector(t(E.new[c("src", "dest")]))), directed = TRUE);
		}
		
		igraph::E(G)$weight = tree$len;
		R[["graph"]] = G;
	}
	return(R);
}

#' Synthesize single-cell data through tumor simulation
#' 
#' This function simulates a evolution in a tumor through two phases: 1) simulation of evolution,
#' 2) sampling. 
#' 
#' The simulation of evolution starts with a single cell. 
#' Then for \code{evolution.step} steps, on each step a cell is selected for duplication. 
#' A new cell as its child is added to 
#' the evolutionary tree. To each node in the evolutionary tree an advantage is assigned
#' representing its relative advantage in replication and in being sampled. Advantage of a node
#' is calculated by increasing (decreasing) its parents advantage by \code{advantage.increase.step} 
#' (\code{advantage.decrease.step}) with probability proportional to \code{advantage.increase.ratio} 
#' (\code{advantage.decrease.ratio}). 
#' With a probability proportional to \code{advantage.keep.ratio} the advantage of a node
#' is equal to its parent's advantage.
#' 
#' Sequences for each node is build based on its parent's sequence by adding some mutations.
#' Mutations are added for each locus independently with rate \code{mutation.rate}.
#' 
#' In the sampling phase, \code{sample} cells are selected from the evolutionary tree nodes.
#' Result of the sequencing process for a cell is determined by the sequence of the node in the evolutionary tree
#' with addition of some random errors. Errors are result of applying some false positives with rate \code{fp.rate},
#' applying some false negatives with rate \code{fn.rate}, and adding some missing values
#' with rate \code{mv.rate}.
#' 
#' @param sample Number of samples.
#' @param site number of sites (loci)
#' @param evolution.step Number of evolutionary steps in the process of production of 
#'   the evolutionary tree.
#' @param mutation.rate The rate of mutation on each evolutionary step in evolutionary tree synthesis.
#' @param advantage.increase.ratio,advantage.decrease.ratio,advantage.keep.ratio A child node
#'   in the evolutionary tree is chosen for increase/decrease/keep its parent advantage with
#'   probabilities proportional to \code{advantage.increase.ratio}/\code{advantage.decrease.ratio}/\code{advantage.keep.ratio}.
#' @param advantage.increase.step,advantage.decrease.step The amount of 
#'   increasing or decreasing the advantage of a cell relative to its parent.
#' @param mv.rate Rate of missing value to be added to the resulting sequences.
#' @param fp.rate,fn.rate Rate of false positive (0 -> 1) and false negative (1 -> 0)
#'   in the sequences.
#' @param seed The seed for randomization.
#' @return The function returns a list. The list consists of 
#'   \itemize{
#'     \item \code{sequence}: A data frame representing
#'       result of sequencing. The data frame has a row for each locus and a column for each sample.
#'     \item \code{true.sequence}: The actual sequence for the sample before adding errors and missing values.
#'     \item \code{true.clone}: A list that stores index of sampled cells for each node in the evolutionary tree.
#'     \item \code{true.tree}: The evolutionary tree that the samples are sampled from. It is a data frame
#'       with \code{src}, \code{dest}, and \code{len} columns representing source, destination and weight of edges of the tree,
#'       respectively.
#'  }
#' 
#' @examples 
#' ## generating a data set with 10 samples and 5 loci through simulation of
#' ## 20-step evolution.
#' synthesis(10, 5, 20, seed=7)
#' ## The result is
#' # $seqeunce
#' #     C1 C2 C3 C4 C5
#' # L1   1  1  1  1  1
#' # L2   3  1  3  3  0
#' # L3   3  1  3  3  1
#' # L4   3  0  1  0  0
#' # L5   1  3  0  3  3
#' # L6   3  1  3  1  0
#' # L7   3  3  1  0  3
#' # L8   3  1  1  3  3
#' # L9   3  3  1  3  1
#' # L10  0  3  0  3  0
#' #
#' # $true.sequence
#' #     C1 C2 C3 C4 C5
#' # L1   0  1  1  1  1
#' # L2   0  1  0  0  1
#' # L3   0  1  0  0  1
#' # L4   0  1  1  1  1
#' # L5   1  1  0  1  0
#' # L6   0  1  0  1  0
#' # L7   0  1  0  0  1
#' # L8   0  1  1  1  1
#' # L9   0  1  1  1  1
#' # L10  0  0  0  0  0
#' #
#' # $true.clone
#' # $true.clone[[1]]
#' # [1] 4
#' #
#' # $true.clone[[2]]
#' # [1] 1
#' #
#' # $true.clone[[3]]
#' # [1] 6
#' #
#' # $true.clone[[4]]
#' # [1] 10
#' #
#' # $true.clone[[5]]
#' # [1] 2
#' #
#' # $true.clone[[6]]
#' # [1] 3
#' #
#' # $true.clone[[7]]
#' # [1] 8 9
#' #
#' # $true.clone[[8]]
#' # [1] 7
#' #
#' # $true.clone[[9]]
#' # [1] 5
#' #
#' #
#' # $true.tree
#' #   src dest len
#' # 1   1    5   3
#' # 2   5    7   1
#' # 3   5   10   2
#' # 4   1   11   3
#' # 5   1   12   2
#' # 6   1   13   3
#' # 7   7   14   2
#' # 8  12   19   1
#' # 9  10   20   1
#' #
synthesis <- function(sample, site, evolution.step, 
					  mutation.rate = 1, 
					  advantage.increase.ratio = 1, advantage.decrease.ratio = 10, advantage.keep.ratio = 100, 
					  advantage.increase.step = 0.01, advantage.decrease.step = 0.01, 
					  mv.rate = 0.5, fp.rate = 0.2, fn.rate = 0.1,
					  seed = -1) {
	X <- .synthesis(sample, site, evolution.step, 
					mutation.rate, 
					advantage.increase.ratio, advantage.decrease.ratio, advantage.keep.ratio, 
					advantage.increase.step, advantage.decrease.step, 
					mv.rate, 
					fp.rate, fn.rate,
					seed);
	
	true.clone = X$true.clone;
	sequence = as.data.frame(X$sequence, stringsAsFactors = TRUE);
	true.sequence = as.data.frame(X$true.sequence, stringsAsFactors = TRUE);
	row.names(sequence) = row.names(true.sequence) =X$locus.names;
	true.tree = data.frame("src" = X$true.tree$d, "dest" = X$true.tree$s, "len" = X$true.tree$w, stringsAsFactors = TRUE);
	Y = list(
		"seqeunce" = sequence,
		"true.sequence" = true.sequence,
		"true.clone" = true.clone,
		"true.tree" = true.tree
	);
	# Y = X
	Y
}

#' Plotting the tree 
#' 
#' Plotting the igraph tree created by scelestial.
#' 
#' @param graph Output of scelestial or the G element of the scelestial output.
#' @param ... Parameters passing to the plot function
#' 
tree.plot <- function(graph, ...) {
	G <- NULL;
	if (igraph::is.igraph(graph)) {
		G <- graph;
	} else if(!is.null(graph[["graph"]]) && igraph::is.igraph(graph[["graph"]])) {
		G <- graph[["graph"]];
	} else {
		stop("Unknown input format. Either pass ourput of the scelestial function or an igraph.")
	}
	#V(G)$size = 50;
	# plot.new()
	# opar <- par()$mar;
	# par(mar=rep(0, 4)) #Give the graph lots of room
	# plot.new()
	# igraph::V(G)$label = V(G)$name;
	l = igraph::layout_as_tree(G)
	igraph::plot.igraph(G, layout = l, 
				# xlim = c(-1, 1),
				# rescale=FALSE, xlim=range(l[,1]), ylim=range(l[,2]),			
				# asp = 0,
				...
	);
	# par(mar=opar)
}

#' Conversion of ten-state sequencing matrix to 0/1-mutation matrix.
#' 
#' @param seq A dataframe representing the ten-state sequencing matrix. Elements of the matrix
#'   are the from "X/Y" for X and Y being nucleotides or "./." for missing value.
#'   Rows represent loci and columns represent samples.
#'   
#' @return A data frame with exactly the same size as the input \code{seq} matrix.
#'   The most abundant state in each loci (row) translated to 0, and 
#'   the others are translated to 1. Missing values are translated to 3.
#' 
#' @examples 
#' ## A small 10-state matrix
#' seq = data.frame("C1" = c("C/C", "C/C"), "C2" = c("A/A", NA), 
#'         "C3" = c("C/C", "A/A"), stringsAsFactors = TRUE)
#' ## Convert it to mutation matrix
#' as.mutation.matrix(seq)
#' #   C1 C2 C3
#' # 1  0  1  0
#' # 2  1  3  0
#' 
as.mutation.matrix <- function(seq) {
	if (!is.data.frame(seq)) {
		stop("Input should be a matrix with loci as rows and cells as columns.")
	}
	X.mutation.matrix = data.frame(t(apply(t(seq), 2, function(R) {
		R.new = factor(rep(1, length(R)), levels=c(0,1,3))
		R[R == "./."] = NA
		R.fact <- factor(R);
		R.tab <- tabulate(R.fact);
		max.level = levels(R.fact)[which.max(R.tab)];
		R.new[R == max.level] = 0
		R.new[is.na(R)] = 3
		R.new
	})), stringsAsFactors = TRUE);
	rownames(X.mutation.matrix) = rownames(seq)
	colnames(X.mutation.matrix) = colnames(seq)
	X.mutation.matrix
}

#' Conversion of 0/1 matrix to 10-state matrix
#' 
#' It converts 0 to A/A and 1 to C/C. 3 that represents missing values are converted
#' to "./.".
#' 
#' @param mut A dataframe representing the mutation matrix.
#' 
#' @note 
#' Note that following function does not provide inverse of as.mutation.matrix.
#' It could be used to generate input for scelestial.
#' 
#' @return A data frame with the exact size as \code{mut}, in which
#'   0, 1 and 3 (or NAs) are replaced with "A/A", "C/C", and "./.",
#'   respectively.
#' 
#' @examples 
#' ## A small 0/1/NA mutation matrix
#' mut = data.frame("C1" = c(0, 0), "C2" = c(0, 3), "C3" = c(1, 0), 
#'         stringsAsFactors = TRUE)
#' ## Convert it to 10-state matrix
#' as.ten.state.matrix(mut)
#' #    C1  C2  C3
#' # 1 A/A A/A C/C
#' # 2 A/A ./. A/A
as.ten.state.matrix <- function(mut) {
	X.new = as.matrix(mut);
	X.new[mut == 0] = "A/A";
	X.new[mut == 1] = "C/C";
	X.new[is.na(mut) | mut == 3] = "./.";
	as.data.frame(X.new, stringsAsFactors = TRUE)
}

#' Generates 10-state sequence matrix from name/10-char string
#' matrix. 
#' 
#' This function is used for conversion of results of internal
#' scelestial result to 10-state sequence matrices.
#' 
#' @param n.seq A two column data frame. First column is the name of a node
#'   and the second column is a string representation of the sequencing result.
#'   Each element of the sequencing result is from a 10-state representation in which
#'   each state represented as a character according to the following encoding:
#'    \tabular{cc}{
#'      One character representation \tab 10-state representation  \cr
#'      "A" \tab "A/A",	\cr
#'      "T" \tab "T/T", \cr
#'      "C" \tab "C/C", \cr
#'      "G" \tab "G/G", \cr
#'      "K" \tab "A/C", \cr
#'      "L" \tab "A/G", \cr
#'      "M" \tab "C/T", \cr
#'      "N" \tab "C/G", \cr
#'      "O" \tab "T/G", \cr
#'      "P" \tab "T/A", \cr
#'      "X" \tab "./."
#'    }
#' 
#' @return A 10-state sequence data frame with samples as columns
#'   and loci as rows. Elements of \code{n.seq} are translated
#'   to their 10-state representations.
#' 
#' @examples 
#' ## A node sequence data framce
#' n.seq = data.frame("node" = c("C1", "C2"), "seq" = c("AKLTCXAAC", "AKKOCXAPC"), 
#'           stringsAsFactors = TRUE)
#' ## Convert it to ten state matrix
#' as.ten.state.matrix.from.node.seq(n.seq)
#' #     V1  V2  V3  V4  V5  V6  V7  V8  V9
#' # C1 A/A A/C A/G T/T C/C ./. A/A A/A C/C
#' # C2 A/A A/C A/C T/G C/C ./. A/A T/A C/C
#' 
as.ten.state.matrix.from.node.seq <- function(n.seq) {
	# print(n.seq);
	X = data.frame(matrix(0, nrow=nrow(n.seq), ncol=0));
	if (nrow(n.seq) == 0 || n.seq[,2][1] == "")
		X
	else {
		nuc = c("A", "C", "T", "G");
		states = c(as.vector(outer(nuc, nuc, function(a,b) {paste(a,b,sep="/")})), "./.");
		
		char.to.ten.state = unlist(list(
			"A" = "A/A",		
			"T" = "T/T",
			"C" = "C/C",
			"G" = "G/G",
			"K" = "A/C",
			"L" = "A/G",
			"M" = "C/T",
			"N" = "C/G",
			"O" = "T/G",
			"P" = "T/A",
			"X" = "./."
		));
		ten.state.to.char = sapply(states, function(s, char.to.ten.state) {
			if (s %in% char.to.ten.state)
				names(which(char.to.ten.state == s))
			else
				names(which(char.to.ten.state == paste(rev(strsplit(s, NULL)[[1]]), collapse='')))
		}, char.to.ten.state, USE.NAMES = TRUE, simplify = TRUE);

		for (i in 1:nchar(as.character(n.seq[,2][1]))) {
			X[[paste("V", i, sep = '')]] <-
			apply(n.seq, 1, function(r, i, char.to.ten.state) char.to.ten.state[substr(r[2], i, i)], i, char.to.ten.state)	
		}
		rownames(X) <- n.seq[,1]
		X
	}
}

#' Running depth first search on a tree and calling
#' functions on entrance/exit events
#' 
#' It is used for internal purposes.
#' 
#' @param nei Neighbor list for each vertex
#' @param v Starting node
#' @param f Parent node
#' @param extra the shared object for the whole DFS
#' @param in.call First function to call
#' @param mid.call.before Function to call before calling child DFS
#' @param mid.call.after Function to call after calling child DFS
#' @param out.call Last function to call
#' 
#' @return the \code{extra} parameter modified with \code{in.call}, \code{mid.call.before}, \code{mid.call.after}, and \code{out.call} functions
#' 
my.general.dfs <- function(nei, v, f, extra, in.call, mid.call.before, mid.call.after, out.call) {
	extra = in.call(nei, v, f, extra);
	for (u in nei[[v]]) {
		if (u != f) {
			extra = mid.call.before(nei, u, v, f, extra);
			extra = my.general.dfs(nei, u, v, extra, in.call, mid.call.before, mid.call.after, out.call);
			extra = mid.call.after(nei, u, v, f, extra);
		}
	}
	extra = out.call(nei, v, f, extra);
	extra
}

#' Runs DFS on tree and calculates parent of each node
#' as well as depth and upper-depth of nodes.
#' 
#' It is used for internal purposes.
#' 
#' @param graph The tree
#' @param root The starting node of DFS.
#' @return a list with \code{father} representing the parent node, and
#'   \code{balance.depth} representing the distance between the node and
#'   the farthest node to it, as the elements.
#' 
my.dfs <- function(graph, root = NULL) {
	if (nrow(graph) == 0) {
		stop("graph is empty");
	}
	vertices <- levels(unlist(graph[1:2]));
	nei <- sapply(vertices, function(x) {
		# print(x);
		c(as.character(graph[graph["src"] == x, "dest"]), as.character(graph[graph["dest"] == x, "src"]))
	}, USE.NAMES = TRUE, simplify = FALSE);
	if (is.null(root)) {
		if (length(vertices) < 1) {
			stop("Graph has no vertex")
		}
		root = vertices[1];
	}
	
	in.call.df <- function(nei, v, f, extra) {
		# print(paste(v, f));
		if (f != -1) {
			extra[["father"]][[v]] = f;
		}
		extra$depth[[v]] = 0;
		extra
	}
	
	mid.call.df <- function(nei, u, v, f, extra) {
		extra$depth[[v]] = max(extra$depth[[v]], extra$depth[[u]]+1);
		extra
	}
	
	ret = my.general.dfs(nei, root, -1, list(father=list(), depth=list(), udepth=list()), 
						 in.call = in.call.df, 
						 mid.call.before = function(nei, u, v, f, extra) extra,
						 mid.call.after = mid.call.df,
						 out.call = function(nei, v, f, extra) extra)
	
	
	mid.call.ud <- function(nei, u, v, f, extra) {
		ud = sapply(nei[[v]], function(x, u, v, f, extra) {
				if (x == -1 || x == f || x == u) 0 else extra$depth[[x]]+2
			}, u, v, f, extra)
		extra$udepth[[u]] = max(1+extra$udepth[[v]], ud);
		
		# print(paste(" m.ud ", "u", u, "v", v, "f", f, "u-ud:", extra$udepth[[u]], "fud", extra$udepth[[v]]))
		extra
	}

	ret$udepth[[root]] = 0
	ret = my.general.dfs(nei, root, -1, ret, 
						 in.call = function(nei, v, f, extra) extra, 
						 mid.call.before = mid.call.ud,
						 mid.call.after = function(nei, u, v, f, extra) extra,
						 out.call = function (nei, v, f, extra) extra)
	list(father=unlist(ret$father), 
		 balance.depth = apply(rbind(ret$udepth[vertices], ret$depth[vertices]), 2, function(a) max(unlist(a)) )
		 )
}

#' Calculates distance matrix for a synthetized data
#' 
#' @param D Output of synthesis function
#' @param normalize If true, sum of all elements of resulting
#'   table is added up to one.
#' @return The distance matrix of the true tree.
#' 
#' @examples 
#' ## Synthesise an evolution
#' S = synthesis(10, 5, 20, seed=7)
#' ## Calculating the distance matrix of the true tree.
#' distance.matrix.true.tree(S)
#' #              C3          C6          C4          C2          C7
#' # C3  0.000000000 0.004587156 0.006880734 0.009174312 0.013761468
#' # C6  0.004587156 0.000000000 0.002293578 0.009174312 0.013761468
#' # C4  0.006880734 0.002293578 0.000000000 0.011467890 0.016055046
#' # C2  0.009174312 0.009174312 0.011467890 0.000000000 0.004587156
#' # C7  0.013761468 0.013761468 0.016055046 0.004587156 0.000000000
#' # C10 0.006880734 0.006880734 0.009174312 0.011467890 0.016055046
#' # C8  0.006880734 0.011467890 0.013761468 0.016055046 0.020642202
#' # C9  0.006880734 0.011467890 0.013761468 0.016055046 0.020642202
#' # C1  0.011467890 0.011467890 0.013761468 0.002293578 0.006880734
#' # C5  0.011467890 0.011467890 0.013761468 0.002293578 0.006880734
#' # C10          C8          C9          C1          C5
#' # C3  0.006880734 0.006880734 0.006880734 0.011467890 0.011467890
#' # C6  0.006880734 0.011467890 0.011467890 0.011467890 0.011467890
#' # C4  0.009174312 0.013761468 0.013761468 0.013761468 0.013761468
#' # C2  0.011467890 0.016055046 0.016055046 0.002293578 0.002293578
#' # C7  0.016055046 0.020642202 0.020642202 0.006880734 0.006880734
#' # C10 0.000000000 0.013761468 0.013761468 0.013761468 0.013761468
#' # C8  0.013761468 0.000000000 0.000000000 0.018348624 0.018348624
#' # C9  0.013761468 0.000000000 0.000000000 0.018348624 0.018348624
#' # C1  0.013761468 0.018348624 0.018348624 0.000000000 0.000000000
#' # C5  0.013761468 0.018348624 0.018348624 0.000000000 0.000000000
#' 
distance.matrix.true.tree <- function(D, normalize = TRUE) {
	vertices = unlist(D$true.clone);
	cellClone = list();
	for (i in 1:length(D$true.clone)) {
		for (l in D$true.clone[[i]]) {
			cellClone[[l]] = names(D$true.clone)[i]
		}
	}
	cellClone = unlist(cellClone);
	distance.matrix.tree(D$true.tree, vertices, cellClone, normalize)
}

#' Calculates distance matrix for result of scelestial
#' 
#' @param SP Output of scelestial function
#' @param normalize If true, sum of all elements of resulting
#'   table is added up to one.
#' @return The distance matrix
#' 
#' @examples 
#' ## Synthesise an evolution
#' S = synthesis(10, 5, 20, seed=7)
#' ## Run Scelestial
#' SC = scelestial(as.ten.state.matrix(S$seqeunce))
#' ## Calculate the distance matrix
#' distance.matrix.scelestial(SC)
#' #              C1         C10          C2          C3          C4
#' # C1  0.000000000 0.003512891 0.015222451 0.014051472 0.008196692
#' # C10 0.003512891 0.000000000 0.011709560 0.010538580 0.004683800
#' # C2  0.015222451 0.011709560 0.000000000 0.010538627 0.007025759
#' # C3  0.014051472 0.010538580 0.010538627 0.000000000 0.005854780
#' # C4  0.008196692 0.004683800 0.007025759 0.005854780 0.000000000
#' # C5  0.011709560 0.008196668 0.003512891 0.007025736 0.003512868
#' # C6  0.023419213 0.019906322 0.019906368 0.009367741 0.015222521
#' # C7  0.018735342 0.015222451 0.015222498 0.004683871 0.010538651
#' # C8  0.015222474 0.011709583 0.014051542 0.012880562 0.007025783
#' # C9  0.010538627 0.007025736 0.009367695 0.008196715 0.002341935
#' # C5          C6          C7          C8          C9
#' # C1  0.011709560 0.023419213 0.018735342 0.015222474 0.010538627
#' # C10 0.008196668 0.019906322 0.015222451 0.011709583 0.007025736
#' # C2  0.003512891 0.019906368 0.015222498 0.014051542 0.009367695
#' # C3  0.007025736 0.009367741 0.004683871 0.012880562 0.008196715
#' # C4  0.003512868 0.015222521 0.010538651 0.007025783 0.002341935
#' # C5  0.000000000 0.016393477 0.011709606 0.010538651 0.005854803
#' # C6  0.016393477 0.000000000 0.004683871 0.022248304 0.017564457
#' # C7  0.011709606 0.004683871 0.000000000 0.017564433 0.012880586
#' # C8  0.010538651 0.022248304 0.017564433 0.000000000 0.004683847
#' # C9  0.005854803 0.017564457 0.012880586 0.004683847 0.000000000
#' 
distance.matrix.scelestial <- function(SP, normalize = TRUE) {
	vertices <- rownames(SP$input);
	distance.matrix.tree(SP$tree, vertices, vertices, normalize)
}

#' Calculates distance matrix for a nodes on a tree.
#' 
#' It is used for internal purposes.
#' 
#' @param graph The tree
#' @param cell.names Name of the cells to be the row and column
#'   name of the resulting matrix
#' @param tree.nodes For each cell.names a tree node is stored
#'   in tree.nodes.
#' @param normalize If TRUE the resulting matrix is normalized.
#' 
#' @return A matrix with equal number of rows and columns, a row/column
#'   for each cell. Elements of matrix represent distance between cells
#'   on the \code{graph}.
#' 
#' @examples 
#' ## Synthesise an evolution
#' S = synthesis(10, 5, 20, seed=7)
#' ## Run Scelestial
#' SC = scelestial(as.ten.state.matrix(S$seqeunce))
#' ## Calculate the distance matrix
#' vertices <- rownames(SC$input);
#' distance.matrix.tree(SC$tree, vertices, vertices, normalize = TRUE)
#' #              C1         C10          C2          C3          C4
#' # C1  0.000000000 0.003512891 0.015222451 0.014051472 0.008196692
#' # C10 0.003512891 0.000000000 0.011709560 0.010538580 0.004683800
#' # C2  0.015222451 0.011709560 0.000000000 0.010538627 0.007025759
#' # C3  0.014051472 0.010538580 0.010538627 0.000000000 0.005854780
#' # C4  0.008196692 0.004683800 0.007025759 0.005854780 0.000000000
#' # C5  0.011709560 0.008196668 0.003512891 0.007025736 0.003512868
#' # C6  0.023419213 0.019906322 0.019906368 0.009367741 0.015222521
#' # C7  0.018735342 0.015222451 0.015222498 0.004683871 0.010538651
#' # C8  0.015222474 0.011709583 0.014051542 0.012880562 0.007025783
#' # C9  0.010538627 0.007025736 0.009367695 0.008196715 0.002341935
#' # C5          C6          C7          C8          C9
#' # C1  0.011709560 0.023419213 0.018735342 0.015222474 0.010538627
#' # C10 0.008196668 0.019906322 0.015222451 0.011709583 0.007025736
#' # C2  0.003512891 0.019906368 0.015222498 0.014051542 0.009367695
#' # C3  0.007025736 0.009367741 0.004683871 0.012880562 0.008196715
#' # C4  0.003512868 0.015222521 0.010538651 0.007025783 0.002341935
#' # C5  0.000000000 0.016393477 0.011709606 0.010538651 0.005854803
#' # C6  0.016393477 0.000000000 0.004683871 0.022248304 0.017564457
#' # C7  0.011709606 0.004683871 0.000000000 0.017564433 0.012880586
#' # C8  0.010538651 0.022248304 0.017564433 0.000000000 0.004683847
#' # C9  0.005854803 0.017564457 0.012880586 0.004683847 0.000000000
#' 
distance.matrix.tree <- function(graph, cell.names, tree.nodes, normalize = TRUE) {
	if (nrow(graph) == 0) {
		stop("graph is empty");
	}
	vertices <- levels(as.factor(unlist(graph[1:2])));

	nei <- sapply(vertices, function(x) {
		# print(x);
		c(as.character(graph[graph["src"] == x, "dest"]), as.character(graph[graph["dest"] == x, "src"]))
	}, USE.NAMES = TRUE, simplify = FALSE);
	
	.edgeslen <- unlist(apply(graph, 1, function(r) {
		# print(r);
		n1 = paste(r["src"], r["dest"], sep=":");
		n2 = paste(r["dest"], r["src"], sep=":");
		l = list(r["len"], r["len"])
		names(l) = c(n1, n2)
		l
	}));
	edgeslen = as.numeric(.edgeslen)
	names(edgeslen) = names(.edgeslen)
	
	X = matrix(0, nrow=length(vertices), ncol=0)
	rownames(X) = vertices
	for (root in vertices) {
		dist = list(0);
		names(dist) = c(root);
		ret = my.general.dfs(nei, root, -1, list(dist = dist,
												 len=edgeslen), 
							 in.call = function(nei, v, f, extra) extra, 
							 mid.call.before = function(nei, u, v, f, extra) {
							 	extra$dist[[u]] = extra$dist[[v]] + extra$len[[paste(paste(v,u,sep=":"),"len", sep=".")]];
							 	extra
							 },
							 mid.call.after = function(nei, u, v, f, extra) extra,
							 out.call = function(nei, v, f, extra) extra)
		X.colnames.prev = colnames(X)
		X <- cbind(X, unlist(ret$dist[vertices]))
		colnames(X) = c(X.colnames.prev, root)
	}
	
	R = cbind(cell.names, tree.nodes);
	O = outer(1:nrow(R), 1:nrow(R), function(a,b, R, X) {
		n1=R[a, "tree.nodes"]; 
		n2=R[b, "tree.nodes"]; 
		# print(cbind(n1,n2));
		X[cbind(n1, n2)]
	}, R, X);
	rownames(O) = colnames(O) = cell.names;
	if (normalize)
		O = O / sum(O);
	O
}


