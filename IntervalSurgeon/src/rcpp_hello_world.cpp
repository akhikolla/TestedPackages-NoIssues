#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

IntegerVector rcpp_hyper_cube_sizes(IntegerMatrix lengths) {
	int k = lengths.ncol();
	int n = lengths.nrow();
	IntegerVector cube_sizes(n);
	for (int i = 0; i < n; i++) {
		int arr_len = 1;
		for (int j = 0; j < k; j++)
			arr_len *= lengths(i, j);
		cube_sizes[i] = arr_len;
	}
	return cube_sizes;
}

// [[Rcpp::export]]
IntegerMatrix rcpp_hyper_cubes(IntegerVector vals, IntegerVector starts, IntegerMatrix lengths) {
	int k = lengths.ncol();
	int n = lengths.nrow();
	IntegerVector cube_sizes = rcpp_hyper_cube_sizes(lengths);
	int total_tupes = 0;
	for (int i = 0; i < n; i++)
		total_tupes += cube_sizes[i];
	IntegerMatrix result(total_tupes, k+1);
	int place = 0;
	IntegerVector curVals(k, 0);
	for (int i = 0; i < n; i++) {
		if (cube_sizes[i] > 0) {
			for (int label_cube = 0; label_cube < cube_sizes[i]; label_cube++)
				result(place + label_cube, k) = i;
			int give_each = cube_sizes[i];
			int batches = 1;
			for (int j = 0; j < k; j++) {
				give_each /= lengths(i, j);
				batches *= lengths(i, j);
				for (int batch = 0; batch < batches; batch++) {
					for (int row = 0; row < give_each; row++)
						result(place + batch * give_each + row, j) = vals[starts[j] + curVals[j] + batch % lengths(i, j)];
				}
			}
			place += cube_sizes[i];
		}
		for (int j = 0; j < k; j++)
			curVals[j] += lengths(i, j);
	}
	return result;
}

// [[Rcpp::export]]
IntegerVector rcpp_depth(IntegerVector sorted_starts, IntegerVector sorted_ends, IntegerVector pts) {
	int n = pts.length()-1;
	int intervals = sorted_starts.length();
	IntegerVector depth(n);
	int current_depth = 0;
	int current_start = 0;
	int current_end = 0;
	for (int i = 0; i < n; i++) {
		while (current_start < intervals && sorted_starts[current_start] == pts[i]) {
		       current_depth++;
		       current_start++;
		}
		while (current_end < intervals && sorted_ends[current_end] == pts[i]) {
		       current_depth--;
		       current_end++;
		}
		depth[i] = current_depth;
	}
	return depth;
}

// [[Rcpp::export]]
IntegerVector rcpp_pile(IntegerVector starts, IntegerVector ends, IntegerVector pts, int total_members) {
	IntegerVector membership(total_members);
	int n = pts.length()-1;
	int intervals = starts.length();
	int current_depth = 0;
	int current_start = 0;
	int current_member = 0;
	for (int i = 0; i < n; i++) {
		int last_round_members = current_depth;
		int last_round_first_member = current_member - current_depth;
		for (int mem = last_round_first_member; mem < last_round_first_member + last_round_members; mem++) {
			if (ends[membership[mem]] != pts[i]) {
				membership[current_member] = membership[mem];
				current_member++;
			} else {
				current_depth--;
			}
		}
		while (current_start < intervals && starts[current_start] == pts[i]) {
			if (ends[current_start] != pts[i]) {
				membership[current_member] = current_start;
				current_member++;
				current_depth++;
			}
			current_start++;
		}
	}
	return membership;
}

// [[Rcpp::export]]
LogicalVector dash_set_overlaps(IntegerVector starts1, IntegerVector ends1, IntegerVector starts2, IntegerVector ends2, bool state1, bool state2, bool op_is_and, IntegerVector pts) {
	int n = pts.length() - 1;
	int n1 = starts1.length();
	int n2 = starts2.length();
	LogicalVector result(n, false);
	bool open1 = false;
	bool open2 = false;
	int current1 = 0;
	int current2 = 0;
	int i = 0;
	while (i < n)  {
		if (current1 < n1 && pts[i] == starts1[current1]) {
			open1 = true;
		}
		if (current2 < n2 && pts[i] == starts2[current2]) {
			open2 = true;
		}
		if (current1 < n1 && pts[i] == ends1[current1]) {
			open1 = false;
			current1++;
		}
		if (current2 < n2 && pts[i] == ends2[current2]) {
			open2 = false;
			current2++;
		}
		if (op_is_and) {
			if (open1 == state1 && open2 == state2) {
				result[i] = true;
			}
		}
		else {
			if (open1 == state1 || open2 == state2) {
				result[i] = true;
			}
		}
		i++;
	}
	return result;
}
