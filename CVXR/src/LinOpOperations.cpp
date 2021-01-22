//   Copyright 2017 Steven Diamond
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "LinOp.h"
#include "LinOpOperations.h"
#include "Utils.h"
#include <cassert>
#include <map>
#include <iostream>

/***********************
 * FUNCTION PROTOTYPES *
 ***********************/
std::vector<Matrix> build_vector(Matrix &mat);
std::vector<Matrix> get_sum_coefficients(LinOp &lin);
std::vector<Matrix> get_sum_entries_mat(LinOp &lin);
std::vector<Matrix> get_trace_mat(LinOp &lin);
std::vector<Matrix> get_neg_mat(LinOp &lin);
std::vector<Matrix> get_div_mat(LinOp &lin);
std::vector<Matrix> get_promote_mat(LinOp &lin);
std::vector<Matrix> get_mul_mat(LinOp &lin);
std::vector<Matrix> get_mul_elemwise_mat(LinOp &lin);
std::vector<Matrix> get_rmul_mat(LinOp &lin);
std::vector<Matrix> get_index_mat(LinOp &lin);
std::vector<Matrix> get_transpose_mat(LinOp &lin);
std::vector<Matrix> get_reshape_mat(LinOp &lin);
std::vector<Matrix> get_diag_vec_mat(LinOp &lin);
std::vector<Matrix> get_diag_matrix_mat(LinOp &lin);
std::vector<Matrix> get_upper_tri_mat(LinOp &lin);
std::vector<Matrix> get_conv_mat(LinOp &lin);
std::vector<Matrix> get_hstack_mat(LinOp &lin);
std::vector<Matrix> get_vstack_mat(LinOp &lin);
std::vector<Matrix> get_kron_mat(LinOp &lin);

/**
 * Computes a vector of coefficient matrices for the linOp LIN based on the
 * type of linOp.
 *
 * Note: This function assumes LIN has been initialized with the correct
 * data, size, and arguments for each linOp type. No error-checking or
 * error-handling for these types of errors is performed.
 *
 * Parameters: LinOp node LIN
 *
 * Returns: std::vector of sparse coefficient matrices for LIN
 */
std::vector<Matrix> get_func_coeffs(LinOp& lin) {
	std::vector<Matrix> coeffs;
	switch (lin.type) {
	case PROMOTE:
		coeffs = get_promote_mat(lin);
		break;
	case MUL_EXPR:
		coeffs = get_mul_mat(lin);
		break;
	case RMUL_EXPR:
		coeffs = get_rmul_mat(lin);
		break;
	case MUL_ELEM:
		coeffs = get_mul_elemwise_mat(lin);
		break;
	case DIV:
		coeffs = get_div_mat(lin);
		break;
	case SUM:
#ifdef _R_DEBUG_
	  Rcpp::Rcout << "SUM" << std::endl;    
#endif  
		coeffs = get_sum_coefficients(lin);
		break;
	case NEG:
#ifdef _R_DEBUG_
	  Rcpp::Rcout << "NEG" << std::endl;    
#endif  
		coeffs = get_neg_mat(lin);
		break;
	case INDEX:
#ifdef _R_DEBUG_
	  Rcpp::Rcout << "INDEX" << std::endl;    
#endif  
		coeffs = get_index_mat(lin);
		break;
	case TRANSPOSE:
		coeffs = get_transpose_mat(lin);
		break;
	case SUM_ENTRIES:
#ifdef _R_DEBUG_
	  Rcpp::Rcout << "SUM_ENTRIES" << std::endl;    
#endif  
		coeffs = get_sum_entries_mat(lin);
		break;
	case TRACE:
		coeffs = get_trace_mat(lin);
		break;
	case RESHAPE:
		coeffs = get_reshape_mat(lin);
		break;
	case DIAG_VEC:
		coeffs = get_diag_vec_mat(lin);
		break;
	case DIAG_MAT:
		coeffs = get_diag_matrix_mat(lin);
		break;
	case UPPER_TRI:
		coeffs = get_upper_tri_mat(lin);
		break;
	case CONV:
		coeffs = get_conv_mat(lin);
		break;
	case HSTACK:
		coeffs = get_hstack_mat(lin);
		break;
	case VSTACK:
		coeffs = get_vstack_mat(lin);
		break;
	case KRON:
		coeffs = get_kron_mat(lin);
		break;
	default:
#ifdef _R_INTERFACE_
	        Rcpp::stop("Error: linOp type invalid");
#else
	        std::cerr << "Error: linOp type invalid." << std::endl;
		exit(-1);
#endif
	}
	return coeffs;
}

/*******************
 * HELPER FUNCTIONS
 *******************/

/**
 * Returns a vector containing the sparse matrix MAT
 */
std::vector<Matrix> build_vector(Matrix &mat) {
	std::vector<Matrix> vec;
	vec.push_back(mat);
	return vec;
}

/**
 * Returns an N x N sparse identity matrix.
 */
Matrix sparse_eye (int n) {
	Matrix eye_n(n, n);
	eye_n.setIdentity();
	return eye_n;
}

/**
 * Returns a sparse ROWS x COLS matrix of all ones.
 *
 * TODO: This function returns a sparse representation of a dense matrix,
 * which might not be extremely efficient, but does make it easier downstream.
 */
Matrix sparse_ones(int rows, int cols)
{
	Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(rows, cols);
	Matrix mat = ones.sparseView();
	return mat;
}

/**
 * Reshapes the input matrix into a single column vector that preserves
 * columnwise ordering. Equivalent to Matlab's (:) operation.
 *
 * Params: sparse Eigen matrix MAT of size ROWS by COLS.
 * Returns: sparse Eigen matrix OUT of size ROWS * COLS by 1
 */

Matrix sparse_reshape_to_vec(Matrix &mat) {
	int rows = mat.rows();
	int cols = mat.cols();
	Matrix out(rows * cols, 1);
	std::vector<Triplet> tripletList;
	tripletList.reserve(rows * cols);
	for ( int k = 0; k < mat.outerSize(); ++k) {
		for (Matrix::InnerIterator it(mat, k); it; ++it) {
			tripletList.push_back(Triplet(it.col() * rows + it.row(), 0,
			                              it.value()));
		}
	}
	out.setFromTriplets(tripletList.begin(), tripletList.end());
	out.makeCompressed();
	return out;
}

/**
 * Builds and returns the stacked coefficient matrices for both horizontal
 * and vertical stacking linOps.
 *
 * Constructs coefficient matrix COEFF for each argument of LIN. COEFF is
 * essentially an identity matrix with an offset to account for stacking.
 *
 * If the stacking is vertical, the columns of the matrix for each argument
 * are interleaved with each other. Otherwise, if the stacking is horizontal,
 * the columns are laid out in the order of the arguments.
 *
 * Parameters: linOP LIN that performs a stacking operation (HSTACK or VSTACK)
 * 						 boolean VERTICAL: True if vertical stack. False otherwise.
 *
 * Returns: vector COEFF_MATS containing the stacked coefficient matrices
 * 					for each argument.
 *
 */
std::vector<Matrix> stack_matrices(LinOp &lin, bool vertical) {
	std::vector<Matrix> coeffs_mats;
	int offset = 0;
	int num_args = lin.args.size();
	for (int idx = 0; idx < num_args; idx++) {
		LinOp arg = *lin.args[idx];

		/* If VERTICAL, columns that are interleaved. Otherwise, they are
			 laid out in order. */
		int column_offset;
		int offset_increment;
		if (vertical) {
			column_offset = lin.size[0];
			offset_increment = arg.size[0];
		} else {
			column_offset = arg.size[0];
			offset_increment = arg.size[0] * arg.size[1];
		}

		std::vector<Triplet> tripletList;
		tripletList.reserve(arg.size[0] * arg.size[1]);
		for (int i = 0; i < arg.size[0]; i++) {
			for (int j = 0; j < arg.size[1]; j++) {
				int row_idx = i + (j * column_offset) + offset;
				int col_idx = i + (j * arg.size[0]);
				tripletList.push_back(Triplet(row_idx, col_idx, 1));
			}
		}

		Matrix coeff(lin.size[0] * lin.size[1], arg.size[0] * arg.size[1]);
		coeff.setFromTriplets(tripletList.begin(), tripletList.end());
		coeff.makeCompressed();
		coeffs_mats.push_back(coeff);
		offset += offset_increment;
	}
	return coeffs_mats;
}

/******************
 * The remaining helper functions are all used to retrieve constant
 * data from the linOp object. Depending on the interface to the calling
 * package, the implementation of these functions will have to change
 * accordingly!
 ******************/


/**
 * Returns the matrix stored in the data field of LIN as a sparse eigen matrix
 * If COLUMN is true, the matrix is reshaped into a column vector which
 * preserves the columnwise ordering of the elements, equivalent to
 * matlab (:) operator.
 *
 * Note all matrices are returned in a sparse representation to force
 * sparse matrix operations in build_matrix.
 *
 * Params: LinOp LIN with DATA containing a 2d vector representation of a
 * 				 matrix. boolean COLUMN
 *
 * Returns: sparse eigen matrix COEFFS
 *
 */
Matrix get_constant_data(LinOp &lin, bool column) {
	Matrix coeffs;
	if (lin.sparse) {
		if (column) {
			coeffs = sparse_reshape_to_vec(lin.sparse_data);
		} else {
			coeffs = lin.sparse_data;
		}
	} else {
		if (column) {
			Eigen::Map<Eigen::MatrixXd> column(lin.dense_data.data(),
			                                   lin.dense_data.rows() *
			                                   lin.dense_data.cols(), 1);
			coeffs = column.sparseView();
		} else {
			coeffs = lin.dense_data.sparseView();
		}
	}
	coeffs.makeCompressed();
	return coeffs;
}

/**
 * Interface for the INDEX linOp to retrieve slice data. Assumes that the
 * INDEX linOp stores slice data in the following format
 *
 * 		vector(row_data, col_data),
 *
 * where row_data = vector(start_idx, end_idx, step_size) and
 * col_data = vector(start_idx, end_idx, step_size).
 *
 * Implements Python slice semantics, e.g. if ROW_START (or ROW_END) is
 * negative, it is assumed to refer to ROWS + ROW_START (ROWS + ROW_END).
 * Same for columns.
 *
 * Parameters: linOp LIN with type INDEX and slice data.
 *
 * Returns: a std::vector containing 2 std::vector of ints.
 * 					The first vector is the row slice data in the form
 * 							(start, end, step_size)
 * 					The second the vector is the column slice data in the form
 * 							(start, end, step_size)
 */
std::vector<std::vector<int> > get_slice_data(LinOp &lin, int rows, int cols) {
	assert(lin.type == INDEX);
	std::vector<int> row_slice = lin.slice[0];
	std::vector<int> col_slice = lin.slice[1];
	assert(row_slice.size() == 3);
	assert(col_slice.size() == 3);
  
	std::vector<std::vector<int> > slices;
	slices.push_back(row_slice);
	slices.push_back(col_slice);
	return slices;
}

/**
 * Interface for the DIV linOp to retrieve the constant divisor.
 *
 * Parameters: linOp LIN of type DIV with a scalar divisor stored in the
 * 							0,0 component of the DENSE_DATA matrix.
 *
 * Returns: scalar divisor
 */
double get_divisor_data(LinOp &lin) {
	assert(lin.type == DIV);
	return lin.dense_data(0, 0);
}

/**
 * Interface for the VARIABLE linOp to retrieve its variable ID.
 *
 * Parameters: linOp LIN of type VARIABLE with a variable ID in the
 * 							0,0 component of the DENSE_DATA matrix.
 *
 * Returns: integer variable ID
 */
int get_id_data(LinOp &lin) {
	assert(lin.type == VARIABLE);
	return int(lin.dense_data(0, 0));
}

/*****************************
 * LinOP -> Matrix FUNCTIONS
 *****************************/
/**
 * Return the coefficients for KRON.
 *
 * Parameters: linOp LIN with type KRON
 * Returns: vector containing the coefficient matrix for the Kronecker
 						product.
 */
std::vector<Matrix> get_kron_mat(LinOp &lin) {
	assert(lin.type == KRON);
	Matrix constant = get_constant_data(lin, false);
	int lh_rows = constant.rows();
	int lh_cols = constant.cols();
	int rh_rows =  lin.args[0]->size[0];
	int rh_cols =  lin.args[0]->size[1];

	int rows = rh_rows * rh_cols * lh_rows * lh_cols;
	int cols = rh_rows * rh_cols;
	Matrix coeffs(rows, cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(rh_rows * rh_cols * constant.nonZeros());
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			int row = (rh_rows * rh_cols * (lh_rows * it.col())) + (it.row() * rh_rows);
			int col = 0;
			for(int j = 0; j < rh_cols; j++){
				for(int i = 0; i < rh_rows; i++) {
					tripletList.push_back(Triplet(row + i, col, it.value()));
					col++;
				}
				row += lh_rows * rh_rows;
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for VSTACK.
 *
 * Parameters: linOp LIN with type VSTACK
 * Returns: vector of coefficient matrices for each argument.
 */
std::vector<Matrix> get_vstack_mat(LinOp &lin) {
	assert(lin.type == VSTACK);
	return stack_matrices(lin, true);
}

/**
 * Return the coefficients for HSTACK
 *
 * Parameters: linOp LIN with type HSTACK
 * Returns: vector of coefficient matrices for each argument.
 */
std::vector<Matrix> get_hstack_mat(LinOp &lin) {
	assert(lin.type == HSTACK);
	return stack_matrices(lin, false);
}

/**
 * Return the coefficients for CONV operator. The coefficient matrix is
 * constructed by creating a toeplitz matrix with the constant vector
 * in DATA as the columns. Multiplication by this matrix is equivalent
 * to convolution.
 *
 * Parameters: linOp LIN with type CONV. Data should should contain a
 *						 column vector that the variables are convolved with.
 *
 * Returns: vector of coefficients for convolution linOp
 */
std::vector<Matrix> get_conv_mat(LinOp &lin) {
	assert(lin.type == CONV);
	Matrix constant = get_constant_data(lin, false);
	int rows = lin.size[0];
	int nonzeros = constant.rows();
	int cols = lin.args[0]->size[0];

	Matrix toeplitz(rows, cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(nonzeros * cols);
	for (int col = 0; col < cols; col++) {
		int row_start = col;
		for ( int k = 0; k < constant.outerSize(); ++k ) {
			for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
				int row_idx = row_start + it.row();
				tripletList.push_back(Triplet(row_idx, col, it.value()));
			}
		}
	}
	toeplitz.setFromTriplets(tripletList.begin(), tripletList.end());
	toeplitz.makeCompressed();
	return build_vector(toeplitz);
}

/**
 * Return the coefficients for UPPER_TRI: an ENTRIES by ROWS * COLS matrix
 * where the i, j entry in the original matrix has a 1 in row COUNT and
 * corresponding column if j > i and 0 otherwise.
 *
 * Parameters: LinOp with type UPPER_TRI.
 * Returns: vector of coefficients for upper triangular matrix linOp
 */
std::vector<Matrix> get_upper_tri_mat(LinOp &lin) {
	assert(lin.type == UPPER_TRI);
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];

	int entries = lin.size[0];
	Matrix coeffs(entries, rows * cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(entries);
	int count = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (j > i) {
				// index in the extracted vector
				int row_idx = count;
				count++;
				// index in the original matrix
				int col_idx = j * rows + i;
				tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for DIAG_MAT (diagonal matrix to vector): a
 * N by N^2 matrix where each row has a 1 in the row * N + row entry
 * corresponding to the diagonal and 0 otherwise.
 *
 * Parameters: LinOp of type DIAG_MAT
 *
 * Returns: vector containing coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_diag_matrix_mat(LinOp &lin) {
	assert(lin.type == DIAG_MAT);
	int rows = lin.size[0];

	Matrix coeffs(rows, rows * rows);
	std::vector<Triplet> tripletList;
	tripletList.reserve(rows);
	for (int i = 0; i < rows; i++) {
		// index in the extracted vector
		int row_idx = i;
		// index in the original matrix
		int col_idx = i * rows + i;
		tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
	}

	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for DIAG_VEC (vector to diagonal matrix): a
 * N^2 by N matrix where each column I has a 1 in row I * N + I
 * corresponding to the diagonal entry and 0 otherwise.
 *
 * Parameters: linOp of type DIAG_VEC
 *
 * Returns: vector containing coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_diag_vec_mat(LinOp &lin) {
	assert(lin.type == DIAG_VEC);
	int rows = lin.size[0];

	Matrix coeffs(rows * rows, rows);
	std::vector<Triplet> tripletList;
	tripletList.reserve(rows);
	for (int i = 0; i < rows; i++) {
		// index in the diagonal matrix
		int row_idx = i * rows + i;
		//index in the original vector
		int col_idx = i;
		tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for TRANSPOSE: a ROWS*COLS by ROWS*COLS matrix
 * such that element ij in the vectorized matrix is mapped to ji after
 * multiplication (i.e. entry (rows * j + i, i * cols + j) = 1 and else 0)
 *
 * Parameters: linOp of type TRANSPOSE
 *
 * Returns: vector containing coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_transpose_mat(LinOp &lin) {
	assert(lin.type == TRANSPOSE);
	int rows = lin.size[0];
	int cols = lin.size[1];

	Matrix coeffs(rows * cols, rows * cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(rows * cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int row_idx = rows * j + i;
			int col_idx = i * cols + j;
			tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for INDEX: a N by ROWS*COLS matrix
 * where N is the number of total elements in the slice. Element i, j
 * is 1 if element j in the vectorized matrix is the i-th element of the
 * slice and 0 otherwise.
 *
 * Parameters: LinOp of type INDEX
 *
 * Returns: vector containing coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_index_mat(LinOp &lin) {
	assert(lin.type == INDEX);
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];
	Matrix coeffs (lin.size[0] * lin.size[1], rows * cols);

	/* If slice is empty, return empty matrix */
	if (coeffs.rows () == 0 ||  coeffs.cols() == 0) {
		return build_vector(coeffs);
	}

#ifdef _R_INTERFACE_	
	
	/* Set the index coefficients by looping over the column selection
	 * first to remain consistent with CVXPY. */
	std::vector<Triplet> tripletList;
	std::vector<int> col_slice = lin.slice[1];
	std::vector<int> row_slice = lin.slice[0];
	int counter = 0;
	for (int j = 0; j < col_slice.size(); j++) {
	  for (int i = 0; i < row_slice.size(); i++) {
#ifdef _R_DEBUG_
	    Rcpp::Rcout << "i, j: "  << col_slice[j] << ", " << row_slice[i] << std::endl;	      
#endif
	    int row_idx = counter;
	    int col_idx = col_slice[j] * rows + row_slice[i];
	    tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
	    counter++;
	  }
	}
#else
	std::vector<std::vector<int> > slices = get_slice_data(lin, rows, cols);

	/* Row Slice Data */
	int row_start = slices[0][0];
	int row_end = slices[0][1];
	int row_step = slices[0][2];

	/* Column Slice Data */
	int col_start = slices[1][0];
	int col_end = slices[1][1];
	int col_step = slices[1][2];

	/* Set the index coefficients by looping over the column selection
	 * first to remain consistent with CVXPY. */
	std::vector<Triplet> tripletList;
	int col = col_start;
	int counter = 0;
	while (true) {
		if (col < 0 || col >= cols) {
			break;
		}
		int row = row_start;
		while (true) {
			if (row < 0 || row >= rows) {
				break;
			}
			int row_idx = counter;
			int col_idx = col * rows + row;
			tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
			counter++;
			row += row_step;
			if ((row_step > 0 && row >= row_end) || (row_step < 0 && row <= row_end)) {
				break;
			}
		}
		col += col_step;
		if ((col_step > 0 && col >= col_end) || (col_step < 0 && col <= col_end)) {
			break;
		}
	}
#endif
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
#ifdef _R_DEBUG_
	Rcpp::Rcout << Eigen::MatrixXd(coeffs) << std::endl;
#endif	
	return build_vector(coeffs);
}

/**
 * Return the coefficients for MUL_ELEM: an N x N diagonal matrix where the
 * n-th element on the diagonal corresponds to the element n = j*rows + i in
 * the data matrix CONSTANT.
 *
 * Parameters: linOp of type MUL_ELEM
 *
 * Returns: vector containing the coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_mul_elemwise_mat(LinOp &lin) {
	assert(lin.type == MUL_ELEM);
	Matrix constant = get_constant_data(lin, true);
	int n = constant.rows();

	// build a giant diagonal matrix
	std::vector<Triplet> tripletList;
	tripletList.reserve(n);
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			tripletList.push_back(Triplet(it.row(), it.row(), it.value()));
		}
	}
	Matrix coeffs(n, n);
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for RMUL_EXPR (right multiplication): a ROWS * N
 * by COLS * N matrix given by the kronecker product between the
 * transpose of the constant matrix CONSTANT and a N x N identity matrix.
 *
 * Parameters: linOp of type RMUL_EXPR
 *
 * Returns: vector containing the corresponding coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_rmul_mat(LinOp &lin) {
	assert(lin.type == RMUL_EXPR);
	Matrix constant = get_constant_data(lin, false);
	int rows = constant.rows();
	int cols = constant.cols();
	int n = lin.size[0];

	Matrix coeffs(cols * n, rows * n);
	std::vector<Triplet> tripletList;
	tripletList.reserve(n * constant.nonZeros());
	for ( int k = 0; k < constant.outerSize(); ++k ) {
		for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
			double val = it.value();

			// each element of CONSTANT occupies an N x N block in the matrix
			int row_start = it.col() * n;
			int col_start = it.row() * n;
			for (int i = 0; i < n; i++) {
				int row_idx = row_start + i;
				int col_idx = col_start + i;
				tripletList.push_back(Triplet(row_idx, col_idx, val));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for MUL_EXPR (left multiplication): a NUM_BLOCKS * ROWS
 * by NUM_BLOCKS * COLS block diagonal matrix where each diagonal block is the
 * constant data BLOCK.
 *
 * Parameters: linOp with type MUL_EXPR
 *
 * Returns: vector containing coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_mul_mat(LinOp &lin) {
	assert(lin.type == MUL_EXPR);
	Matrix block = get_constant_data(lin, false);
	int block_rows = block.rows();
	int block_cols = block.cols();

	// Don't replicate scalars
	if(block_rows == 1 && block_cols == 1){
		return build_vector(block);
	}

	int num_blocks = lin.size[1];
	Matrix coeffs (num_blocks * block_rows, num_blocks * block_cols);

	std::vector<Triplet> tripletList;
	tripletList.reserve(num_blocks * block.nonZeros());
	for (int curr_block = 0; curr_block < num_blocks; curr_block++) {
		int start_i = curr_block * block_rows;
		int start_j = curr_block * block_cols;
		for ( int k = 0; k < block.outerSize(); ++k ) {
			for ( Matrix::InnerIterator it(block, k); it; ++it ) {
				tripletList.push_back(Triplet(start_i + it.row(), start_j + it.col(),
				                              it.value()));
			}
		}
	}
	coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for PROMOTE: a column vector of size N with all
 * entries 1. Note this is treated as sparse for consistency of later
 * multiplications with sparse matrices.
 *
 * Parameters: linOP with type PROMOTE
 *
 * Returns: vector containing coefficient matrix ONES.
 *
 */
std::vector<Matrix> get_promote_mat(LinOp &lin) {
	assert(lin.type == PROMOTE);
	int num_entries = lin.size[0] * lin.size[1];
	Matrix ones = sparse_ones(num_entries, 1);
	ones.makeCompressed();
	return build_vector(ones);
}

/**
 * Return the coefficients for RESHAPE: a 1x1 matrix [1]. In Eigen, this
 * requires special case handling to multiply against an arbitrary m x n
 * matrix.
 *
 * Parameters: LinOp with type RESHAPE
 *
 * Returns: vector containing the coefficient matrix ONE.
 *
 */
std::vector<Matrix> get_reshape_mat(LinOp &lin) {
	assert(lin.type == RESHAPE);
	Matrix one(1, 1);
	one.insert(0, 0) = 1;
	one.makeCompressed();
	return build_vector(one);
}

/**
 * Return the coefficients for DIV: a diagonal matrix where each diagonal
 * entry is 1 / DIVISOR.
 *
 * Parameters: linOp with type DIV
 *
 * Returns: vector containing the coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_div_mat(LinOp &lin) {
	assert(lin.type == DIV);
	// assumes scalar divisor
	double divisor = get_divisor_data(lin);
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs /= divisor;
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for NEG: -I, where I is an identity of size m * n.
 *
 * Parameters: linOp with type NEG
 *
 * Returns: vector containing the coefficient matrix COEFFS
 */
std::vector<Matrix> get_neg_mat(LinOp &lin) {
	assert(lin.type == NEG);
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs *= -1;
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficients for TRACE: A single row vector v^T \in R^(n^2)
 * with 1 if v_{i}  corresponds to a diagonal entry (i.e. i * n + i) and 0
 * otherwise.
 *
 * Parameters: LinOp with type TRACE
 *
 * Returns: vector containing the coefficient matrix COEFFS
 *
 */
std::vector<Matrix> get_trace_mat(LinOp &lin) {
	assert(lin.type == TRACE);
	int rows = lin.args[0]->size[0];
	Matrix coeffs (1, rows * rows);
	for (int i = 0; i < rows; i++) {
		coeffs.insert(0, i * rows + i) = 1;
	}
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficient matrix for SUM_ENTRIES. A single row vector of 1's
 * of size 1 by (data.rows x data.cols).
 *
 * Parameters: LinOp with type SUM_ENTRIES
 *
 * Returns: vector containing the coefficient matrix COEFFS
 */
std::vector<Matrix> get_sum_entries_mat(LinOp &lin) {
	assert(lin.type == SUM_ENTRIES);
	// assumes all args have the same size
	int rows = lin.args[0]->size[0];
	int cols = lin.args[0]->size[1];
	Matrix coeffs = sparse_ones(1, rows * cols);
	coeffs.makeCompressed();
	return build_vector(coeffs);
}

/**
 * Return the coefficient matrix for SUM. Note each element is a 1x1 matrix,
 * i.e. a scalar, and requires special case handling Eigen to multiply with
 * a general m x n matrix.
 *
 * Parameters: LinOp with type SUM
 *
 * Returns: A vector of length N where each element is a 1x1 matrix
 */
std::vector<Matrix> get_sum_coefficients(LinOp &lin) {
	assert(lin.type == SUM);
	int n = lin.args.size();
	std::vector<Matrix> coeffs;
	Matrix scalar(1, 1);
	scalar.insert(0, 0) = 1;
	scalar.makeCompressed();
	for (int i = 0; i < n; i++) {
		coeffs.push_back(scalar);
	}
	return coeffs;
}

/**
 * Return a map from the variable ID to the coefficient matrix for the
 * corresponding VARIABLE linOp, which is an identity matrix of total
 * linop size x total linop size.
 *
 * Parameters: VARIABLE Type LinOp LIN
 *
 * Returns: Map from VARIABLE_ID to coefficient matrix COEFFS for LIN
 *
 */
std::map<int, Matrix> get_variable_coeffs(LinOp &lin) {
	assert(lin.type == VARIABLE);
	std::map<int, Matrix> id_to_coeffs;
	int id = get_id_data(lin);

	// create a giant identity matrix
	int n = lin.size[0] * lin.size[1];
	Matrix coeffs = sparse_eye(n);
	coeffs.makeCompressed();
	id_to_coeffs[id] = coeffs;
	return id_to_coeffs;
}

/**
 * Returns a map from CONSTANT_ID to the data matrix of the corresponding
 * CONSTANT type LinOp. The coefficient matrix is the data matrix reshaped
 * as a ROWS * COLS by 1 column vector.
 * Note the data is treated as sparse regardless of the underlying
 * representation.
 *
 * Parameters: CONSTANT linop LIN
 *
 * Returns: map from CONSTANT_ID to the coefficient matrix COEFFS for LIN.
 */
std::map<int, Matrix> get_const_coeffs(LinOp &lin) {
	assert(lin.has_constant_type());
	std::map<int, Matrix> id_to_coeffs;
	int id = CONSTANT_ID;

	// get coeffs as a column vector
	Matrix coeffs = get_constant_data(lin, true);
	coeffs.makeCompressed();
	id_to_coeffs[id] = coeffs;
	return id_to_coeffs;
}
