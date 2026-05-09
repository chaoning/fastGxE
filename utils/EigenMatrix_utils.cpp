#include <cstdint>


/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-06-22 21:08:46
 * LastEditors: Chao Ning
 * LastEditTime: 2026-03-31 15:13:45
 */


#include <cassert>
#include <set>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "EigenMatrix_utils.hpp"
#include "mkl.h"


using namespace Eigen;
using namespace std;



/**
 * @brief Compute the Pearson correlation matrix between columns of a dense matrix.
 *
 * Each column of @p mat is treated as one variable and each row as one sample.
 * The implementation centers the matrix once, reuses the centered values to
 * compute both column standard deviations and the covariance matrix, and then
 * rescales the covariance matrix into a correlation matrix.
 *
 * Zero-variance columns do not have defined off-diagonal correlations. For
 * those columns this function keeps the corresponding row/column entries at 0
 * and forces the diagonal entry to 1.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp` to compute the correlation matrix of `m_bye_mat`
 * - `test/eigenmatrix_utils_example.cpp` as a small runnable example
 */
Eigen::MatrixXd computeCorrelationMatrix(const Eigen::MatrixXd& mat) {
        if (mat.rows() == 0 || mat.cols() == 0) {
            return Eigen::MatrixXd(mat.cols(), mat.cols());
        }

        const double n = static_cast<double>(mat.rows());
        const Eigen::VectorXd mean = mat.colwise().mean();
        // Center once and reuse the result to avoid building multiple large temporaries.
        Eigen::MatrixXd centered = mat.rowwise() - mean.transpose();
        Eigen::VectorXd stdDev = (centered.array().square().colwise().sum() / n).sqrt();

        Eigen::MatrixXd correlationMatrix = (centered.transpose() * centered) / n;
        Eigen::ArrayXd invStd = Eigen::ArrayXd::Zero(stdDev.size());
        for (Eigen::Index i = 0; i < stdDev.size(); ++i) {
            if (stdDev(i) > 0.0) {
                invStd(i) = 1.0 / stdDev(i);
            }
        }

        correlationMatrix =
            invStd.matrix().asDiagonal() * correlationMatrix * invStd.matrix().asDiagonal();

        // Zero-variance columns have undefined off-diagonal correlations; keep
        // their rows/columns at 0 but preserve a unit diagonal entry.
        for (Eigen::Index i = 0; i < stdDev.size(); ++i) {
            if (stdDev(i) == 0.0) {
                correlationMatrix(i, i) = 1.0;
            }
        }

        return correlationMatrix;
}


/**
 * @brief Remove a set of columns from a matrix in place.
 *
 * The input index vector may be unsorted and may contain duplicates. This
 * function normalizes the index list first, validates the indices, builds the
 * list of columns to keep, and copies the kept columns into a reduced matrix.
 * That avoids repeated left-shifts and repeated `conservativeResize` calls.
 *
 * Used in:
 * - `utils/phen.cpp` when dropping covariate columns flagged for removal
 * - `test/eigenmatrix_utils_example.cpp` as a small runnable example
 */
void remove_col(MatrixXd& mat, const vector<std::int64_t>& col_to_remove_in){
    if (col_to_remove_in.empty()) {
        return;
    }

    vector<std::int64_t> col_to_remove = col_to_remove_in;
    std::sort(col_to_remove.begin(), col_to_remove.end());
    col_to_remove.erase(std::unique(col_to_remove.begin(), col_to_remove.end()), col_to_remove.end());

    for (std::int64_t idx : col_to_remove) {
        eigen_assert(idx >= 0 && idx < mat.cols() && "Column index out of range");
    }

    std::vector<Eigen::Index> keep_cols;
    keep_cols.reserve(mat.cols() - static_cast<Eigen::Index>(col_to_remove.size()));

    Eigen::Index remove_pos = 0;
    for (Eigen::Index col = 0; col < mat.cols(); ++col) {
        if (remove_pos < static_cast<Eigen::Index>(col_to_remove.size()) &&
            col == static_cast<Eigen::Index>(col_to_remove[remove_pos])) {
            ++remove_pos;
        } else {
            keep_cols.push_back(col);
        }
    }

    Eigen::MatrixXd reduced(mat.rows(), keep_cols.size());
    for (Eigen::Index new_col = 0; new_col < static_cast<Eigen::Index>(keep_cols.size()); ++new_col) {
        reduced.col(new_col) = mat.col(keep_cols[new_col]);
    }

    mat.swap(reduced);
}

/**
 * @brief Scale each matrix column by the matching entry in a vector.
 *
 * This is an in-place column-wise multiply:
 * `mat.col(j) *= vec(j)` for all columns `j`.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp` when applying SNP weights to genotype blocks
 * - `test/eigenmatrix_utils_example.cpp` as a small runnable example
 */
void scale_cols_inplace(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec) {
    assert(mat.cols() == vec.size() && "Vector size must match the number of matrix columns!");
    // Eigen implements this column-wise scaling via row-wise broadcasting of
    // the transposed scaling vector.
    mat.array().rowwise() *= vec.transpose().array();
}


/**
 * @brief MKL/LAPACK-backed LLT helper for dense symmetric positive-definite matrices.
 *
 * The class stores the Cholesky factorization of the last matrix passed to
 * `compute()` and then reuses that factorization to derive quantities needed by
 * the mixed-model code, namely the log-determinant and the dense inverse.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mmsusie.cpp`
 * - `gmatrix/processGRM.cpp`
 * - `test/eigenmatrix_utils_example.cpp`
 */

void CustomLLT::compute(const Eigen::MatrixXd &A){
        info_computer = -1, info_inverse = -1;
        // LAPACK factorization overwrites the input buffer, so keep a private copy.
        A_tmp = A;
        ncols = A.cols();
        // Store the Cholesky factor in the lower triangle.
        info_computer = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols);
        eigen_assert(info_computer == 0 && "Fail to factorize. Matrix is not positive define");
}


/**
 * @brief Return log(det(A)) from the stored Cholesky factor.
 *
 * Requires a successful prior call to `compute()`.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mmsusie.cpp`
 * - `test/eigenmatrix_utils_example.cpp`
 */
double CustomLLT::logDeterminant(){
        eigen_assert(info_computer == 0 && "You must first call compute()");
        // For A = L L^T, log|A| = 2 * sum(log(diag(L))).
        double logdet = A_tmp.diagonal().array().square().log().sum();
        return logdet;
}

/**
 * @brief Return the dense inverse reconstructed from the stored Cholesky factor.
 *
 * Requires a successful prior call to `compute()`.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mmsusie.cpp`
 * - `gmatrix/processGRM.cpp`
 * - `test/eigenmatrix_utils_example.cpp`
 */
Eigen::MatrixXd CustomLLT::inverse(){
        eigen_assert(info_computer == 0 && "You must first call compute()");
        // Invert from the stored Cholesky factor.
        info_inverse = LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols);
        eigen_assert(info_inverse == 0 && "Fail to inverser");
        // dpotri fills only the lower triangle; copy to upper to recover the full symmetric inverse.
        A_tmp.triangularView<Eigen::Upper>() = A_tmp.transpose();
        return A_tmp;
}
