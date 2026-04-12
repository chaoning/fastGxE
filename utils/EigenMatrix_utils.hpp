#include <cstdint>


/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-06-22 21:08:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2026-03-31 11:54:26
 */
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "mkl.h"


/**
 * @brief Compute the column-wise correlation matrix.
 * @param mat Input matrix, with rows as samples and columns as variables.
 * @return Pearson correlation matrix between columns.
 */
Eigen::MatrixXd computeCorrelationMatrix(const Eigen::MatrixXd& mat);

/**
 * @brief Remove a set of columns from a matrix in place.
 * @param mat Matrix to modify.
 * @param col_to_remove Zero-based column indices to remove.
 */
void remove_col(Eigen::MatrixXd& mat, std::vector<std::int64_t> col_to_remove);

/**
 * @brief Scale matrix columns in place.
 * @param mat Matrix to scale.
 * @param vec Vector of column-wise scaling factors.
 */
void scale_cols_inplace(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec);


/**
 * @brief MKL/LAPACK-based LLT factorization helper for dense SPD matrices.
 * 
 */
class CustomLLT {

public:
    /**
     * @brief Factorize a dense symmetric positive-definite matrix.
     * @param A Input matrix.
     */
    void compute(const Eigen::MatrixXd &A);

    /**
     * @brief Return the log-determinant of the factored matrix.
     */
    double logDeterminant();

    /**
     * @brief Return the dense inverse from the stored factorization.
     */
    Eigen::MatrixXd inverse();

private:
    /// Internal factorization buffer.
    Eigen::MatrixXd A_tmp;
    /// Matrix dimension used by LAPACK.
    std::int64_t ncols;
    /// LAPACK status codes.
    lapack_int info_computer, info_inverse;
};
