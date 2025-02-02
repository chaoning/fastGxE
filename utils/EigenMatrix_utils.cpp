
/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-06-22 21:08:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-31 19:35:40
 */


#include <iostream>
#include <cassert>
#include <set>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>

#include "EigenMatrix_utils.hpp"
#include "mkl.h"


using namespace Eigen;
using namespace std;


/**
 * @brief remove given column NO.
 * 
 * @param mat 
 * @param col_to_remove 
 */
void remove_col(MatrixXd& mat, vector<long long> col_to_remove){
    long long nrows = mat.rows();
    long long ncols;
    set<long long> st(col_to_remove.begin(), col_to_remove.end());
    col_to_remove.assign(st.begin(), st.end()); // remove duplication, sort
    for(long long i = col_to_remove.size()-1; i >= 0; i--){
        ncols = mat.cols() - 1;
        if(col_to_remove[i] < ncols){
            mat.block(0, col_to_remove[i], nrows, ncols - col_to_remove[i]) = 
                mat.block(0, col_to_remove[i] + 1, nrows, ncols - col_to_remove[i]);
        }
        mat.conservativeResize(nrows, ncols);
    }
}



void mat_col_elementwise_dot_vec(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec) {
    assert(mat.rows() == vec.size() && "Vector size must match the number of matrix columns!");
    mat.array().colwise() *= vec.array();
}


void mat_row_elementwise_dot_vec(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec) {
    assert(mat.cols() == vec.size() && "Vector size must match the number of matrix rows!");
    mat.array().rowwise() *= vec.transpose().array();  // 使用 vec 转置以进行逐行广播
}


/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-06-22 21:08:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-30 23:35:40
 */


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_mat;


/**
 * @brief dense perform LLT decomposition, log-transformed determinant, solve, inverse
 * 
 */

void CustomLLT::compute(Eigen::MatrixXd &A){
        info_computer = -1, info_inverse = -1, info_solve = -1;
        A_tmp = A;
        diag = A.diagonal();
        ncols = A.cols();
        info_computer = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols);
        eigen_assert(info_computer == 0 && "Fail to factorize. Matrix is not positive define");
}


double CustomLLT::logDeterminant(){
        eigen_assert(info_computer == 0 && "You must first call compute()");
        double logdet = A_tmp.diagonal().array().square().log().sum();
        return logdet;
}


eigen_mat CustomLLT::solve(eigen_mat b){
        int ncols_b = b.cols();
        long long nrows_b = b.rows();

        eigen_assert(info_computer == 0 && "You must first call compute()");
        eigen_assert(info_inverse != 0 && "You cannot call inverse() before solve()");
        info_solve = LAPACKE_dpotrs (LAPACK_COL_MAJOR , 'L', ncols , ncols_b, A_tmp.data(), ncols, b.data(), nrows_b);
        eigen_assert(info_solve == 0 && "Fail to solve");
        return b;
}


Eigen::MatrixXd CustomLLT::inverse(){
        eigen_assert(info_computer == 0 && "You must first call compute()");
        info_inverse = LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols);
        eigen_assert(info_inverse == 0 && "Fail to inverser");
        A_tmp.triangularView<Eigen::Upper>() = A_tmp.transpose();
        return A_tmp;
}

void CustomLLT::pdelete(){
        A_tmp.resize(0, 0);
}



/**
 * @brief dense perform LDLT decomposition, log-transformed determinant, solve, inverse
 * 
 */

void CustomLDLT::compute(Eigen::MatrixXd &A){
        info_computer = -1, info_inverse = -1, info_solve = -1;
        A_tmp = A;
        diag = A.diagonal();
        ncols = A.cols();
        ipiv = new lapack_int[ncols];
        info_computer = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols, ipiv);
        eigen_assert(info_computer == 0 && "Fail to factorize");
}


double CustomLDLT::logDeterminant(){
        eigen_assert(info_computer == 0 && "You must first call compute()");
        double logdet = A_tmp.diagonal().array().log().sum();
        return logdet;
}


eigen_mat CustomLDLT::solve(eigen_mat b){
        int ncols_b = b.cols();
        long long nrows_b = b.rows();

        eigen_assert(info_computer == 0 && "You must first call compute()");
        eigen_assert(info_inverse != 0 && "You cannot call inverse() before solve()");

        info_solve = LAPACKE_dsytrs (LAPACK_COL_MAJOR , 'L', ncols , ncols_b, A_tmp.data(), ncols, ipiv, b.data(), nrows_b);
        eigen_assert(info_solve == 0 && "Fail to solve");

        return b;
}


Eigen::MatrixXd CustomLDLT::inverse(){
        eigen_assert(info_computer == 0 && "You must first call compute()");

        info_inverse = LAPACKE_dsytri(LAPACK_COL_MAJOR, 'L', ncols, A_tmp.data(), ncols, ipiv);
        eigen_assert(info_inverse == 0 && "Fail to inverser");
        A_tmp.triangularView<Eigen::Upper>() = A_tmp.transpose();

        return A_tmp;
}

void CustomLDLT::pdelete(){
        delete[] ipiv;
        A_tmp.resize(0, 0);
}

/************************/

typedef Eigen::SparseMatrix<double> SpMat;

/**
 * @brief Pardiso LDLT
 * 
 */

CustomPardisoLDLT::CustomPardisoLDLT(){
        for (int i = 0; i < 64; i++ ){
            pt[i] = 0;
        }
        maxfct = 1;        /* Maximum number of numerical factorizations. */
        mnum = 1;         /* Which factorization to use. */
        mtype = -2;       /* Real symmetric matrix */
        phase = 0;
        nrhs = 1;


        m_iparm = Eigen::VectorXi::Zero(64);
        m_iparm[0] = 1;         /* No solver default */
        m_iparm[1] = 2;         /* Fill-in reordering from METIS */
        m_iparm[3] = 0;         /* No iterative-direct algorithm */
        m_iparm[4] = 0;         /* No user fill-in reducing permutation */
        m_iparm[5] = 0;         /* Write solution into x */
        m_iparm[6] = 0;         /* Not in use */
        m_iparm[7] = 2;         /* Max numbers of iterative refinement steps */
        m_iparm[8] = 0;         /* Not in use */
        m_iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
        m_iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
        m_iparm[11] = 0;        /* Not in use */
        m_iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
        m_iparm[13] = 0;        /* Output: Number of perturbed pivots */
        m_iparm[14] = 0;        /* Not in use */
        m_iparm[15] = 0;        /* Not in use */
        m_iparm[16] = 0;        /* Not in use */
        m_iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
        m_iparm[18] = -1;       /* Output: Mflops for LU factorization */
        m_iparm[19] = 0;        /* Output: Numbers of CG Iterations */
        m_iparm[27] = 0;        // double precision
        m_iparm[34] = 1;       // Zero-based indexing: columns and rows indexing in arrays ia, ja, and perm starts from 0 (C-style indexing)
        m_iparm[55] = 1;       // Diagonal and pivoting control.

        
        msglvl = 0;           /* prints statistical information to the screen */
        error = 0;            /* Initialize error flag */
        info_analyze = -10, info_factorize = -10, info_solve = -10, info_getdiag = -10, info_pdelete = -10;
}

CustomPardisoLDLT::~CustomPardisoLDLT(){
        // eigen_assert(info_pdelete == 0 && "You must first call pdelete()");
}


void CustomPardisoLDLT::analyze(SpMat &A, MKL_INT *_ia, MKL_INT *_ja){
        // PARDISO supports only upper, row-major or lower, column major
        A_tmp = A.triangularView<Eigen::Lower>(); // only keep lower part for column major
        A_tmp.makeCompressed();
        n = A_tmp.rows();
        /* .. Reordering and Symbolic Factorization. This step also allocates all memory that is necessary for the factorization. */
        phase = 11;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), _ia, _ja, &idum, &nrhs, m_iparm.data(), &msglvl, &ddum, &ddum, &info_analyze);
        eigen_assert(info_analyze == 0 && "Fail to analyze the matrix");
}


void CustomPardisoLDLT::factorize(SpMat &A, MKL_INT *_ia, MKL_INT *_ja){
        // PARDISO supports only upper, row-major or lower, column major
        A_tmp = A.triangularView<Eigen::Lower>(); // only keep lower part for column major
        A_tmp.makeCompressed();
        n = A_tmp.rows();
        /* Numerical factorization. */
        phase = 22;
        eigen_assert(info_analyze == 0 && "You must first call analyze()");
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), _ia, _ja, &idum, &nrhs, m_iparm.data(), &msglvl, &ddum, &ddum, &info_factorize);
        eigen_assert(info_factorize == 0 && "Fail to factorize");
}


eigen_mat CustomPardisoLDLT::solve(eigen_mat b, MKL_INT *_ia, MKL_INT *_ja){
        eigen_assert(info_factorize == 0 && "You must first call factorize()");
        phase = 33;
        m_iparm[7] = 2;
        int nrhs = b.cols();
        eigen_mat x(n, nrhs);
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), _ia, _ja, &idum, &nrhs, m_iparm.data(), &msglvl, b.data(), x.data(), &info_solve);
        eigen_assert(info_solve == 0 && "Fail to solve");
        return x;
}


Eigen::VectorXd CustomPardisoLDLT::getdiag(){
        eigen_assert(info_factorize == 0 && "You must first call factorize()");
        Eigen::VectorXd df = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd da = Eigen::VectorXd::Zero(n);
        pardiso_getdiag (pt, df.data(), da.data(), &mnum, &info_getdiag);
        eigen_assert(info_getdiag == 0 && "Fail to getdiag");
        return df;
}


/*
Eigen::VectorXd CustomPardisoLDLT::getdiag2(){
        eigen_assert(info_factorize == 0 && "You must first call factorize()");
        m_iparm[7] = 0;
        phase = 332;
        Eigen::VectorXd b_tmp = Eigen::VectorXd::Ones(n);
        Eigen::VectorXd invdiag = Eigen::VectorXd::Zero(n);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
                &n, A_tmp.valuePtr(), A_tmp.outerIndexPtr(), A_tmp.innerIndexPtr(), &idum, &nrhs, m_iparm.data(), &msglvl, b_tmp.data(), invdiag.data(), &info_getdiag);
        eigen_assert(info_getdiag == 0 && "Fail to getdiag");
        invdiag = 1/invdiag.array();
        return invdiag;
}
*/


void CustomPardisoLDLT::pdelete(MKL_INT *_ia, MKL_INT *_ja){
        phase = -1;           /* Release internal memory. */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, _ia, _ja, &idum, &nrhs,
             m_iparm.data(), &msglvl, &ddum, &ddum, &info_pdelete);
        eigen_assert(info_pdelete == 0 && "Fail to release mermory");
        A_tmp.resize(0, 0);
}




/************************/


/**
 * @brief Pardiso LLT
 * 
 */

CustomPardisoLLT::CustomPardisoLLT(){
        for (int i = 0; i < 64; i++ ){
            pt[i] = 0;
        }
        maxfct = 1;        /* Maximum number of numerical factorizations. */
        mnum = 1;         /* Which factorization to use. */
        mtype = 2;       /* real and symmetric positive definite */
        phase = 0;
        nrhs = 1;


        m_iparm = Eigen::VectorXi::Zero(64);
        m_iparm[0] = 1;         /* No solver default */
        m_iparm[1] = 2;         /* Fill-in reordering from METIS */
        m_iparm[3] = 0;         /* No iterative-direct algorithm */
        m_iparm[4] = 0;         /* No user fill-in reducing permutation */
        m_iparm[5] = 0;         /* Write solution into x */
        m_iparm[6] = 0;         /* Not in use */
        m_iparm[7] = 2;         /* Max numbers of iterative refinement steps */
        m_iparm[8] = 0;         /* Not in use */
        m_iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
        m_iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
        m_iparm[11] = 0;        /* Not in use */
        m_iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
        m_iparm[13] = 0;        /* Output: Number of perturbed pivots */
        m_iparm[14] = 0;        /* Not in use */
        m_iparm[15] = 0;        /* Not in use */
        m_iparm[16] = 0;        /* Not in use */
        m_iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
        m_iparm[18] = -1;       /* Output: Mflops for LU factorization */
        m_iparm[19] = 0;        /* Output: Numbers of CG Iterations */
        m_iparm[27] = 0;        // double precision
        m_iparm[34] = 1;       // Zero-based indexing: columns and rows indexing in arrays ia, ja, and perm starts from 0 (C-style indexing)
        m_iparm[55] = 1;       // Diagonal and pivoting control.

        
        msglvl = 0;           /* prints statistical information to the screen */
        error = 0;            /* Initialize error flag */
        info_analyze = -10, info_factorize = -10, info_solve = -10, info_getdiag = -10, info_pdelete = -10;
}

CustomPardisoLLT::~CustomPardisoLLT(){
        // eigen_assert(info_pdelete == 0 && "You must first call pdelete()");
}


void CustomPardisoLLT::analyze(SpMat &A){
        // PARDISO supports only upper, row-major or lower, column major
        A_tmp = A;
        A_tmp = A_tmp.triangularView<Eigen::Lower>(); // only keep lower part for column major
        A_tmp.makeCompressed();
        n = A_tmp.rows();
        /* .. Reordering and Symbolic Factorization. This step also allocates all memory that is necessary for the factorization. */
        phase = 11;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), A_tmp.outerIndexPtr(), A_tmp.innerIndexPtr(), &idum, &nrhs, m_iparm.data(), &msglvl, &ddum, &ddum, &info_analyze);
        eigen_assert(info_analyze == 0 && "Fail to analyze the matrix");
}


void CustomPardisoLLT::factorize(SpMat &A){
        // PARDISO supports only upper, row-major or lower, column major
        A_tmp = A;
        A_tmp = A_tmp.triangularView<Eigen::Lower>(); // only keep lower part for column major
        A_tmp.makeCompressed();
        n = A_tmp.rows();
        /* Numerical factorization. */
        phase = 22;
        eigen_assert(info_analyze == 0 && "You must first call analyze()");
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), A_tmp.outerIndexPtr(), A_tmp.innerIndexPtr(), &idum, &nrhs, m_iparm.data(), &msglvl, &ddum, &ddum, &info_factorize);
        eigen_assert(info_factorize == 0 && "Fail to factorize");
}


eigen_mat CustomPardisoLLT::solve(eigen_mat b){
        eigen_assert(info_factorize == 0 && "You must first call factorize()");
        phase = 33;
        m_iparm[7] = 2;
        int nrhs = b.cols();
        eigen_mat x(n, nrhs);
        pardiso (pt, &maxfct, &mnum, &mtype, &phase,
                &n, A_tmp.valuePtr(), A_tmp.outerIndexPtr(), A_tmp.innerIndexPtr(), &idum, &nrhs, m_iparm.data(), &msglvl, b.data(), x.data(), &info_solve);
        eigen_assert(info_solve == 0 && "Fail to solve");
        return x;
}


Eigen::VectorXd CustomPardisoLLT::getdiag(){
        eigen_assert(info_factorize == 0 && "You must first call factorize()");
        Eigen::VectorXd df = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd da = Eigen::VectorXd::Zero(n);
        pardiso_getdiag (pt, df.data(), da.data(), &mnum, &info_getdiag);
        eigen_assert(info_getdiag == 0 && "Fail to getdiag");
        return df;
}



void CustomPardisoLLT::pdelete(){
        phase = -1;           /* Release internal memory. */
        pardiso(pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, A_tmp.outerIndexPtr(), A_tmp.innerIndexPtr(), &idum, &nrhs,
             m_iparm.data(), &msglvl, &ddum, &ddum, &info_pdelete);
        eigen_assert(info_pdelete == 0 && "Fail to release mermory");
        A_tmp.resize(0, 0);
}










/**
 * @brief KroneckerProduct
 * 
 * @param A 
 * @param B 
 * @return 
 */
Eigen::MatrixXd custom_eigen_kron(Eigen::MatrixXd& A, Eigen::MatrixXd& B){
    typedef Eigen::KroneckerProduct<Eigen::MatrixXd, Eigen::MatrixXd> kron_dense_dense;
    Eigen::MatrixXd ret = kron_dense_dense(A, B);
    return ret;
}


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::MatrixXd& A, Eigen::SparseMatrix<double>& B){
    typedef Eigen::KroneckerProductSparse<Eigen::MatrixXd, Eigen::SparseMatrix<double>> kron_dense_sparse;
    Eigen::SparseMatrix<double> ret = kron_dense_sparse(A, B);
    return ret;
}


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B){
    typedef Eigen::KroneckerProductSparse<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> kron_sparse_sparse;
    Eigen::SparseMatrix<double> ret = kron_sparse_sparse(A, B);
    return ret;
}


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& B){
    typedef Eigen::KroneckerProductSparse<Eigen::SparseMatrix<double>, Eigen::MatrixXd> kron_sparse_dense;
    Eigen::SparseMatrix<double> ret = kron_sparse_dense(A, B);
    return ret;
}



/**
 * @brief  matrix dot
 * 
 */

Eigen::MatrixXd custom_matrix_dot(Eigen::MatrixXd& A, Eigen::MatrixXd& B){
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(A.rows(), B.cols());
        cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, A.rows(), B.cols(), A.cols(), 1, 
                A.data(), A.rows(), 
                B.data(), B.rows(), 
                0, C.data(), A.rows());
        return C;
}
