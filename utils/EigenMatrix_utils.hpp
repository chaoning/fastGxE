
/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-06-22 21:08:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-01 14:07:01
 */
#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "mkl.h"


void remove_col(Eigen::MatrixXd& mat, std::vector<long long> col_to_remove);

void mat_col_elementwise_dot_vec(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec);

void mat_row_elementwise_dot_vec(Eigen::MatrixXd& mat, const Eigen::VectorXd& vec);


/**
 * @brief perform LLT decomposition, log-transformed determinant, solve, inverse
 * 
 */
class CustomLLT {

public:
    
    void compute(Eigen::MatrixXd &A);
    
    double logDeterminant();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b);

    Eigen::MatrixXd inverse();

    void pdelete();

private:
    Eigen::MatrixXd A_tmp;
    long long ncols;
    lapack_int info_computer, info_inverse, info_solve;
    Eigen::VectorXd diag;
};


/**
 * @brief perform LDLT decomposition, log-transformed determinant, solve, inverse
 * 
 */
class CustomLDLT {

public:
    
    void compute(Eigen::MatrixXd &A);
    
    double logDeterminant();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b);

    Eigen::MatrixXd inverse();

    void pdelete();

private:
    Eigen::MatrixXd A_tmp;
    long long ncols;
    lapack_int info_computer, info_inverse, info_solve;
    Eigen::VectorXd diag;
    lapack_int *ipiv;
};



/*********************/


class CustomPardisoLDLT{
public:
    CustomPardisoLDLT();

    ~CustomPardisoLDLT();

    void analyze(Eigen::SparseMatrix<double> &A, MKL_INT *_ia, MKL_INT *_ja);


    void factorize(Eigen::SparseMatrix<double> &A, MKL_INT *_ia, MKL_INT *_ja);


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b, MKL_INT *_ia, MKL_INT *_ja);


    Eigen::VectorXd getdiag();



    Eigen::VectorXd getdiag2();

    void pdelete(MKL_INT *_ia, MKL_INT *_ja);


private:
    mutable void *pt[64];
    MKL_INT maxfct;
    MKL_INT mnum;
    MKL_INT mtype;
    MKL_INT phase;
    MKL_INT n;
    Eigen::SparseMatrix<double> A_tmp;
    MKL_INT idum;
    MKL_INT nrhs;
    Eigen::VectorXi m_iparm;
    MKL_INT msglvl;
    double ddum;
    MKL_INT error;
    MKL_INT info_analyze, info_factorize, info_solve, info_getdiag, info_pdelete;
};


/*********************/


class CustomPardisoLLT{
public:
    CustomPardisoLLT();

    ~CustomPardisoLLT();

    void analyze(Eigen::SparseMatrix<double> &A);


    void factorize(Eigen::SparseMatrix<double> &A);


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b);


    Eigen::VectorXd getdiag();


    void pdelete();


private:
    mutable void *pt[64];
    MKL_INT maxfct;
    MKL_INT mnum;
    MKL_INT mtype;
    MKL_INT phase;
    MKL_INT n;
    Eigen::SparseMatrix<double> A_tmp;
    MKL_INT idum;
    MKL_INT nrhs;
    Eigen::VectorXi m_iparm;
    MKL_INT msglvl;
    double ddum;
    MKL_INT error;
    MKL_INT info_analyze, info_factorize, info_solve, info_getdiag, info_pdelete;
};


Eigen::MatrixXd custom_eigen_kron(Eigen::MatrixXd& A, Eigen::MatrixXd& B);


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::MatrixXd& A, Eigen::SparseMatrix<double>& B);


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B);


Eigen::SparseMatrix<double> custom_eigen_kron(Eigen::SparseMatrix<double>& A, Eigen::MatrixXd& B);



Eigen::MatrixXd custom_matrix_dot(Eigen::MatrixXd& A, Eigen::MatrixXd& B);
