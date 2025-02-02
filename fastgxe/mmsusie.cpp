/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2023-09-08 11:17:09
 * @LastEditTime: 2025-02-01 22:14:58
 * @LastEditors: Chao Ning
 */
#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen



#include <getopt.h>
#include "mkl.h"
#include "omp.h"
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <gsl/gsl_cdf.h>
#include <LBFGSB.h>

#include <CLI/CLI.hpp>  // Include CLI11
#include "../utils/geno.hpp"
#include "../gmatrix/processGRM.hpp"
#include "../utils/phen.hpp"
#include "../utils/EigenMatrix_utils.hpp"
#include "../utils/chi2score.hpp"
#include "../utils/string_utils.hpp"
#include "../utils/iterator_utils.hpp"


using std::set;
using std::map;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::endl;
using Eigen::HouseholderQR;
using Eigen::ColPivHouseholderQR;
using namespace LBFGSpp;

#include "mmsusie.hpp"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


MMSUSIE::MMSUSIE(){

}

MMSUSIE::~MMSUSIE() {
}


class optimize_prior_variance
{
private:
    int n;
    VectorXd betahat_Vec, shat2_Vec, prior_weights_Vec;
public:
    optimize_prior_variance(int n_) : n(n_) {}
    void set_value(VectorXd betahat_Vec_, VectorXd shat2_Vec_, VectorXd prior_weights_Vec_){
        betahat_Vec = betahat_Vec_;
        shat2_Vec = shat2_Vec_;
        prior_weights_Vec = prior_weights_Vec_;
    }

    double neg_loglik(double V0){
        int p = betahat_Vec.size();
        double fx = 0.0;
        VectorXd zscore2 = betahat_Vec.array() * betahat_Vec.array() / shat2_Vec.array();
        VectorXd lbf_Vec = 0.5 * (shat2_Vec.array() / (shat2_Vec.array() + V0)).log() 
                +  0.5 * zscore2.array() * V0 / (shat2_Vec.array() + V0);
        
        double maxlbf = lbf_Vec.maxCoeff();
        VectorXd w_Vec = (lbf_Vec.array() - maxlbf).exp();
        VectorXd w_weighted_Vec = w_Vec.cwiseProduct(prior_weights_Vec);
        double weighted_sum_w = (w_weighted_Vec).sum();
        fx = -(std::log(weighted_sum_w) + maxlbf);
        return fx;
    }
    
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        int p = betahat_Vec.size();
        double fx = 0.0;
        VectorXd zscore2 = betahat_Vec.array() * betahat_Vec.array() / shat2_Vec.array();
        VectorXd lbf_Vec = 0.5 * (shat2_Vec.array() / (shat2_Vec.array() + x(0))).log() 
                +  0.5 * zscore2.array() * x(0) / (shat2_Vec.array() + x(0));
        
        double maxlbf = lbf_Vec.maxCoeff();
        VectorXd w_Vec = (lbf_Vec.array() - maxlbf).exp();
        VectorXd w_weighted_Vec = w_Vec.cwiseProduct(prior_weights_Vec);
        double weighted_sum_w = (w_weighted_Vec).sum();
        VectorXd alpha_Vec = w_weighted_Vec / weighted_sum_w;
        VectorXd lgrad_Vec = 0.5 * 1 / (shat2_Vec.array() + x(0)) * (zscore2.array() * shat2_Vec.array() / (shat2_Vec.array() + x(0)) - 1); 
        grad(0) = -(alpha_Vec.cwiseProduct(lgrad_Vec)).sum();
        fx = -(std::log(weighted_sum_w) + maxlbf);
        return fx;
    }
};

class optimize_sigma{
private:
    int n;
    VectorXd y, Xr_Vec;
    MatrixXd alpha_Mat, mu_Mat, mu2_Mat, X;
    double vmat_logdet;
    SparseMatrix < double > Vi;
    SparseMatrix < double > gxe_rm, grm_mat;
    long long num_id;

    vector<MatrixXd> grm_mat_group_vec;
    vector<long long> grm_index_vec;
    MatrixXd bye_mat;

public:
    optimize_sigma(int n_) : n(n_) {}

    void set_value(VectorXd& y_, MatrixXd& X_, SparseMatrix < double >& gxe_rm_, SparseMatrix < double >& grm_mat_,
        VectorXd& Xr_Vec_, MatrixXd& alpha_Mat_, MatrixXd& mu_Mat_, MatrixXd& mu2_Mat_,
        vector<MatrixXd>& grm_mat_group_vec_, vector<long long>& grm_index_vec_,
                    MatrixXd& bye_mat_){
        y = y_;
        num_id = y.size();
        X = X_;

        gxe_rm = gxe_rm_;
        grm_mat = grm_mat_;
        
        Xr_Vec = Xr_Vec_;
        alpha_Mat = alpha_Mat_;
        mu_Mat = mu_Mat_;
        mu2_Mat = mu2_Mat_;

        grm_mat_group_vec = grm_mat_group_vec_;
        grm_index_vec = grm_index_vec_;
        bye_mat = bye_mat_;
    }

    double get_V_logdet(){
        return vmat_logdet;
    }
    SparseMatrix < double > get_Vi(){
        return Vi;
    }

    void calVi(const VectorXd& varcom){
        std::vector<Eigen::Triplet<double>> tripletList;
        long long num_bye = bye_mat.cols();
        long long num_id = bye_mat.rows();
        CustomLLT llt_solver;
        vmat_logdet = 0;
        for(long long i = 0; i < grm_mat_group_vec.size(); i++){
            MatrixXd tmp_grm_mat = grm_mat_group_vec[i];
            long long start_index = grm_index_vec[i];
            long long num_element = tmp_grm_mat.rows();
            MatrixXd bye_mat_part = bye_mat.middleRows(start_index, num_element);
            if(i == 0){
                if(tmp_grm_mat.size() != 0){
                    MatrixXd tmp_gxe_mat = (bye_mat_part.cwiseProduct(bye_mat_part)).rowwise().sum() / num_bye;
                    tmp_gxe_mat = (tmp_grm_mat.cwiseProduct(tmp_gxe_mat)).eval();
                    Eigen::ArrayXXd tmp_mat = tmp_grm_mat.array() * varcom(0) + 
                                tmp_gxe_mat.array() * varcom(1) + varcom(2);
                    vmat_logdet += tmp_mat.log().sum();
                    tmp_mat = 1 / tmp_mat.array();
                    for(long long m = 0; m < num_element; m++){
                        tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat(m, 0)));
                    }
                }
            }else{
                MatrixXd tmp_gxe_mat = bye_mat_part * bye_mat_part.transpose() / num_bye;
                tmp_gxe_mat = (tmp_grm_mat.cwiseProduct(tmp_gxe_mat)).eval();

                MatrixXd tmp_mat = MatrixXd::Identity(num_element, num_element) * varcom(2);
                tmp_mat += tmp_grm_mat * varcom(0) + tmp_gxe_mat * varcom(1);
                llt_solver.compute(tmp_mat);
                vmat_logdet += llt_solver.logDeterminant();
                tmp_mat = llt_solver.inverse();
                for(long long m = 0; m < num_element; m++){
                    for(long long n = 0; n < num_element; n++)
                        tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat(m, n)));
                }
            }
        }
        
        Vi.resize(num_id, num_id);
        Vi.setZero();
        Vi.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.shrink_to_fit();
    }
    

    SparseMatrix < double > calGxEGRM(vector<MatrixXd>& grm_mat_group_vec, vector<long long>& grm_index_vec,
                    MatrixXd& bye_mat){
        long long num_bye = bye_mat.cols();
        long long num_id = bye_mat.rows();
        std::vector<Eigen::Triplet<double>> tripletList;
        vector<MatrixXd> gxe_rm_group_vec(grm_mat_group_vec.size());
        for(long long i = 0; i < grm_mat_group_vec.size(); i++){
            MatrixXd tmp_mat = grm_mat_group_vec[i];
            gxe_rm_group_vec[i] = tmp_mat; // Elements index 0 will not be null
            long long start_index = grm_index_vec[i];
            long long num_element = tmp_mat.rows();
            MatrixXd bye_mat_part = bye_mat.middleRows(start_index, num_element);
            if(i == 0){
                if(tmp_mat.size() != 0){
                    MatrixXd tmp_mat2 = (bye_mat_part.cwiseProduct(bye_mat_part)).rowwise().sum() / num_bye;
                    tmp_mat2 = (tmp_mat2.cwiseProduct(tmp_mat)).eval();
                    gxe_rm_group_vec[i] = tmp_mat2;
                    for(long long m = 0; m < num_element; m++){
                        tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat2(m, 0)));
                    }
                }
            }else{
                MatrixXd tmp_mat2 = bye_mat_part * bye_mat_part.transpose() / num_bye;
                tmp_mat2 = (tmp_mat2.cwiseProduct(tmp_mat)).eval();
                gxe_rm_group_vec[i] = tmp_mat2;
                for(long long m = 0; m < num_element; m++){
                    for(long long n = 0; n < num_element; n++)
                        tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat2(m, n)));
                }
            }
        }
        tripletList.shrink_to_fit();
        SparseMatrix < double > gxe_rm(num_id, num_id);
        gxe_rm.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.shrink_to_fit();
        return gxe_rm;
    }

    double neg_loglik(){
        long long p = X.cols();
        VectorXd X_xtVix(p);
        for(long long i = 0; i < p; i++){
            X_xtVix(i) = X.col(i).transpose() * this->Vi * X.col(i);
        }
        double fx = - 0.5 * n * std::log(2 * M_PI) - 0.5 * this->vmat_logdet 
                - 0.5 * ( (y - Xr_Vec).transpose() * this->Vi * (y - Xr_Vec) + 
                    (alpha_Mat.cwiseProduct(mu2_Mat).colwise().sum()).cwiseProduct(X_xtVix.transpose()).sum() -
                    (alpha_Mat.cwiseProduct(mu_Mat).array().square().matrix().colwise().sum()).cwiseProduct(X_xtVix.transpose()).sum() 
                    );
        return -fx;
    }

    void neg_gradVec(VectorXd& grad){
        long long p = X.cols();
        grad.resize(3);
        for(long long i = 0; i < 3; i++){
            SparseMatrix < double > spFD(num_id, num_id);
            if(i == 0)
                spFD = grm_mat;
            else if (i == 1)
                spFD = gxe_rm;
            else
                spFD.setIdentity();
            
            VectorXd X_xtVi3x(p);
            for(long long i = 0; i < p; i++){
                X_xtVi3x(i) = X.col(i).transpose() * this->Vi * spFD * this->Vi * X.col(i);
            }
            grad(i) = -0.5 * this->Vi.cwiseProduct(spFD).sum() 
                + 0.5 * (y - Xr_Vec).transpose() * this->Vi * spFD * this->Vi * (y - Xr_Vec) + 
                0.5 * (alpha_Mat.cwiseProduct(mu2_Mat).colwise().sum()).cwiseProduct(X_xtVi3x.transpose()).sum() -
                0.5 * (alpha_Mat.cwiseProduct(mu_Mat).array().square().matrix().colwise().sum()).cwiseProduct(X_xtVi3x.transpose()).sum() ;
        }
        grad = -grad;
    }

    double operator()(const VectorXd& x, VectorXd& grad){
        calVi(x);
        double fx = neg_loglik();
        neg_gradVec(grad);
        return fx;
    }

};


void MMSUSIE::pre_data_GxE(bool standardize_env, int phen_correct){
    spdlog::info("Standardize the interaction environment covariates");
    // rank of bye_mat
    ColPivHouseholderQR<MatrixXd> qr(this->m_bye_mat);
    if(qr.rank() < m_num_bye){
        spdlog::error("There are dependent columns in the interaction environment covariates!");
        exit(1);
    }

    if(standardize_env){
        // scale to (0, 1)
        VectorXd meanVec = this->m_bye_mat.colwise().mean();
        MatrixXd centered = this->m_bye_mat.rowwise() - meanVec.transpose();
        VectorXd stdVec = (centered.cwiseProduct(centered).colwise().sum().array() / centered.rows()).sqrt();
        this->m_bye_mat = centered.array().rowwise() * (1 / stdVec.transpose().array());
    }
    // std::cout << _bye_mat.transpose() << std::endl << std::endl;

    if(phen_correct == 1)
    {
        this->m_y = (this->m_y - m_xmat * ((m_xmat.transpose() * m_xmat).inverse() * (m_xmat.transpose() * m_y))).eval();
        long long nrows = this->m_xmat.rows();
        this->m_xmat = MatrixXd::Ones(nrows, 1);
        this->m_xmat.conservativeResize(nrows, 1 + this->m_bye_mat.cols());
        this->m_xmat.rightCols(this->m_bye_mat.cols()) = this->m_bye_mat; // interaction covariate
    }else if (phen_correct == 2)
    {
        long long nrows = this->m_xmat.rows();
        long long ncols = this->m_xmat.cols();
        this->m_xmat.conservativeResize(nrows, ncols + this->m_bye_mat.cols());
        this->m_xmat.rightCols(this->m_bye_mat.cols()) = this->m_bye_mat; // interaction covariate
        qr.compute(this->m_xmat);
        if(qr.rank() < m_xmat.cols()){
            spdlog::error("There are dependent columns between interaction environment covariates and other covariates!");
            exit(1);
        }

        this->m_y = (this->m_y - m_xmat * ((m_xmat.transpose() * m_xmat).inverse() * (m_xmat.transpose() * m_y))).eval();

        this->m_xmat = MatrixXd::Ones(nrows, 1);
    }else{
        spdlog::error("--phen-correct must be 1 or 2");
        exit(1);
    }
}



Eigen::MatrixXd MMSUSIE::pre_data_mmsusie(string bed_file, vector<string> snp_vec, int phen_correct){
    if(snp_vec.empty()){
        if(phen_correct == 1){
            return this->m_bye_mat;
        }else{
            spdlog::error("For empty SNPs, --phen-correct must be 1");
            exit(1);
        }
    }else{
        GENO GenoA(bed_file);
        vector<long long> id_index_in_bed_vec = GenoA.find_fam_index_pro(m_id_in_data_vec);
        vector<long long> snp_index_vec = GenoA.find_bim_index_pro(snp_vec);

        MatrixXd snp_mat;
        VectorXd freq_arr, missing_rate_arr;
        GenoA.read_bed_by_snp_index(bed_file + ".bed", snp_mat, freq_arr, missing_rate_arr, snp_index_vec, id_index_in_bed_vec);
        
        // scale to zero mean and 1 std
        snp_mat.rowwise() -= 2*freq_arr.transpose();
        VectorXd snp_std_Vec = ((snp_mat.cwiseProduct(snp_mat)).colwise().sum().array() / snp_mat.rows()).sqrt();
        snp_mat = snp_mat.array().rowwise() / snp_std_Vec.transpose().array();

        // marginal and GxE effects
        Eigen::MatrixXd X;
        if(phen_correct == 1){
            X.resize(id_index_in_bed_vec.size(), snp_index_vec.size() + this->m_bye_mat.cols() + snp_index_vec.size() * this->m_bye_mat.cols());
            X.leftCols(snp_mat.cols()) = snp_mat;
            X.middleCols(snp_mat.cols(), this->m_bye_mat.cols()) = this->m_bye_mat;
            long long k = snp_mat.cols() + this->m_bye_mat.cols();
            for(long long i = 0; i < snp_index_vec.size(); i++){
                for(long long j = 0; j < this->m_bye_mat.cols(); j++){
                    X.col(k) = snp_mat.col(i).cwiseProduct(this->m_bye_mat.col(j));
                    k++;
                }
            }
            spdlog::info("Number of susie variables: {}", k);
        }else if (phen_correct == 2)
        {
            X.resize(id_index_in_bed_vec.size(), snp_index_vec.size() + snp_index_vec.size() * this->m_bye_mat.cols());
            X.leftCols(snp_mat.cols()) = snp_mat;
            long long k = snp_mat.cols();
            for(long long i = 0; i < snp_index_vec.size(); i++){
                for(long long j = 0; j < this->m_bye_mat.cols(); j++){
                    X.col(k) = snp_mat.col(i).cwiseProduct(this->m_bye_mat.col(j));
                    k++;
                }
            }
            spdlog::info("Number of susie variables: {}", k);
        }
        return X;
    }
}


void MMSUSIE::mmsusiefun(MatrixXd X, long long L, long long maxiter, double tol, 
               double coverage, double min_abs_corr, bool estimate_sigma){
    this->mmsusiefun2(X, this->m_y, L, maxiter, tol, coverage, min_abs_corr, estimate_sigma);
}

/**
 * long long L: the number of pre-defined casual variables
*/
void MMSUSIE::mmsusiefun2(MatrixXd X, VectorXd y, long long L, long long maxiter, double tol, 
               double coverage, double min_abs_corr, bool estimate_sigma){
    long long p = X.cols(), n = X.rows();
    if(p < L) L = p;
    double varY = (y.array() - y.mean()).square().sum() / y.size();

    // X.col(i).T * Vi * X.col(i) for each element, simple linear regression
    VectorXd X_xtVix(p);
    for(long long i = 0; i < p; i++){
        X_xtVix(i) = X.col(i).transpose() * this->m_Vi0 * X.col(i);
    }
    VectorXd shat2_Vec = 1 / X_xtVix.array(); // (x'V(-1)x)(-1) for simple linear regression

    // Initialize susie fit
    VectorXd prior_weights_Vec = VectorXd::Ones(p) / p;
    MatrixXd alpha_Mat = MatrixXd::Ones(L, p) / p;
    MatrixXd mu_Mat = MatrixXd::Zero(L, p);
    MatrixXd mu2_Mat = MatrixXd::Zero(L, p);
    VectorXd Xr_Vec = VectorXd::Zero(n);
    VectorXd KL_Vec = VectorXd::Ones(L) * NAN;
    VectorXd lbf_Vec = VectorXd::Ones(L) * NAN;
    MatrixXd lbf_variable_Mat = MatrixXd::Ones(L, p) * NAN;
    VectorXd V0_Vec = VectorXd::Ones(L) * varY * 0.2;

    // Initialize elbo to NA.
    VectorXd elbo_Vec = VectorXd::Ones(maxiter + 1) * NAN;
    elbo_Vec[0] = -std::numeric_limits<double>::infinity();

    optimize_sigma optimize_sigmaA(3);
    SparseMatrix < double > gxe_rm;
    if(estimate_sigma)
        gxe_rm = optimize_sigmaA.calGxEGRM(this->m_grm_mat_group_vec, this->m_grm_index_vec, this->m_bye_mat);
    VectorXd varcom = this->m_varcom_null;
    SparseMatrix <double> Vi = this->m_Vi0;
    double V_logdet = this->m_V0_logdet;
    for(long long iter = 0; iter < maxiter; iter++){
        spdlog::info("Iteration: {}", iter + 1);
        // update each effect once
        for(long long l = 0; l < L; l++){
            // Remove lth effect from fitted values
            Xr_Vec = Xr_Vec - X * (alpha_Mat.row(l).cwiseProduct(mu_Mat.row(l))).transpose();

            // Compute residuals
            VectorXd R = y - Xr_Vec;

            /*** Bayesian single-effect linear regression using R outcomes ***/
            VectorXd XtViy_Vec = X.transpose() * (Vi * R);
            VectorXd betahat_Vec = shat2_Vec.array() * XtViy_Vec.array(); // betas for p simple least-squares
            
            // optimize the prior variance
            double V0 = V0_Vec(l);
            LBFGSBParam<double> param;
            LBFGSBSolver<double> solver(param);
            optimize_prior_variance fun(1);
            fun.set_value(betahat_Vec, shat2_Vec, prior_weights_Vec);
            VectorXd lb = VectorXd::Constant(1, 1.0e-10);
            VectorXd ub = VectorXd::Constant(1, 1.0e10);
            // Initial values
            VectorXd x = VectorXd::Constant(1, V0);
            double fx;
            int niter = solver.minimize(fun, x, fx, lb, ub);

            // std::cout << niter << " iterations" << std::endl;
            // std::cout << "x = " << x.transpose() << std::endl;
            // std::cout << "f(x) = " << fx << std::endl;
            // std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
            // std::cout << "projected grad norm = " << solver.final_grad_norm() << std::endl;
            if(fun.neg_loglik(x(0)) < fun.neg_loglik(V0)) V0 = x(0);
            V0_Vec(l) = V0;

            // V0_Vec << 0.0329366,  0.00108375, 0.000713443, 0.000129242, 0.000129264, 0.000129265, 0.000129248, 0.000129224, 0.000129208, 0.000129211;
            // V0 = V0_Vec(l);
            
            // double V0 = 0.03;
            // log(bf) for each SNP
            VectorXd zscore2 = betahat_Vec.array() * betahat_Vec.array() / shat2_Vec.array();
            VectorXd lbf_snp_Vec = 0.5 * (shat2_Vec.array() / (shat2_Vec.array() + V0)).log() 
                    +  0.5 * zscore2.array() * V0 / (shat2_Vec.array() + V0);
            
            // Posterior prob for each SNP.
            double maxlbf = lbf_snp_Vec.maxCoeff();
            VectorXd w_Vec = (lbf_snp_Vec.array() - maxlbf).exp();
            VectorXd w_weighted_Vec = w_Vec.cwiseProduct(prior_weights_Vec);
            double weighted_sum_w = (w_weighted_Vec).sum();
            VectorXd alpha_Vec = w_weighted_Vec / weighted_sum_w;

            VectorXd post_var = 1 / (1 / V0 + 1 / shat2_Vec.array()); // Posterior variance.
            VectorXd post_mean = betahat_Vec.array() / shat2_Vec.array() * post_var.array();
            VectorXd post_mean2 = post_var + post_mean.cwiseProduct(post_mean); // Second moment.

            // BF for single effect model.
            double lbf_model = maxlbf + std::log(weighted_sum_w);
            double loglik = lbf_model - 0.5 * n * std::log(2 * M_PI) - 0.5 * V_logdet 
                            - 0.5 * R.transpose() * Vi * R;
            
            // update
            mu_Mat.row(l) = post_mean;
            alpha_Mat.row(l) = alpha_Vec;
            mu2_Mat.row(l) = post_mean2;
            lbf_Vec(l) = lbf_model;
            lbf_variable_Mat.row(l) = lbf_snp_Vec;
            double SER_posterior_e_loglik = - 0.5 * n * std::log(2 * M_PI) - 0.5 * V_logdet 
                            - 0.5 * ( R.transpose() * Vi * R - 
                                      2 * (R.transpose() * Vi * X * alpha_Vec.cwiseProduct(post_mean)).sum() + 
                                      (X_xtVix.cwiseProduct(alpha_Vec.cwiseProduct(post_mean2))).sum() );
            KL_Vec(l) = -loglik + SER_posterior_e_loglik;

            Xr_Vec = Xr_Vec + X * (alpha_Mat.row(l).cwiseProduct(mu_Mat.row(l))).transpose();
        }

        spdlog::info("Estimated prior variances: {}", join_string(V0_Vec));

        elbo_Vec[iter + 1] = - 0.5 * n * std::log(2 * M_PI) - 0.5 * V_logdet 
                       - 0.5 * ( (y - Xr_Vec).transpose() * Vi * (y - Xr_Vec) + 
                                 (alpha_Mat.cwiseProduct(mu2_Mat).colwise().sum()).cwiseProduct(X_xtVix.transpose()).sum() -
                                 (alpha_Mat.cwiseProduct(mu_Mat).array().square().matrix().colwise().sum()).cwiseProduct(X_xtVix.transpose()).sum() 
                       ) - KL_Vec.sum();
        spdlog::info("ELBO: {}", elbo_Vec[iter + 1]);
        if(std::fabs(elbo_Vec[iter + 1] - elbo_Vec[iter]) < tol){
            break;
        }

        if(estimate_sigma){
            // estimate sigma
            optimize_sigma optimize_sigmaA(3);
            optimize_sigmaA.set_value(y, X, gxe_rm, m_grm_mat, Xr_Vec, alpha_Mat, mu_Mat, mu2_Mat,
                    m_grm_mat_group_vec, m_grm_index_vec, m_bye_mat);
            LBFGSBParam<double> param;
            LBFGSBSolver<double> solver(param);
            VectorXd lb = VectorXd::Constant(3, 1.0e-10);
            VectorXd ub = VectorXd::Constant(3, 1.0e10);
            // Initial values
            VectorXd x = varcom;
            double fx;
            int niter = solver.minimize(optimize_sigmaA, x, fx, lb, ub);
            varcom = x;
            optimize_sigmaA.calVi(varcom);
            Vi = optimize_sigmaA.get_Vi();
            V_logdet = optimize_sigmaA.get_V_logdet();
            
            for(long long i = 0; i < p; i++){
                X_xtVix(i) = X.col(i).transpose() * Vi * X.col(i);
            }
            shat2_Vec = 1 / X_xtVix.array();
            spdlog::info("Estimated variances: {}", join_string(varcom)) ;
        }
    }

    // alpha
    ofstream fout(this->m_out_file + ".alpha");
    if(!fout.is_open()) {
        spdlog::error("Fail to open {}.alpha", this->m_out_file);
        exit(1);
    }
    fout << alpha_Mat << std::endl;
    fout.close();

    // mu
    fout.open(this->m_out_file + ".mu");
    if(!fout.is_open()) {
        spdlog::error("Fail to open {}.mu", this->m_out_file);
        exit(1);
    }
    fout << mu_Mat << std::endl;
    fout.close();

    // get PIP
    Eigen::VectorXd pipVec = getPIP(alpha_Mat);
    fout.open(this->m_out_file + ".PIP");
    if(!fout.is_open()) {
        spdlog::error("Fail to open {}.PIP", this->m_out_file);
        exit(1);
    }
    for(auto val:pipVec){
        fout << val << std::endl;
    }
    fout.close();

    // get CS
    Eigen::MatrixXi status = in_CS(alpha_Mat, coverage);
    std::vector<std::vector<int>> cs_vec = getCS(status);
    Eigen::VectorXd claimed_coverage = computeClaimedCoverage(cs_vec, alpha_Mat);

    // purity according to min abs corr
    std::vector<std::vector<int>> cs_purity_vec = cs_vec;
    Eigen::VectorXd claimed_coverage_purity = claimed_coverage;
    getCSpurity(cs_purity_vec, claimed_coverage_purity, X, min_abs_corr);

    fout.open(this->m_out_file + ".claimed_coverage");
    if(!fout.is_open()) {
        spdlog::error("Fail to open {}.claimed_coverage", this->m_out_file);
        exit(1);
    }
    for(auto val:claimed_coverage_purity){
        fout << val << std::endl;
    }
    fout.close();

    fout.open(this->m_out_file + ".CS");
    if(!fout.is_open()) {
        spdlog::error("Fail to open {}.CS", this->m_out_file);
        exit(1);
    }
    for(auto vec:cs_purity_vec){
        for(auto x:vec){
            fout << x << " ";
        }
        fout << std::endl;
    }
    fout.close();
}



Eigen::VectorXi MMSUSIE::in_CS_x(const Eigen::VectorXd& x, double coverage) {

    // the number of elements for cumulative sum of coverage > 0.9
    Eigen::VectorXd sorted_x = x;
    std::sort(sorted_x.data(), sorted_x.data() + sorted_x.size(), std::greater<double>());

    double cumulative_sum = 0.0;
    int count = 0;
    for (int i = 0; i < sorted_x.size(); ++i) {
        cumulative_sum += sorted_x(i);
        count++;
        if (cumulative_sum >= coverage) break;
    }
    
    // set to 1 for elements to calcualte the cumulative sum
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < x.size(); ++i) {
        pairs.push_back({x(i), i});
    }

    std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    Eigen::VectorXi result(x.size());
    result.setZero();
    for (int i = 0; i < count; ++i) {
        result(pairs[i].second) = 1;
    }
    return result;
}


Eigen::MatrixXi MMSUSIE::in_CS(Eigen::MatrixXd& alpha, double coverage) {
    Eigen::MatrixXi status(alpha.rows(), alpha.cols());
    for (long long i = 0; i < alpha.rows(); i++) {
        Eigen::VectorXd x = alpha.row(i);
        status.row(i) = in_CS_x(x, coverage);
    }
    return status;
}


std::vector<std::vector<int>> MMSUSIE::getCS(const Eigen::MatrixXi& status) {
    std::vector<std::vector<int>> cs;
    
    for (int i = 0; i < status.rows(); ++i) {
        std::vector<int> currentRowIndices;
        for (int j = 0; j < status.cols(); ++j) {
            if (status(i, j) != 0) {
                currentRowIndices.push_back(j);
            }
        }
        cs.push_back(currentRowIndices);
    }

    return cs;
}


Eigen::VectorXd MMSUSIE::computeClaimedCoverage(const std::vector<std::vector<int>>& cs, const Eigen::MatrixXd& alpha) {
    Eigen::VectorXd claimed_coverage(cs.size());
    for(int i = 0; i < alpha.rows(); i++){
        auto currentSet = cs[i];
        double sum = 0.0;
        for (int index : currentSet) {
            sum += alpha(i, index);
        }
        claimed_coverage(i) = sum;
    }
    return claimed_coverage;
}


double MMSUSIE::computeMinCorrelation(const Eigen::MatrixXd& matrix) {
    // Number of rows (data points)
    int n = matrix.rows();
    
    // Mean center the data
    Eigen::MatrixXd centered = matrix.rowwise() - matrix.colwise().mean();
    
    // Compute the covariance matrix
    Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(n - 1);
    
    // Compute the correlation matrix from the covariance matrix
    Eigen::MatrixXd correlation = cov.array() / 
        (cov.diagonal().array().sqrt().replicate(1, cov.cols())).array() / 
        (cov.diagonal().array().sqrt().replicate(1, cov.cols())).transpose().array();
    
    // low tri elements
    long long num_row = correlation.rows();
    Eigen::VectorXd correlation_trilVec(num_row * (num_row - 1) / 2);
    long long k = 0;
    for(long long i = 0; i < num_row; i++){
        for(long long j = 0; j < i; j++){
            correlation_trilVec(k++) = std::fabs(correlation(i, j));
        }
    }
    return correlation_trilVec.minCoeff();
}


void MMSUSIE::getCSpurity(std::vector<std::vector<int>>& cs, Eigen::VectorXd& claimed_coverage, Eigen::MatrixXd& X, double& min_abs_corr){
    long long numCS = cs.size();
    std::vector<int> isPurity;
    for(long long i = 0; i < numCS; i++){
        std::vector<int> csi = cs[i];
        if(csi.size() == 1){
            isPurity.push_back(i);
        }else{
            Eigen::MatrixXd Xsub = X(Eigen::all, csi);
            double minCorr = computeMinCorrelation(Xsub);
            if(minCorr > min_abs_corr){
                isPurity.push_back(i);
            }
        }
    }

    std::vector<std::vector<int>> csPurity;
    Eigen::VectorXd claimed_coveragePurity(isPurity.size());
    for(long long i = 0; i < isPurity.size(); i++){
        csPurity.push_back(cs[isPurity[i]]);
        claimed_coveragePurity[i] = claimed_coverage[isPurity[i]];
    }
    cs = csPurity;
    claimed_coverage = claimed_coveragePurity;
}

Eigen::VectorXd MMSUSIE::getPIP(Eigen::MatrixXd& alpha_Mat){
    Eigen::MatrixXd alpha_Mat_tmp = 1 - alpha_Mat.array();
    Eigen::VectorXd pipVec(alpha_Mat_tmp.cols());
    for(long long j = 0; j < alpha_Mat_tmp.cols(); j++){
        double mac = 1;
        for(long long i = 0; i < alpha_Mat_tmp.rows(); i++){
            mac = mac * alpha_Mat_tmp(i, j);
        }
        pipVec(j) = 1 - mac;
    }
    return pipVec;
}


int MMSUSIE::run(int argc, char *argv[]) {
    CLI::App app{"mmSuSiE - Mixed model sum of single effect regression"};

    app.description(R"(
    Quick Start:
        fastgxe --mmsusie --grm grm_file --bfile bed_file --snp-focus snp_name --data data_file --trait BMI --env-int age smoking:alcohol --out test_mmsusie
    )");

    int threads = 10;
    int phen_correct = 2;
    bool standardize_env = true;
    
    std::string data_file, agrm_file, bed_file, out_file;
    std::vector<std::string> trait_vec, covariate_vec, class_vec, bye_vec, snp_focus_vec;
    std::vector<std::string> missing_in_data_vec = {"NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan", "<NA>", "<na>", "N/A", "n/a"};
    long long num_causal = 10;
    double tol = 1e-3;
    double coverage = 0.95;
    double min_abs_cor = 0.5;
    bool estimate_sigma = false;
    int maxiter = 100;
    double cc_par = 1.0e-7, cc_gra = 1.0e-6, cc_logL = 5.0e-5;

    // Register CLI11 options
    bool mmsusie=false;
    app.add_flag("--mmsusie", mmsusie, "Mixed model SuSiE");

    app.add_option("-p,--threads", threads, 
    "Number of threads to use (default: 10).")
    ->default_val(10);

    app.add_option("--data", data_file, 
        "Path to input data file (required).")
        ->required();

    app.add_option("--trait", trait_vec, 
        "Trait to analyze (required, expects 1 value).")
        ->expected(1)->required();
    
    app.add_option("--snp-focus", snp_focus_vec, 
        "SNP of focus.")
        ->expected(1)->required();
    
    app.add_option("--env-int", bye_vec, 
        "List of interacting environmental covariates (required).\n"
        "  - Supports multiple values (e.g., --env-int age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --env-int age:BMI).\n"
        "  - Expands to include all covariates in the specified range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1)->required();
    
    app.add_option("--out", out_file, 
        "Path to output file (required).")
        ->required();

    app.add_option("--grm", agrm_file, 
        "Path to genetic relationship matrix (GRM) file (required).")
        ->required();

    app.add_option("--bfile", bed_file, 
        "Path to binary PLINK BED file (optional).");

    app.add_option("--covar", covariate_vec, 
        "List of covariates for the analysis (required).\n"
        "  - Supports multiple values (e.g., --covar age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --covar age:BMI).\n"
        "  - Expands to include all covariates in the range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);

    app.add_option("--class", class_vec, 
        "List of categorical class variables (required).\n"
        "  - Supports multiple values (e.g., --class gender ethnicity region).\n"
        "  - Use ':' to specify a range (e.g., --class gender:region).\n"
        "  - Expands to include all variables in the range.\n"
        "  - Example order: {gender, ethnicity, region, education, income}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);

    app.add_flag("--no-standardize-env", [&standardize_env](int count) {
        if (count > 0) standardize_env = false;  // Flip standardize_env to false when flag is used
    }, "Disable standardization of interacting environmental covariates.\n"
       "  - Standardization is ENABLED by default (mean 0, std 1).\n"
       "  - Use --no-standardize-env to turn it OFF.");
    
    app.add_option("--phen-correct", phen_correct, "How to correct the phenotypes")
       ->check(CLI::Range(1, 2)); 
    
    app.add_option("--missing-data", missing_in_data_vec, 
        "List of missing value indicators for phenotype/covariates.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-data . -999 \"?\").")
        ->expected(-1);
    
    app.add_option("--num-causal", num_causal, "Number of causal variables")->check(CLI::PositiveNumber);
    app.add_option("--tol", tol, "Tolerance value for mmSuSiE")->check(CLI::PositiveNumber);
    app.add_option("--coverage", coverage, "Coverage value")->check(CLI::NonNegativeNumber);
    app.add_option("--min-abs-cor", min_abs_cor, "Minimum absolute correlation")->check(CLI::NonNegativeNumber);
    app.add_flag("--estimate-sigma", estimate_sigma, "Enable sigma estimation");
    
    app.add_option("--maxiter", maxiter, 
        "Maximum number of optimization iterations (default: 200).")
        ->default_val(200);

    app.add_option("--cc-par", cc_par, 
        "Convergence threshold for parameter updates (default: 1e-7).")
        ->default_val(1e-7);

    app.add_option("--cc-gra", cc_gra, 
        "Convergence threshold for gradient norm (default: 1e-6).")
        ->default_val(1e-6);

    app.add_option("--cc-logL", cc_logL, 
        "Convergence threshold for log-likelihood change (default: 5e-5).")
        ->default_val(5e-5);


    // Parse command-line arguments
    CLI11_PARSE(app, argc, argv);

    // Log parsed arguments
    spdlog::info("=== Parsed Arguments ===");
    spdlog::info("Threads: {}", threads);
    spdlog::info("Input Data File: {}", data_file);
    spdlog::info("Output File: {}", out_file);
    spdlog::info("GRM File: {}", agrm_file);
    spdlog::info("BED File: {}", bed_file.empty() ? "Not Provided" : bed_file);
    spdlog::info("Trait: {}", trait_vec[0]);
    spdlog::info("SNP of focus: {}", join_string(snp_focus_vec));
    spdlog::info("Interacting environmental covariates: {}", join_string(bye_vec, ", "));
    spdlog::info("Standardization of Interacting environments: {}", standardize_env ? "ENABLED" : "DISABLED");
    
    
    spdlog::info("Covariates: {}", covariate_vec.empty() ? "None" : join_string(covariate_vec, ", "));
    spdlog::info("Class Variables: {}", class_vec.empty() ? "None" : join_string(class_vec, ", "));
    spdlog::info("Missing Data Indicators: {}", join_string(missing_in_data_vec, ", "));

    spdlog::info("Max Iterations: {}", maxiter);
    spdlog::info("CC Par: {}, CC Gra: {}, CC LogL: {}", cc_par, cc_gra, cc_logL);

    spdlog::info("========================");
    
    // Set MKL and OpenMP threads
    if (threads <= 0) threads = 10;
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);

    // Obtain head row
    std::ifstream fin(data_file);
    if (!fin.is_open()) {
        spdlog::error("Failed to open data file: {}", data_file);
        exit(1);  // Ensure early exit to prevent further errors
    }

    std::string line;
    if (std::getline(fin, line)) {
        process_line(line);
        if (line.empty()) {
            spdlog::error("Head line is empty or starts with #");
            exit(1);
        }
    } else {
        spdlog::error("Failed to read the head line from file: {}", data_file);
        exit(1);
    }

    std::vector<std::string> head_vec = split_string(line);
    fin.close();

    // Trait index
    
    std::vector<std::string> strNoFound_vec;
    std::vector<long long> trait_index_vec = find_index(head_vec, trait_vec, strNoFound_vec);

    if (!strNoFound_vec.empty()) {
        spdlog::error("Trait names not found in the header: {}", join_string(strNoFound_vec));
        exit(1);
    }
    strNoFound_vec.clear();

    // covariate index
    covariate_vec = expand_variable_ranges(covariate_vec, head_vec);
    if(!covariate_vec.empty()){
        spdlog::info("Number of covariates: {}, Covariates: {}", 
             covariate_vec.size(), join_string(covariate_vec, ", "));
    }
    vector <long long> covariate_index_vec = find_index(head_vec, covariate_vec, strNoFound_vec);
    strNoFound_vec.clear();

    // class index
    class_vec = expand_variable_ranges(class_vec, head_vec);
    if(!class_vec.empty()){
        spdlog::info("Number of class variables: {}, class variables: {}", 
             class_vec.size(), join_string(class_vec, ", "));
    }
    vector <long long> class_index_vec = find_index(head_vec, class_vec, strNoFound_vec);
    strNoFound_vec.clear();


    // interacting environment
    bye_vec = expand_variable_ranges(bye_vec, head_vec);
    spdlog::info("Number of interacting environmental covariates: {}, Interacting environmental covariates: {}", 
             bye_vec.size(), join_string(bye_vec, ", "));
    vector <long long> bye_index_vec = find_index(head_vec, bye_vec, strNoFound_vec);
    if(strNoFound_vec.size() != 0){
        spdlog::error("Interacting environments not found in the header: {}", join_string(strNoFound_vec));
        exit(1);
    }
    strNoFound_vec.clear();
    
    if (!bye_index_vec.empty()) {
        std::vector<long long> tmp1 = find_index(covariate_vec, bye_vec, strNoFound_vec);
        std::vector<long long> tmp2 = find_index(class_vec, bye_vec, strNoFound_vec);
        
        if (!tmp1.empty() || !tmp2.empty()) {
            spdlog::error("Interacting environments should not be included in --covar or --class. They are treated as covariates by default.");
            exit(1);
        }
    }
    strNoFound_vec.clear();



    // prepare data
    this->pre_data(out_file, data_file, agrm_file, covariate_index_vec, class_index_vec, 
            bye_index_vec, trait_index_vec, missing_in_data_vec);
    
    VectorXd init_varcom;

    bool test_main = false;
    this->process_grm(test_main);
    Eigen::VectorXd varcom;
    phen_correct = true; // must be true, save resources
    this->pre_data_GxE(standardize_env, phen_correct);
    varcom = this->varcom_GxE(init_varcom, maxiter, cc_par, cc_gra, cc_logL);

    Eigen::MatrixXd X = this->pre_data_mmsusie(bed_file, snp_focus_vec, phen_correct);
    this->mmsusiefun(X, num_causal, maxiter, tol, coverage, min_abs_cor, estimate_sigma);
    
    return 0;
}

