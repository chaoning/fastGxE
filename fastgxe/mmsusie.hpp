/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2023-09-08 11:17:28
 * @LastEditTime: 2025-02-01 21:41:38
 * @LastEditors: Chao Ning
 */
#pragma once

#include "fastgxe.hpp"


class MMSUSIE:public fastGxE{
public:
        MMSUSIE();
        ~MMSUSIE();
        
        void pre_data_GxE(bool standardize_env, int phen_correct);
        Eigen::MatrixXd pre_data_mmsusie(string bed_file, vector<string> snp_vec, int phen_correct);

        void mmsusiefun(MatrixXd X, long long L, long long maxiter, double tol, 
                        double coverage, double min_abs_corr, bool estimate_sigma);
        void mmsusiefun2(MatrixXd X, VectorXd y, long long L, long long maxiter, double tol, 
                   double coverage, double min_abs_corr, bool estimate_sigma);

        void getCSpurity(std::vector<std::vector<int>>& cs, Eigen::VectorXd& claimed_coverage, Eigen::MatrixXd& X, double& min_abs_corr);
        double computeMinCorrelation(const Eigen::MatrixXd& matrix);
        Eigen::VectorXd computeClaimedCoverage(const std::vector<std::vector<int>>& cs, const Eigen::MatrixXd& alpha);
        std::vector<std::vector<int>> getCS(const Eigen::MatrixXi& status);
        Eigen::MatrixXi in_CS(Eigen::MatrixXd& res, double coverage = 0.9);
        Eigen::VectorXi in_CS_x(const Eigen::VectorXd& x, double coverage = 0.9);
        Eigen::VectorXd getPIP(Eigen::MatrixXd& alpha_Mat);

        int run(int argc, char* argv[]);

private:

};

