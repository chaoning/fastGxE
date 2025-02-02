/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:50
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-01 20:32:22
 */
#pragma once

#include <string>
#include <vector>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

using std::string;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;


class fastGxE{


protected:
    string m_out_file;

    long long m_num_trait;
    long long m_num_id;

    vector<string>  m_id_in_data_vec;

    MatrixXd m_y, m_xmat, m_bye_mat;
    long long m_num_bye;
    VectorXd m_enviW_Vec;
    MatrixXd m_y_trans, m_xmat_trans;
    VectorXd m_varcom_null;

    vector<MatrixXd> m_grm_mat_group_vec;
    vector<long long> m_grm_index_vec;
    SparseMatrix < double > m_grm_mat;

    SparseMatrix < double > m_Vi0;
    SparseMatrix < double > m_grm_eigenvecs;
    VectorXd m_grm_eigenvals;
    double m_logL_null, m_V0_logdet;

    vector<long long> m_condition_snp_index_vec;


public:
    fastGxE();
    ~fastGxE();
    int run(int argc, char **argv);
    int pre_data(string out_file, string data_file, string agrm_file, vector<long long> covariate_arr, 
            vector<long long> class_arr, vector<long long> bye_arr, vector<long long> trait, 
            vector<string> missing_in_data_vec);
    void process_grm(bool use_eigen);
    void pre_data_GxE(bool standardize_env, bool phen_correct);

    void reset_output_prefix(string out_file, long long first, long long second);

    bool isValidSnp(double af, double missing_rate, double maf_cut, double missing_rate_cut);

    // main
    VectorXd varcom_main(VectorXd& init_varcom, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0);
    void test_main(string bed_file, 
                long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut, 
                int maxiter, double cc_par, double cc_gra);
    void output(std::ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr, 
                double* eff_arr, double* se_arr, VectorXd& p_arr,
                long long start_snp, long long num_snp_read, double p_cut);

    // GxE
    VectorXd varcom_GxE(VectorXd& init_varcom, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0);

    void test_GxE_multi(string bed_file, long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut);
    void output_GxE_multi(std::ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr,
            VectorXd& beta_main_Vec, VectorXd& se_main_Vec, VectorXd& p_main_Vec, 
            VectorXd& score_Vec, VectorXd& p_Vec, long long start_snp, long long num_snp_read);

    void test_GxE(string bed_file, 
                long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut);
    void output_GxE(std::ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr,
            MatrixXd& beta_mat, MatrixXd& se_mat, MatrixXd& p_mat, MatrixXd& p_combined_Mat,
            long long start_snp, long long num_snp_read);

    vector<long long> get_start_end_pos(string bed_file, 
            vector<string> snp_range_vec, vector<long long> split_task_vec);
};
