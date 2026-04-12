/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:50
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-01 20:32:22
 */
#include <cstdint>

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Sparse>


using std::string;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

class PHEN;

class fastGxE{

protected:
    // Output prefix and sample-level bookkeeping shared across the analysis flow.
    string m_out_file;

    std::int64_t m_num_trait;
    std::int64_t m_num_id;

    vector<string>  m_id_in_data_vec;

    // Phenotype, fixed-effect design and interacting environments.
    MatrixXd m_y, m_xmat, m_bye_mat;
    std::int64_t m_num_bye;
    VectorXd m_enviW_Vec;
    MatrixXd m_y_trans, m_xmat_trans;
    VectorXd m_varcom_null;

    // GRM blocks and the sparse/full representations derived from them.
    vector<MatrixXd> m_grm_mat_group_vec;
    vector<std::int64_t> m_grm_index_vec;
    SparseMatrix < double > m_grm_mat;

    // Null-model covariance caches reused by downstream SNP scans.
    SparseMatrix < double > m_Vi0;
    SparseMatrix < double > m_grm_eigenvecs;
    VectorXd m_grm_eigenvals;
    double m_logL_null, m_V0_logdet;

    vector<std::int64_t> m_condition_snp_index_vec;

private:
    struct GroupRecord {
        std::int64_t grm_index = 0;
        std::int64_t group_id = 0;
        std::int64_t group_size = 0;
    };

    struct SampleBlockLocation {
        std::int64_t block_id = 0;
        std::int64_t row_index = 0;
    };

    struct GroupLayout {
        vector<MatrixXd> grm_blocks;
        std::unordered_map<std::int64_t, SampleBlockLocation> sample_location_map;
    };

    // Data preparation helpers: align phenotype rows with the GRM order, then
    // build the fixed-effect and environment matrices once for all later steps.
    void initialize_analysis_metadata(const string& out_file,
            const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
            const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait);
    void load_filtered_phenotype(PHEN& phenoA, const string& data_file,
            const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
            const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait,
            const vector<string>& grm_sample_ids, const vector<string>& missing_in_data_vec);
    vector<string> extract_retained_sample_ids(const PHEN& phenoA) const;
    vector<GroupRecord> load_and_sort_group_records(const string& agrm_file,
            const vector<string>& retained_sample_ids, std::int64_t expected_num_samples) const;
    vector<string> build_reordered_sample_ids(const vector<GroupRecord>& group_records,
            const vector<string>& grm_sample_ids) const;
    vector<GroupRecord> reorder_phenotype_by_group(PHEN& phenoA, const string& agrm_file,
            const vector<string>& grm_sample_ids);
    void log_removed_design_columns(const vector<std::int64_t>& removed_column_indices) const;
    void build_fixed_and_environment_matrices(PHEN& phenoA,
            const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
            const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait);
    GroupLayout build_group_layout(const vector<GroupRecord>& group_records) const;
    void load_grouped_grm_blocks(const string& agrm_file,
            const std::unordered_map<std::int64_t, SampleBlockLocation>& sample_location_map,
            vector<MatrixXd>& grm_blocks) const;

public:
    fastGxE();
    ~fastGxE();

    // Top-level pipeline orchestration.
    int run(int argc, char **argv);
    int pre_data(const string& out_file, const string& data_file, const string& agrm_file,
            const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
            const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait,
            const vector<string>& missing_in_data_vec);
    void process_grm(bool use_eigen);
    void pre_data_GxE(bool standardize_env, bool phen_correct);

    // Utility helpers reused by different scan modes.
    void reset_output_prefix(const string& out_file, std::int64_t first, std::int64_t second);

    bool isValidSnp(double af, double missing_rate, double maf_cut, double missing_rate_cut);

    // Null-model fitting and genome-wide scans for SNP main effects.
    VectorXd varcom_main(VectorXd& init_varcom, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0);
    void test_main(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut, 
                int maxiter, double cc_par, double cc_gra);
    void output(std::ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
                const double* eff_arr, const double* se_arr, const VectorXd& p_arr,
                std::int64_t start_snp, std::int64_t num_snp_read, double p_cut);

    // Null-model fitting and genome-wide scans for SNP-by-environment effects.
    VectorXd varcom_GxE(VectorXd& init_varcom, bool no_noisebye, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0);

    void test_GxE_multi(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut);
    void output_GxE_multi(std::ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
            const VectorXd& beta_main_Vec, const VectorXd& se_main_Vec, const VectorXd& p_main_Vec,
            const VectorXd& score_Vec, const VectorXd& p_Vec, std::int64_t start_snp, std::int64_t num_snp_read);

    void test_GxE(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut);
    void output_GxE(std::ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
            const MatrixXd& beta_mat, const MatrixXd& se_mat, const MatrixXd& p_mat, const MatrixXd& p_combined_Mat,
            std::int64_t start_snp, std::int64_t num_snp_read);

    vector<std::int64_t> get_start_end_pos(const string& bed_file,
            const vector<string>& snp_range_vec, const vector<std::int64_t>& split_task_vec);
};
