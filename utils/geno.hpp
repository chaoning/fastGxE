/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:13:31
 * @LastEditTime: 2025-01-30 13:47:51
 * @LastEditors: Chao Ning
 */
#include <cstdint>

#ifndef GENO_HPP
#define GENO_HPP

#define EIGEN_USE_MKL_ALL  // Must be before including Eigen

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::ifstream;
using std::vector;
using std::string;
using std::set;
using std::map;

class GENO {
private:
    string _geno_file;
    std::int64_t _num_id;
    std::int64_t _num_snp;
    vector<string> _iid_cache;
    vector<string> _sid_cache;
    bool _iid_cache_loaded = false;
    bool _sid_cache_loaded = false;

public:
    explicit GENO(const string& geno_file);

    std::int64_t get_num_iid();
    std::int64_t get_num_sid();
    vector<string> iid_vec();
    vector<string> sid_vec();
    vector<string> snp_anno(std::int64_t start_snp, std::int64_t num_snp_read);
    vector<string> snp_anno_by_snp_index(const vector<std::int64_t>& snp_index_vec);

    void bed_allele_freq(const string& in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr);
    void fs_feature_mean(const string& in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr, const vector<string>& missing_in_geno_vec);
    void allele_freq(int input_geno_fmt, VectorXd& freq_arr, VectorXd& nobs_geno_arr, const vector<string>& missing_in_geno_vec);
    
    vector<std::int64_t> find_fam_index(const vector<string>& id_in_need_vec);
    vector<std::int64_t> find_bim_index(const vector<string>& snp_in_need_vec);

    void validate_bed_size();

    void read_bed_centered_to_buffer_omp(const string& in_file, double* snp_mat_part_pt, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                                         std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec);

    void read_bed_omp(const string& in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                      std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec);

    void read_fs(const string& in_file, MatrixXd& geno_mat_part, VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                 std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec);

    void read_geno(int input_geno_fmt, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                   std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec);

    void read_bed_by_snp_indices(const string& in_file, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                                 const vector<std::int64_t>& snp_index_vec, const vector<std::int64_t>& index_vec);

    void read_fs_by_feature_indices(const string& in_file, MatrixXd& geno_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                                    const vector<std::int64_t>& feature_index_vec, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec);

    void read_geno_by_index(int input_geno_fmt, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                            const vector<std::int64_t>& snp_index_vec, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec);

};

#endif // GENO_HPP
