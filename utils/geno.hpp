/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:13:31
 * @LastEditTime: 2025-01-30 13:47:51
 * @LastEditors: Chao Ning
 */
#ifndef GENO_HPP
#define GENO_HPP

#define EIGEN_USE_MKL_ALL  // Must be before including Eigen

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
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
    long long _num_id;
    long long _num_snp;

public:
    explicit GENO(string geno_file);

    long long get_num_iid();
    long long get_num_sid();
    vector<string> iid_vec();
    vector<string> sid_vec();
    vector<string> snp_anno(long long start_snp, long long num_snp_read);
    vector<string> snp_anno_by_snp_index(vector<long long> snp_index_vec);

    void bed_allele_freq(string in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr);
    void dgeno_allele_freq(string in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr, vector<string> missing_in_geno_vec);
    void allele_freq(int input_geno_fmt, VectorXd& freq_arr, VectorXd& nobs_geno_arr, vector<string> missing_in_geno_vec);
    
    vector<long long> find_fam_index(vector<string> id_in_need_vec);
    vector<long long> find_fam_index_pro(vector<string> id_in_need_vec);
    vector<long long> find_bim_index(vector<string> snp_in_need_vec);
    vector<long long> find_bim_index_pro(vector<string> snp_in_need_vec);

    void check_geno_id_order(string geno_file);
    void check_bed_numChar();

    void read_bed(string in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                  long long start_snp, long long num_snp_read, vector<long long> index_vec);

    void read_bed_omp_spMVLMM(string in_file, double* snp_mat_part_pt, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                              long long start_snp, long long num_snp_read, vector<long long> index_vec);

    void read_bed_omp(string in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                      long long start_snp, long long num_snp_read, vector<long long> index_vec);

    void read_dgeno(string in_file, MatrixXd& geno_mat_part, VectorXd& freq_arr, VectorXd& missing_rate_arr, 
                    long long start_snp, long long num_snp_read, vector<long long> index_vec, vector<string> missing_in_geno_vec);

    void read_geno(int input_geno_fmt, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                   long long start_snp, long long num_snp_read, vector<long long> index_vec, vector<string> missing_in_geno_vec);

    void read_bed_by_snp_index(string in_file, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                               vector<long long> snp_index_vec, vector<long long> index_vec);

    void read_dgeno_by_geno_index(string in_file, MatrixXd& geno_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                                  vector<long long> geno_index_vec, vector<long long> index_vec, vector<string> missing_in_geno_vec);

    void read_geno_by_index(int input_geno_fmt, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                            vector<long long> snp_index_vec, vector<long long> index_vec, vector<string> missing_in_geno_vec);

    void read_bed_2Dpart(MatrixXd& snp_mat_part, long long start_snp, long long num_snp_read, long long start_id, long long num_id_read);
};

#endif // GENO_HPP
