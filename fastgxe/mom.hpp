/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:50
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-01 19:44:40
 */
#include <cstdint>

#pragma once

#include <string>
#include <vector>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using std::string;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;


class MoM{


private:
      double m_maf_cut, m_missing_rate_cut;

public:
    MoM();
    ~MoM();
    void MoMV(bool no_noisebye, const string &out_file, const string &data_file, const vector<std::int64_t> &covariate_arr,
          const vector<std::int64_t> &class_arr, const vector<std::int64_t> &bye_arr, const vector<std::int64_t> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, std::int64_t start_snp, std::int64_t num_snp_read,
          std::int64_t num_randomB);
    void MoMMV(const string &out_file, const string &data_file, const vector<std::int64_t> &covariate_arr,
          const vector<std::int64_t> &class_arr, const vector<std::int64_t> &bye_arr, const vector<std::int64_t> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, std::int64_t start_snp, std::int64_t num_snp_read,
          std::int64_t num_randomB);
    int run(int argc, char* argv[]);

    void process_snps(std::int64_t num_snp_read, std::int64_t num_iid_bed, std::int64_t num_used_id,
                  std::int64_t num_randomB, const char* bytes_vec, const vector<std::int64_t> &index_vec,
                  MatrixXd &GB, MatrixXd &VB, const MatrixXd &B, const MatrixXd &bye_mat);
};
