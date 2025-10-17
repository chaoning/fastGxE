/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:50
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-01 19:44:40
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


class MoM{


private:
      double m_maf_cut, m_missing_rate_cut;

public:
    MoM();
    ~MoM();
    void MoMV(bool no_noisebye, const string &out_file, const string &data_file, const vector<long long> &covariate_arr,
          const vector<long long> &class_arr, const vector<long long> &bye_arr, const vector<long long> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, long long start_snp, long long num_snp_read,
          long long num_randomB);
    void MoMMV(const string &out_file, const string &data_file, const vector<long long> &covariate_arr,
          const vector<long long> &class_arr, const vector<long long> &bye_arr, const vector<long long> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, long long start_snp, long long num_snp_read,
          long long num_randomB);
    int run(int argc, char* argv[]);

    void process_snps(long long num_snp_read, long long num_iid_bed, long long num_used_id,
                  long long num_randomB, const char* bytes_vec, const vector<long long> &index_vec,
                  MatrixXd &GB, MatrixXd &VB, const MatrixXd &B, const MatrixXd &bye_mat);
};
