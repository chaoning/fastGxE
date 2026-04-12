/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-08-12 17:23:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-31 12:05:33
 */
#include <cstdint>

#pragma once
#define EIGEN_USE_MKL_ALL  // !must be before include Eigen

#include <vector>
#include <string>
#include <map>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/SparseCore>

class ProcessGRM{

public:
    ProcessGRM();
    /**
     * @brief Destructor for Gmatrix class.
     */
    ~ProcessGRM();
    int run(int argc, char* argv[]);
    std::vector<std::string> read_grm_id(const std::string& grm_file);

    void out_grm_id(const std::string& out_file, const std::vector<std::string>& grm_id_vec);
    
    void read_grm_bin(const std::string& grm_file, const std::map<std::string, std::int64_t>& grm_id_map, Eigen::MatrixXd& mat, double val = 0);

    void read_grm_sp_bin(const std::string& grm_sparse_file, const std::map<std::string, std::int64_t>& grm_id_map, Eigen::SparseMatrix<double>& mat, double val = 0);


    void merge_grm(const std::string& grm_file, int npart, const std::string& out_file);
    void merge_file(const std::vector<std::string>& file_vec, const std::string& out_file, bool bin);
    

    void out_gmat(const Eigen::MatrixXd& gmat, const std::vector<std::string>& grm_id_vec, int out_fmt, const std::string& out_file);
    void out_gmat(const Eigen::SparseMatrix<double>& gmat, const std::vector<std::string>& grm_id_vec, int out_fmt, const std::string& out_file);

    void pca(const Eigen::MatrixXd& gmat, const std::string& out_file);

    void epistatic_grm(const std::string& prefix1, const std::string& prefix2, const std::string& out_file);

    void remove_close_id(const std::string& grm_file, int sparse, double cutoff, const std::string& out_file);

    void group_related_samples(const std::string& grm_file, double cutoff);

};
