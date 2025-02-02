/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-08-12 17:23:46
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-31 12:05:33
 */
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
    std::vector<std::string> read_grm_id(std::string grm_file);
    std::vector<std::string> read_gcta_grm_id(std::string grm_file);

    void out_grm_id(std::string out_file, std::vector<std::string> grm_id_vec);
    
    void read_grm_bin(std::string grm_file, std::map<std::string, long long> grm_id_map, Eigen::MatrixXd& mat, double val = 0);
    void read_gcta_grm_bin(std::string grm_file, std::map<std::string, long long> grm_id_map, Eigen::MatrixXd& mat, double val = 0);

    void read_grm_sp_bin(std::string grm_sparse_file, std::map<std::string, long long> grm_id_map, Eigen::SparseMatrix<double>& mat, double val = 0);


    void merge_grm(std::string grm_file, int sparse, int npart, std::string out_file);
    void merge_file(std::vector<std::string> file_vec, std::string out_file, bool bin);
    

    void out_gmat(Eigen::MatrixXd& gmat, std::vector<std::string> grm_id_vec, int out_fmt, std::string out_file);
    void out_gmat(Eigen::SparseMatrix<double>& gmat, std::vector<std::string> grm_id_vec, int out_fmt, std::string out_file);

    void pca(Eigen::MatrixXd& gmat, std::string out_file);

    void epistatic_grm(std::string grm_file, int sparse, std::string grm_type, std::string out_file);

    void remove_id(std::string grm_file, int sparse, std::string id_file, std::string out_file);

    void keep_id(std::string grm_file, int sparse, std::string id_file, std::string out_file);

    void remove_close_id(std::string grm_file, int sparse, double cutoff, std::string out_file);

    void group(std::string grm_file, int sparse, double cutoff, std::string out_file);

};

