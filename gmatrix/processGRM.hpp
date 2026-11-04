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

#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// GRM files hold float32 values and use the same layouts as fastGMA, so a GRM from
// either program works in both. Files from fastGxE 1.0 hold float64 and stay readable.
//   .grm.bin     lower triangle by row, as in GCTA; the width follows from the file size
//   .grm.sp.bin  the 8 bytes FGMASPF4, then records of (int64 row, int64 col, float32 value).
//                Without those bytes the values are float64. The number of records is not
//                known, so the file size cannot identify the width as it does for .grm.bin.

// Returns the bytes per value of a dense GRM, or stops if the size fits neither width.
int dense_grm_value_bytes(const std::string& bin_file, std::int64_t num_id);

// Leaves `fin` at the first record of a .grm.sp.bin and returns the bytes per value.
int read_sparse_grm_header(std::ifstream& fin);

// Reads one record of a .grm.sp.bin, returning false at a clean end of file. Stops on a
// partial record, an index outside [0, num_id) or a non-finite value, because a file
// read with the wrong layout would otherwise give wrong results without any error.
bool read_sparse_grm_record(std::ifstream& fin, int value_bytes, std::int64_t& index0, std::int64_t& index1,
                            double& value, const std::string& input_file, std::int64_t num_id);

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


    void merge_grm(const std::string& grm_file, int npart);
    void merge_file(const std::vector<std::string>& file_vec, const std::string& out_file, bool bin);
    

    void out_gmat(const Eigen::MatrixXd& gmat, const std::vector<std::string>& grm_id_vec, int out_fmt, const std::string& out_file);
    void out_gmat(const Eigen::SparseMatrix<double>& gmat, const std::vector<std::string>& grm_id_vec, int out_fmt, const std::string& out_file);

    void pca(const Eigen::MatrixXd& gmat, const std::string& out_file);

    void epistatic_grm(const std::string& prefix1, const std::string& prefix2, const std::string& out_file);

    void remove_close_id(const std::string& grm_file, int sparse, double cutoff, const std::string& out_file);

    void group_related_samples(const std::string& grm_file, double cutoff);

};
