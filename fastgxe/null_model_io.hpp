/*
 * @Description: Files that carry results between the steps of the GxE workflow
 * @Author: Chao Ning
 *
 *   <prefix>.null.bin   written by step 1, read by steps 2 and 3
 *   <prefix>.gamma      written by step 2, read by step 3
 */
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

// Everything steps 2 and 3 need from the fitted null model. All sample-level
// arrays follow `sample_ids`, the relatedness-group order set up in step 1.
struct NullModel {
    std::string trait_name;
    std::vector<std::string> sample_ids;
    std::vector<std::string> env_names;
    Eigen::VectorXd varcom;
    Eigen::VectorXd vinv_y;                 // V^{-1} y
    Eigen::MatrixXd env_mat;                // n x q, as used in the model
    Eigen::SparseMatrix<double> vinv;       // V^{-1}; left empty if not loaded
    bool no_noisebye = false;
    std::uint64_t fingerprint = 0;
};

// Calibration factors from the random SNP panel.
struct GxeGamma {
    double main = 0.0;
    Eigen::VectorXd no_main;                // one value per environment
    std::vector<Eigen::MatrixXd> with_main; // one 2 x 2 matrix per environment
    std::int64_t num_random_used = 0;
    std::int64_t num_random_requested = 0;
    std::uint64_t seed = 0;
    std::uint64_t null_fingerprint = 0;
};

// Identifies a null model, so that step 3 can reject a .gamma file that came from
// a different one. It depends only on the fitted model, not on when it was run.
std::uint64_t compute_null_fingerprint(const std::vector<std::string>& sample_ids,
                                       const std::vector<std::string>& env_names,
                                       const Eigen::VectorXd& varcom,
                                       const Eigen::VectorXd& vinv_y,
                                       const Eigen::MatrixXd& env_mat);

void write_null_model(const std::string& prefix, const NullModel& model);

// V^{-1} is the largest part of the file and is only needed by step 2 and the
// exact refit in step 3, so it can be skipped with load_vinv = false.
NullModel read_null_model(const std::string& prefix, bool load_vinv);

void write_gamma_file(const std::string& prefix, const GxeGamma& gamma,
                      const std::vector<std::string>& env_names);
GxeGamma read_gamma_file(const std::string& prefix, const std::vector<std::string>& env_names);
