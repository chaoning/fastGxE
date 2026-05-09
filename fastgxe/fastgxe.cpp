/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:40
 * LastEditors: Chao Ning
 * LastEditTime: 2026-04-11 20:42:12
 */

#include <cstdint>

#define EIGEN_USE_MKL_ALL  // must be before include Eigen

#include <getopt.h>
#include "mkl.h"
#include "omp.h"

#include <limits>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cctype>
#include <string>
#include <utility>
#include <random>
#include <fstream>
#include <chrono>
#include <gsl/gsl_cdf.h>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>  // Include CLI11


#include "../utils/geno.hpp"
#include "../gmatrix/processGRM.hpp"
#include "../utils/phen.hpp"
#include "../utils/EigenMatrix_utils.hpp"
#include "../utils/chi2score.hpp"
#include "../utils/iterator_utils.hpp"
#include "../utils/string_utils.hpp"
#include "../utils/someCommonFun.hpp"
#include "../utils/acat_utils.hpp"
#include "../utils/random_sampling.hpp"
#include "../utils/fatal_error.hpp"
#include "fastgxe.hpp"


using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::endl;
using Eigen::HouseholderQR;
using Eigen::ColPivHouseholderQR;

namespace {

// Numerical constants for internal algorithms.
// These are intentionally not CLI-configurable: changing them would require
// re-validation of statistical calibration.
constexpr double kMinVarianceComponent   = 1e-6;  // smallest allowed variance component
constexpr double kBinaryTraitTolerance   = 1e-8;  // 0/1 phenotype rounding threshold
constexpr double kAiEmFloor              = 1e-12; // minimum AI/EM information diagonal
constexpr double kEmStabilityNudge       = 1e-10; // regularizer added to EM info blocks
constexpr double kCorrMatrixNudge        = 0.001; // ridge added to environment correlation matrix
constexpr double kProbabilityWeightFloor = 1e-6;  // minimum logistic weight to avoid division by zero

std::vector<std::string> read_header_columns(const std::string& data_file) {
    std::ifstream fin(data_file);
    if (!fin.is_open()) {
        fatal_error("Failed to open data file: {}", data_file);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        fatal_error("Failed to read the head line from file: {}", data_file);
    }

    process_line(line);
    if (line.empty()) {
        fatal_error("Head line is empty or starts with #");
    }

    return split_string(line);
}

std::string join_or_none(const std::vector<std::string>& values,
                         const char* delimiter = ", ") {
    return values.empty() ? "None" : join_string(values, delimiter);
}

std::string sanitize_output_label(const std::string& raw_label) {
    std::string sanitized = raw_label;
    for (char& ch : sanitized) {
        const unsigned char uch = static_cast<unsigned char>(ch);
        if (!(std::isalnum(uch) || ch == '_')) {
            ch = '_';
        }
    }
    return sanitized.empty() ? "trait" : sanitized;
}

struct RunOptions {
    int threads = 10;
    bool test_main = false;
    bool test_main_binary = false;
    bool test_main_binary_continuous = false;
    bool test_main_multitrait_continuous = false;
    bool test_gxe = false;
    bool standardize_env = true;
    bool no_noisebye = false;

    std::string data_file;
    std::string agrm_file;
    std::string bed_file;
    std::string out_file;

    std::vector<std::string> trait_vec;
    std::vector<std::string> covariate_vec;
    std::vector<std::string> class_vec;
    std::vector<std::string> bye_vec;
    std::vector<std::string> snp_range_vec;
    std::vector<std::string> missing_in_data_vec = {
        "NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan", "<NA>", "<na>", "N/A", "n/a"
    };
    std::vector<std::int64_t> split_task_vec;

    std::int64_t num_random_snp = 2000;
    int npart_snp = 20;
    int speed = 2;
    int maxiter = 100;
    double p_approx_cut = 1.0e-3;
    double p_cut = 2;
    double maf_cut = 0.01;
    double missing_rate_cut = 0.05;
    double cc_par = 1.0e-7;
    double cc_gra = 1.0e-6;
    double cc_logL = 5.0e-5;
};

struct ResolvedInputColumns {
    std::vector<std::int64_t> trait_index_vec;
    std::vector<std::int64_t> covariate_index_vec;
    std::vector<std::int64_t> class_index_vec;
    std::vector<std::int64_t> bye_index_vec;
};

struct GxeRelationshipData {
    std::vector<MatrixXd> group_blocks;
    SparseMatrix<double> sparse_matrix;
};

struct InverseCovarianceResult {
    SparseMatrix<double> inverse_matrix;
    double logdet = 0.0;
};

struct BinaryWorkingLmmState {
    VectorXd mu;
    VectorXd weight;
    VectorXd residual_diag;
    VectorXd working_response;
};

struct ProjectionCache {
    MatrixXd inverse_times_design;
    MatrixXd design_cross_inverse;
    double design_logdet = 0.0;
};

struct ScalarProjectionCache {
    VectorXd inverse_times_design;
    double design_cross_inverse = 0.0;
};

struct BinaryProjectionState {
    InverseCovarianceResult inverse_result;
    ProjectionCache fixed_effect_projection;
    VectorXd projected_response;
    VectorXd genetic_times_projected_response;
    VectorXd projected_genetic_response;
    double trace_pg = 0.0;
};

struct ScalarVarianceUpdateState {
    double score = 0.0;
    double ai = 0.0;
    double em = 0.0;
    double gamma = 1.0;
    double blended_information = 0.0;
    double delta = 0.0;
    double updated_variance = 0.0;
};

struct BinaryVarianceFitState {
    double genetic_var = 0.0;
    double neg2_reml = 0.0;
    double score = 0.0;
    double ai = 0.0;
    double em = 0.0;
    double gamma = 1.0;
    bool converged = false;
    BinaryProjectionState projection_state;
};

struct JointMainWorkingState {
    VectorXd binary_mu;
    VectorXd binary_weight;
    VectorXd binary_residual_diag;
    VectorXd binary_working_response;
    VectorXd stacked_response;
};

struct JointInverseCovarianceResult {
    SparseMatrix<double> inverse_matrix;
    VectorXd trace_vi_components = VectorXd::Zero(4);
    double logdet = 0.0;
};

struct JointProjectionState {
    JointInverseCovarianceResult inverse_result;
    ProjectionCache fixed_effect_projection;
    VectorXd projected_response;
    VectorXd projected_random_response;
    double neg2_reml = 0.0;
};

struct JointVarianceFitState {
    VectorXd varcom = VectorXd::Zero(4);
    VectorXd score = VectorXd::Zero(4);
    MatrixXd ai = MatrixXd::Zero(4, 4);
    MatrixXd em = MatrixXd::Zero(4, 4);
    double neg2_reml = 0.0;
    bool converged = false;
    JointProjectionState projection_state;
};

struct JointCovarianceFitState {
    VectorXd varcom = VectorXd::Zero(4);
    double score = 0.0;
    double ai = 0.0;
    double em = 0.0;
    double gamma = 1.0;
    double neg2_reml = 0.0;
    bool converged = false;
    JointProjectionState projection_state;
};

struct MultitraitContinuousInverseCovarianceResult {
    SparseMatrix<double> inverse_matrix;
    VectorXd trace_vi_components;
    double logdet = 0.0;
};

struct MultitraitContinuousProjectionState {
    MultitraitContinuousInverseCovarianceResult inverse_result;
    ProjectionCache fixed_effect_projection;
    VectorXd projected_response;
    double neg2_reml = 0.0;
};

struct MultitraitContinuousVarianceFitState {
    VectorXd varcom;
    VectorXd score;
    MatrixXd ai;
    MatrixXd em;
    double neg2_reml = 0.0;
    bool converged = false;
    MultitraitContinuousProjectionState projection_state;
};

struct VarianceUpdateResult {
    VectorXd delta;
    VectorXd updated_varcom;
};

struct MainModelIterationState {
    double logL = 0.0;
    VectorXd fd_mat = VectorXd::Zero(2);
    MatrixXd ai_mat = MatrixXd::Zero(2, 2);
};

struct BedScanContext {
    std::vector<std::int64_t> id_index_in_bed_vec;
    std::int64_t num_id_used = 0;
    std::int64_t num_snp = 0;
    std::int64_t num_id_in_bed = 0;
};

struct SnpPartition {
    std::int64_t start_snp = 0;
    std::int64_t num_snp_read = 0;
};

// Work buffers for scalar (single-trait) SNP scans. Shared by test_main and
// test_main_binary to avoid duplicating the same resize logic.
struct SnpScalarScanBuffers {
    std::vector<double> snp_mat;
    std::vector<double> geno_var;
    std::vector<double> se;
    std::vector<double> prod;
    std::vector<double> eff;
    VectorXd p;

    void resize(std::int64_t num_snp, std::int64_t num_id) {
        const std::size_t n = static_cast<std::size_t>(num_snp);
        snp_mat.resize(n * static_cast<std::size_t>(num_id));
        geno_var.resize(n);
        se.resize(n);
        prod.resize(n);
        eff.resize(n);
        // NaN marks SNPs not yet tested; output() skips rows where p is still NaN.
        p = VectorXd::Constant(num_snp, std::numeric_limits<double>::quiet_NaN());
    }

    bool needs_resize(std::int64_t num_snp) const {
        return num_snp != static_cast<std::int64_t>(p.size());
    }
};

struct SnpMatrixScanBuffers {
    std::vector<double> snp_mat;
    MatrixXd beta_mat;
    MatrixXd se_mat;
    MatrixXd p_mat;

    void resize(std::int64_t num_snp, std::int64_t num_id,
                std::int64_t num_eff_cols, std::int64_t num_p_cols) {
        snp_mat.resize(static_cast<std::size_t>(num_snp) * static_cast<std::size_t>(num_id));
        const double nan = std::numeric_limits<double>::quiet_NaN();
        beta_mat = MatrixXd::Constant(num_snp, num_eff_cols, nan);
        se_mat   = MatrixXd::Constant(num_snp, num_eff_cols, nan);
        p_mat    = MatrixXd::Constant(num_snp, num_p_cols,   nan);
    }

    bool needs_resize(std::int64_t num_snp) const {
        return num_snp != beta_mat.rows();
    }

    void reset_values() {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        beta_mat.setConstant(nan);
        se_mat.setConstant(nan);
        p_mat.setConstant(nan);
    }
};

using GrmTriplet = Eigen::Triplet<double>;

bool is_binary_trait_matrix(const MatrixXd& y, double tolerance);

std::vector<std::int64_t> build_grm_block_offsets(const std::vector<MatrixXd>& grm_blocks) {
    std::vector<std::int64_t> offsets;
    offsets.reserve(grm_blocks.size() + 1);
    offsets.push_back(0);

    for (const auto& block : grm_blocks) {
        offsets.push_back(offsets.back() + block.rows());
    }

    return offsets;
}

std::size_t count_grm_triplets(const std::vector<MatrixXd>& grm_blocks) {
    std::size_t triplet_count = 0;
    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const std::int64_t block_size = grm_blocks[block_idx].rows();
        // Block 0 holds unrelated (singleton) samples stored as a diagonal only, not a full dense block.
        if (block_idx == 0) {
            triplet_count += static_cast<std::size_t>(block_size);
        } else {
            triplet_count += static_cast<std::size_t>(block_size * block_size);
        }
    }
    return triplet_count;
}

void append_singleton_block_triplets(const MatrixXd& singleton_block, std::int64_t start_index,
        std::vector<GrmTriplet>& triplets, bool use_identity) {
    const std::int64_t block_size = singleton_block.rows();
    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        const double value = use_identity ? 1.0 : singleton_block(row_idx, 0);
        triplets.emplace_back(start_index + row_idx, start_index + row_idx, value);
    }
}

void append_dense_block_triplets(const MatrixXd& dense_block, std::int64_t start_index,
        std::vector<GrmTriplet>& triplets) {
    const std::int64_t block_size = dense_block.rows();
    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        for (std::int64_t col_idx = 0; col_idx < block_size; ++col_idx) {
            triplets.emplace_back(
                start_index + row_idx,
                start_index + col_idx,
                dense_block(row_idx, col_idx));
        }
    }
}

void load_singleton_block_eigenvalues(const MatrixXd& singleton_block, std::int64_t start_index,
        VectorXd& eigenvalues) {
    const std::int64_t block_size = singleton_block.rows();
    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        eigenvalues(start_index + row_idx) = singleton_block(row_idx, 0);
    }
}

void append_dense_block_eigendecomposition(const MatrixXd& dense_block, std::int64_t start_index,
        VectorXd& eigenvalues, std::vector<GrmTriplet>& triplets) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(dense_block);
    if (eigensolver.info() != Eigen::Success) {
        fatal_error("Eigen decomposition failed for GRM block starting at index {}.", start_index);
    }

    const VectorXd& block_eigenvalues = eigensolver.eigenvalues();
    const MatrixXd& block_eigenvectors = eigensolver.eigenvectors();
    const std::int64_t block_size = dense_block.rows();

    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        eigenvalues(start_index + row_idx) = block_eigenvalues(row_idx);
        for (std::int64_t col_idx = 0; col_idx < block_size; ++col_idx) {
            triplets.emplace_back(
                start_index + row_idx,
                start_index + col_idx,
                block_eigenvectors(row_idx, col_idx));
        }
    }
}

ColPivHouseholderQR<MatrixXd> require_full_column_rank(
        const MatrixXd& design_matrix,
        std::int64_t expected_rank,
        const char* error_message) {
    ColPivHouseholderQR<MatrixXd> qr(design_matrix);
    if (qr.rank() < expected_rank) {
        fatal_error("{}", error_message);
    }
    return qr;
}

void standardize_environment_matrix(MatrixXd& environment_matrix) {
    const VectorXd mean_vector = environment_matrix.colwise().mean();
    MatrixXd centered_matrix = environment_matrix.rowwise() - mean_vector.transpose();
    const VectorXd std_vector =
        (centered_matrix.cwiseProduct(centered_matrix).colwise().sum().array() / centered_matrix.rows()).sqrt();
    environment_matrix = centered_matrix.array().rowwise() * (1.0 / std_vector.transpose().array());
}

void append_columns(MatrixXd& design_matrix, const MatrixXd& extra_columns) {
    const std::int64_t original_col_count = design_matrix.cols();
    design_matrix.conservativeResize(design_matrix.rows(), original_col_count + extra_columns.cols());
    design_matrix.rightCols(extra_columns.cols()) = extra_columns;
}

std::string format_vector_for_log(const VectorXd& values) {
    const Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "", "", "", "", "");
    std::ostringstream oss;
    oss << values.transpose().format(fmt);
    return oss.str();
}

VectorXd initialize_variance_components(const VectorXd& initial_values, std::int64_t expected_size) {
    if (initial_values.size() == expected_size && initial_values.minCoeff() > 0) {
        return initial_values;
    }
    return VectorXd::Ones(expected_size);
}

std::int64_t num_lower_triangle_elements(std::int64_t matrix_size) {
    return matrix_size * (matrix_size + 1) / 2;
}

std::pair<std::int64_t, std::int64_t> lower_triangle_pair_from_index(
        std::int64_t packed_index, std::int64_t matrix_size) {
    std::int64_t current_index = 0;
    for (std::int64_t row = 0; row < matrix_size; ++row) {
        for (std::int64_t col = 0; col <= row; ++col) {
            if (current_index == packed_index) {
                return {row, col};
            }
            ++current_index;
        }
    }

    fatal_error("Packed covariance index {} is out of range for matrix size {}.",
        packed_index, matrix_size);
}

MatrixXd unpack_lower_triangle_to_symmetric_matrix(
        const VectorXd& packed_values, std::int64_t matrix_size, std::int64_t offset = 0) {
    MatrixXd covariance_matrix = MatrixXd::Zero(matrix_size, matrix_size);
    std::int64_t packed_index = offset;
    for (std::int64_t row = 0; row < matrix_size; ++row) {
        for (std::int64_t col = 0; col <= row; ++col) {
            covariance_matrix(row, col) = packed_values(packed_index);
            covariance_matrix(col, row) = packed_values(packed_index);
            ++packed_index;
        }
    }
    return covariance_matrix;
}

VectorXd pack_symmetric_matrix_lower_triangle(const MatrixXd& symmetric_matrix) {
    const std::int64_t matrix_size = symmetric_matrix.rows();
    VectorXd packed_values(num_lower_triangle_elements(matrix_size));
    std::int64_t packed_index = 0;
    for (std::int64_t row = 0; row < matrix_size; ++row) {
        for (std::int64_t col = 0; col <= row; ++col) {
            packed_values(packed_index) = symmetric_matrix(row, col);
            ++packed_index;
        }
    }
    return packed_values;
}

MatrixXd build_repeated_fixed_effect_design(const MatrixXd& xmat, std::int64_t num_traits) {
    MatrixXd repeated_design = MatrixXd::Zero(xmat.rows() * num_traits, xmat.cols() * num_traits);
    for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
        repeated_design.block(
            trait_index * xmat.rows(),
            trait_index * xmat.cols(),
            xmat.rows(),
            xmat.cols()) = xmat;
    }
    return repeated_design;
}

VectorXd stack_trait_matrix_rows(const MatrixXd& trait_matrix) {
    VectorXd stacked_response(trait_matrix.rows() * trait_matrix.cols());
    for (std::int64_t trait_index = 0; trait_index < trait_matrix.cols(); ++trait_index) {
        stacked_response.segment(
            trait_index * trait_matrix.rows(),
            trait_matrix.rows()) = trait_matrix.col(trait_index);
    }
    return stacked_response;
}

VectorXd fit_linear_fixed_effects(const MatrixXd& xmat, const VectorXd& y) {
    const ColPivHouseholderQR<MatrixXd> qr = require_full_column_rank(
        xmat,
        xmat.cols(),
        "Failed to initialize the continuous-trait fixed-effect model.");
    return qr.solve(y);
}

bool is_valid_joint_main_variance_components(const VectorXd& varcom,
        double min_variance = 1e-6, double min_det = 1e-8) {
    if (varcom.size() != 4) {
        return false;
    }
    if (varcom(0) <= min_variance || varcom(2) <= min_variance || varcom(3) <= min_variance) {
        return false;
    }
    const double det_sigma_g = varcom(0) * varcom(2) - varcom(1) * varcom(1);
    return det_sigma_g > min_det;
}

bool is_valid_multitrait_continuous_variance_components(const VectorXd& varcom,
        std::int64_t num_traits,
        double min_variance = 1e-6,
        double min_eigenvalue = 1e-8) {
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    if (varcom.size() != 2 * num_covariance_terms) {
        return false;
    }

    const MatrixXd genetic_covariance = unpack_lower_triangle_to_symmetric_matrix(
        varcom, num_traits, 0);
    const MatrixXd residual_covariance = unpack_lower_triangle_to_symmetric_matrix(
        varcom, num_traits, num_covariance_terms);

    if ((genetic_covariance.diagonal().array() <= min_variance).any() ||
            (residual_covariance.diagonal().array() <= min_variance).any()) {
        return false;
    }

    Eigen::LDLT<MatrixXd> genetic_ldlt(genetic_covariance);
    Eigen::LDLT<MatrixXd> residual_ldlt(residual_covariance);
    if (genetic_ldlt.info() != Eigen::Success || residual_ldlt.info() != Eigen::Success) {
        return false;
    }

    return genetic_ldlt.vectorD().minCoeff() > min_eigenvalue &&
           residual_ldlt.vectorD().minCoeff() > min_eigenvalue;
}

void require_binary_continuous_trait_matrix(const MatrixXd& y) {
    if (y.cols() != 2 || y.rows() == 0) {
        fatal_error("The binary-continuous joint model requires exactly two traits: one binary and one continuous.");
    }
    if (!is_binary_trait_matrix(MatrixXd(y.leftCols(1)), 1e-8)) {
        fatal_error("The first trait in the binary-continuous joint model must be coded as 0/1.");
    }
    const VectorXd y_binary = y.col(0);
    if (y_binary.mean() <= 0.0 || y_binary.mean() >= 1.0) {
        fatal_error("The binary trait in the binary-continuous joint model requires both cases and controls.");
    }
    const VectorXd y_continuous = y.col(1);
    const double centered_norm =
        (y_continuous.array() - y_continuous.mean()).matrix().squaredNorm();
    if (!(centered_norm > 0.0)) {
        fatal_error("The second trait in the binary-continuous joint model must vary across samples.");
    }
}

void require_multitrait_continuous_trait_matrix(const MatrixXd& y) {
    if (y.cols() < 2 || y.rows() == 0) {
        fatal_error("The multitrait continuous model requires at least two continuous traits.");
    }

    for (int trait_index = 0; trait_index < y.cols(); ++trait_index) {
        const VectorXd trait = y.col(trait_index);
        const double centered_norm =
            (trait.array() - trait.mean()).matrix().squaredNorm();
        if (!(centered_norm > 0.0)) {
            fatal_error("Trait {} in the multitrait continuous model must vary across samples.",
                trait_index + 1);
        }
    }
}

MatrixXd build_gxe_relationship_block(const MatrixXd& grm_block, const MatrixXd& bye_block, std::int64_t num_bye,
        bool is_singleton_block) {
    if (is_singleton_block) {
        MatrixXd gxe_block = (bye_block.cwiseProduct(bye_block)).rowwise().sum() / num_bye;
        return gxe_block.cwiseProduct(grm_block);
    }

    MatrixXd gxe_block = bye_block * bye_block.transpose() / num_bye;
    return gxe_block.cwiseProduct(grm_block);
}

GxeRelationshipData build_gxe_relationship_data(const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets, const MatrixXd& bye_mat,
        std::int64_t num_bye, std::int64_t num_id) {
    GxeRelationshipData relationship_data;
    relationship_data.group_blocks.resize(grm_blocks.size());

    std::vector<GrmTriplet> triplets;
    triplets.reserve(count_grm_triplets(grm_blocks));

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (grm_block.size() == 0) {
            relationship_data.group_blocks[block_idx] = grm_block;
            continue;
        }

        const MatrixXd bye_block = bye_mat.middleRows(start_index, block_size);
        const MatrixXd gxe_block = build_gxe_relationship_block(
            grm_block, bye_block, num_bye, block_idx == 0);
        relationship_data.group_blocks[block_idx] = gxe_block;

        if (block_idx == 0) {
            append_singleton_block_triplets(gxe_block, start_index, triplets, false);
        } else {
            append_dense_block_triplets(gxe_block, start_index, triplets);
        }
    }

    relationship_data.sparse_matrix.resize(num_id, num_id);
    relationship_data.sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return relationship_data;
}

VectorXd build_noise_by_environment_vector(const MatrixXd& bye_mat) {
    return bye_mat.rowwise().squaredNorm() / bye_mat.cols();
}

SparseMatrix<double> build_sparse_diagonal_matrix(const VectorXd& diagonal_values) {
    SparseMatrix<double> diagonal_matrix(diagonal_values.size(), diagonal_values.size());
    diagonal_matrix.reserve(Eigen::VectorXi::Constant(diagonal_values.size(), 1));
    for (std::int64_t i = 0; i < diagonal_values.size(); ++i) {
        diagonal_matrix.insert(i, i) = diagonal_values(i);
    }
    diagonal_matrix.makeCompressed();
    return diagonal_matrix;
}

ProjectionCache build_projection_cache(const SparseMatrix<double>& inverse_matrix,
        const MatrixXd& design_matrix, bool compute_logdet = false) {
    ProjectionCache projection_cache;
    projection_cache.inverse_times_design = inverse_matrix * design_matrix;
    const MatrixXd design_cross_product =
        design_matrix.transpose() * projection_cache.inverse_times_design;

    Eigen::LDLT<MatrixXd> ldlt(design_cross_product);
    if (ldlt.info() != Eigen::Success) {
        fatal_error("Failed to factorize X'V^{{-1}}X for the fixed-effect projection.");
    }

    if (compute_logdet) {
        projection_cache.design_logdet = ldlt.vectorD().array().log().sum();
    }
    projection_cache.design_cross_inverse =
        ldlt.solve(MatrixXd::Identity(design_cross_product.rows(), design_cross_product.cols()));
    return projection_cache;
}

ProjectionCache build_projection_cache(const VectorXd& inverse_diagonal,
        const MatrixXd& design_matrix, bool compute_logdet = false) {
    ProjectionCache projection_cache;
    projection_cache.inverse_times_design =
        (design_matrix.array().colwise() * inverse_diagonal.array()).matrix();
    const MatrixXd design_cross_product =
        design_matrix.transpose() * projection_cache.inverse_times_design;

    CustomLLT llt_solver;
    llt_solver.compute(design_cross_product);
    projection_cache.design_cross_inverse = llt_solver.inverse();
    if (compute_logdet) {
        projection_cache.design_logdet = llt_solver.logDeterminant();
    }
    return projection_cache;
}

ScalarProjectionCache build_scalar_projection_cache(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& design_vector) {
    ScalarProjectionCache projection_cache;
    projection_cache.inverse_times_design = inverse_matrix * design_vector;
    projection_cache.design_cross_inverse =
        1.0 / design_vector.dot(projection_cache.inverse_times_design);
    return projection_cache;
}

// Evaluates P·rhs = V⁻¹rhs - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹rhs without forming the n×n projection matrix P.
VectorXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const VectorXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

MatrixXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const MatrixXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

VectorXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& inverse_times_design, double design_cross_inverse, const VectorXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

MatrixXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& inverse_times_design, double design_cross_inverse, const MatrixXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

VectorXd apply_projection(const VectorXd& inverse_diagonal,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const VectorXd& rhs) {
    return inverse_diagonal.cwiseProduct(rhs) -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

MatrixXd apply_projection(const VectorXd& inverse_diagonal,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const MatrixXd& rhs) {
    MatrixXd projected_rhs = rhs.array().colwise() * inverse_diagonal.array();
    return projected_rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

double append_inverse_singleton_block(const MatrixXd& grm_block, const MatrixXd& gxe_block,
        const VectorXd& nxe_block, const VectorXd& varcom, bool no_noisebye,
        std::int64_t start_index, std::vector<GrmTriplet>& triplets) {
    if (grm_block.size() == 0) {
        return 0.0;
    }

    Eigen::ArrayXd variance_diag =
        (grm_block.array() * varcom(0) + gxe_block.array() * varcom(1) + varcom(varcom.size() - 1)).col(0);
    if (!no_noisebye) {
        variance_diag += nxe_block.array() * varcom(2);
    }

    const double logdet = variance_diag.log().sum();
    const MatrixXd inverse_block = (1.0 / variance_diag).matrix();
    append_singleton_block_triplets(inverse_block, start_index, triplets, false);
    return logdet;
}

double append_inverse_dense_block(const MatrixXd& grm_block, const MatrixXd& gxe_block,
        const VectorXd& nxe_block, const VectorXd& varcom, bool no_noisebye,
        std::int64_t start_index, std::vector<GrmTriplet>& triplets) {
    MatrixXd variance_block = MatrixXd::Identity(grm_block.rows(), grm_block.rows()) * varcom(varcom.size() - 1);
    variance_block += grm_block * varcom(0) + gxe_block * varcom(1);
    if (!no_noisebye) {
        variance_block += nxe_block.asDiagonal() * varcom(2);
    }

    Eigen::LDLT<MatrixXd> ldlt(variance_block);
    if (ldlt.info() != Eigen::Success) {
        fatal_error("Failed to factorize covariance block starting at index {}.", start_index);
    }

    const VectorXd diagonal = ldlt.vectorD();
    const double logdet = diagonal.array().log().sum();
    const MatrixXd inverse_block = ldlt.solve(MatrixXd::Identity(grm_block.rows(), grm_block.rows()));
    append_dense_block_triplets(inverse_block, start_index, triplets);
    return logdet;
}

InverseCovarianceResult build_inverse_covariance_matrix(const std::vector<MatrixXd>& grm_blocks,
        const std::vector<MatrixXd>& gxe_blocks, const std::vector<std::int64_t>& block_offsets,
        const VectorXd& nxe_vec, const VectorXd& varcom, bool no_noisebye, std::int64_t num_id) {
    InverseCovarianceResult result;
    std::vector<GrmTriplet> triplets;
    triplets.reserve(count_grm_triplets(grm_blocks));

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const MatrixXd& gxe_block = gxe_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        const VectorXd nxe_block = nxe_vec.segment(start_index, block_size);

        if (block_idx == 0) {
            result.logdet += append_inverse_singleton_block(
                grm_block, gxe_block, nxe_block, varcom, no_noisebye, start_index, triplets);
        } else {
            result.logdet += append_inverse_dense_block(
                grm_block, gxe_block, nxe_block, varcom, no_noisebye, start_index, triplets);
        }
    }

    result.inverse_matrix.resize(num_id, num_id);
    result.inverse_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

bool is_binary_trait_matrix(const MatrixXd& y, double tolerance = 1e-8) {
    if (y.cols() != 1 || y.rows() == 0) {
        return false;
    }

    for (std::int64_t i = 0; i < y.rows(); ++i) {
        const double value = y(i, 0);
        if (std::abs(value) <= tolerance || std::abs(value - 1.0) <= tolerance) {
            continue;
        }
        return false;
    }
    return true;
}

void require_binary_trait_matrix(const MatrixXd& y) {
    if (!is_binary_trait_matrix(y)) {
        fatal_error("Binary-trait logistic mixed model requires a single 0/1 phenotype column.");
    }
}

double logistic_scalar(double eta) {
    // Two-branch form avoids overflow: exp(-eta) stays finite when eta>=0; exp(eta) when eta<0.
    if (eta >= 0.0) {
        const double exp_neg_eta = std::exp(-eta);
        return 1.0 / (1.0 + exp_neg_eta);
    }
    const double exp_eta = std::exp(eta);
    return exp_eta / (1.0 + exp_eta);
}

VectorXd logistic_mean(const VectorXd& eta) {
    return eta.unaryExpr([](double value) { return logistic_scalar(value); });
}

VectorXd clamp_probability_vector(const VectorXd& mu, double epsilon = 1e-6) {
    return mu.array().min(1.0 - epsilon).max(epsilon).matrix();
}

double bernoulli_loglikelihood(const VectorXd& y, const VectorXd& mu) {
    return (y.array() * mu.array().log() +
            (1.0 - y.array()) * (1.0 - mu.array()).log()).sum();
}

VectorXd fit_logistic_fixed_effects(const MatrixXd& xmat, const VectorXd& y,
        int max_iter = 50, double tolerance = 1e-8) {
    VectorXd beta = VectorXd::Zero(xmat.cols());
    if (xmat.cols() > 0) {
        const double mean_y = y.mean();
        if (mean_y > 0.0 && mean_y < 1.0) {
            beta(0) = std::log(mean_y / (1.0 - mean_y));
        }
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        const VectorXd eta = xmat * beta;
        const VectorXd mu = clamp_probability_vector(logistic_mean(eta));
        const VectorXd weight = mu.array() * (1.0 - mu.array());
        const VectorXd working_response =
            eta.array() + (y.array() - mu.array()) / weight.array();

        MatrixXd weighted_x = xmat.array().colwise() * weight.array();
        MatrixXd xtwx = xmat.transpose() * weighted_x;
        VectorXd xtwz = xmat.transpose() * (weight.array() * working_response.array()).matrix();

        Eigen::LDLT<MatrixXd> ldlt(xtwx);
        if (ldlt.info() != Eigen::Success) {
            fatal_error("Failed to initialize the logistic fixed-effect model.");
        }

        const VectorXd beta_new = ldlt.solve(xtwz);
        if ((beta_new - beta).norm() < tolerance) {
            beta = beta_new;
            break;
        }
        beta = beta_new;
    }

    return beta;
}

double append_binary_inverse_singleton_block(const MatrixXd& grm_block, const VectorXd& residual_diag,
        double genetic_var, std::int64_t start_index, std::vector<GrmTriplet>& triplets) {
    if (grm_block.size() == 0) {
        return 0.0;
    }

    const Eigen::ArrayXd variance_diag =
        genetic_var * grm_block.col(0).array() + residual_diag.array();
    const double logdet = variance_diag.log().sum();
    const MatrixXd inverse_block = (1.0 / variance_diag).matrix();
    append_singleton_block_triplets(inverse_block, start_index, triplets, false);
    return logdet;
}

double append_binary_inverse_dense_block(const MatrixXd& grm_block, const VectorXd& residual_diag,
        double genetic_var, std::int64_t start_index, std::vector<GrmTriplet>& triplets) {
    MatrixXd variance_block = grm_block * genetic_var;
    variance_block.diagonal().array() += residual_diag.array();

    Eigen::LDLT<MatrixXd> ldlt(variance_block);
    if (ldlt.info() != Eigen::Success) {
        fatal_error("Failed to factorize binary-trait covariance block starting at index {}.", start_index);
    }

    const double logdet = ldlt.vectorD().array().log().sum();
    const MatrixXd inverse_block = ldlt.solve(MatrixXd::Identity(grm_block.rows(), grm_block.rows()));
    append_dense_block_triplets(inverse_block, start_index, triplets);
    return logdet;
}

InverseCovarianceResult build_main_binary_inverse_covariance_matrix(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& residual_diag,
        double genetic_var,
        std::int64_t num_id) {
    InverseCovarianceResult result;
    std::vector<GrmTriplet> triplets;
    triplets.reserve(count_grm_triplets(grm_blocks));

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        const VectorXd residual_block = residual_diag.segment(start_index, block_size);

        if (block_idx == 0) {
            result.logdet += append_binary_inverse_singleton_block(
                grm_block, residual_block, genetic_var, start_index, triplets);
        } else {
            result.logdet += append_binary_inverse_dense_block(
                grm_block, residual_block, genetic_var, start_index, triplets);
        }
    }

    result.inverse_matrix.resize(num_id, num_id);
    result.inverse_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

BinaryWorkingLmmState build_binary_working_lmm_state(const VectorXd& y, const VectorXd& eta) {
    BinaryWorkingLmmState state;
    state.mu = clamp_probability_vector(logistic_mean(eta));
    // mu*(1-mu) -> 0 as mu -> {0,1}; floor prevents division by zero in the IRLS weight matrix.
    state.weight = (state.mu.array() * (1.0 - state.mu.array())).max(kProbabilityWeightFloor).matrix();
    state.residual_diag = (1.0 / state.weight.array()).matrix();
    state.working_response =
        (eta.array() + (y.array() - state.mu.array()) / state.weight.array()).matrix();
    return state;
}

BinaryProjectionState evaluate_binary_projection_state(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const SparseMatrix<double>& grm_mat,
        const MatrixXd& xmat,
        const VectorXd& working_response,
        const VectorXd& residual_diag,
        double genetic_var,
        std::int64_t num_id) {
    BinaryProjectionState state;
    state.inverse_result = build_main_binary_inverse_covariance_matrix(
        grm_blocks, block_offsets, residual_diag, genetic_var, num_id);
    state.fixed_effect_projection = build_projection_cache(
        state.inverse_result.inverse_matrix, xmat, true);
    state.projected_response = apply_projection(
        state.inverse_result.inverse_matrix,
        state.fixed_effect_projection.inverse_times_design,
        state.fixed_effect_projection.design_cross_inverse,
        working_response);
    state.genetic_times_projected_response = grm_mat * state.projected_response;

    const double tr_ViG = (state.inverse_result.inverse_matrix.cwiseProduct(grm_mat)).sum();
    const MatrixXd design_cross_genetic_design =
        state.fixed_effect_projection.inverse_times_design.transpose() * grm_mat *
        state.fixed_effect_projection.inverse_times_design;
    state.trace_pg = tr_ViG - (state.fixed_effect_projection.design_cross_inverse
        .cwiseProduct(design_cross_genetic_design)).sum();
    state.projected_genetic_response = apply_projection(
        state.inverse_result.inverse_matrix,
        state.fixed_effect_projection.inverse_times_design,
        state.fixed_effect_projection.design_cross_inverse,
        state.genetic_times_projected_response);
    return state;
}

ScalarVarianceUpdateState find_positive_scalar_variance_update(
        double score,
        double ai,
        double em,
        double variance,
        double min_variance = 1e-6) {
    ScalarVarianceUpdateState result;
    result.score = score;
    result.ai = ai;
    result.em = em;

    const double safe_ai = (std::isfinite(ai) && ai > kAiEmFloor) ? ai : kAiEmFloor;
    const double safe_em = (std::isfinite(em) && em > kAiEmFloor) ? em : kAiEmFloor;

    // Blend AI (Newton) toward EM (Fisher) until the variance update stays positive; gamma=0 is pure AI.
    for (int step = 0; step <= 100; ++step) {
        const double gamma = step * 0.01;
        const double blended_information = (1.0 - gamma) * safe_ai + gamma * safe_em;
        if (!(std::isfinite(blended_information) && blended_information > 0.0)) {
            continue;
        }

        const double delta = score / blended_information;
        const double updated_variance = variance + delta;
        if (std::isfinite(updated_variance) && updated_variance > min_variance) {
            result.gamma = gamma;
            result.blended_information = blended_information;
            result.delta = delta;
            result.updated_variance = updated_variance;
            return result;
        }
    }

    result.gamma = 1.0;
    result.blended_information = safe_em;
    result.delta = score / safe_em;
    result.updated_variance = std::max(variance + result.delta, min_variance);
    return result;
}

double joint_genetic_covariance_bound(const VectorXd& varcom, double min_det = 1e-8) {
    const double margin = std::max(varcom(0) * varcom(2) - min_det, 1e-12);
    return std::sqrt(margin);
}

ScalarVarianceUpdateState find_bounded_scalar_update(
        double score,
        double ai,
        double em,
        double value,
        double lower_bound,
        double upper_bound) {
    ScalarVarianceUpdateState result;
    result.score = score;
    result.ai = ai;
    result.em = em;

    const double safe_ai = (std::isfinite(ai) && ai > kAiEmFloor) ? ai : kAiEmFloor;
    const double safe_em = (std::isfinite(em) && em > kAiEmFloor) ? em : kAiEmFloor;

    for (int step = 0; step <= 100; ++step) {
        const double gamma = step * 0.01;
        const double blended_information = (1.0 - gamma) * safe_ai + gamma * safe_em;
        if (!(std::isfinite(blended_information) && blended_information > 0.0)) {
            continue;
        }

        const double delta = score / blended_information;
        const double updated_value = value + delta;
        if (std::isfinite(updated_value) &&
                updated_value > lower_bound &&
                updated_value < upper_bound) {
            result.gamma = gamma;
            result.blended_information = blended_information;
            result.delta = delta;
            result.updated_variance = updated_value;
            return result;
        }
    }

    const double full_delta = score / safe_em;
    for (int shrink_step = 0; shrink_step <= 60; ++shrink_step) {
        const double shrink = std::pow(0.5, shrink_step);
        const double updated_value = value + shrink * full_delta;
        if (std::isfinite(updated_value) &&
                updated_value > lower_bound &&
                updated_value < upper_bound) {
            result.gamma = 1.0;
            result.blended_information = safe_em;
            result.delta = shrink * full_delta;
            result.updated_variance = updated_value;
            return result;
        }
    }

    result.gamma = 1.0;
    result.blended_information = safe_em;
    result.delta = 0.0;
    result.updated_variance = std::min(std::max(value, lower_bound + 1e-10), upper_bound - 1e-10);
    return result;
}

double compute_binary_working_neg2_reml(
        const VectorXd& working_response,
        const BinaryProjectionState& projection_state) {
    return projection_state.inverse_result.logdet +
           projection_state.fixed_effect_projection.design_logdet +
           working_response.dot(projection_state.projected_response);
}

BinaryVarianceFitState fit_binary_working_genetic_variance(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const SparseMatrix<double>& grm_mat,
        const MatrixXd& xmat,
        const VectorXd& working_response,
        const VectorXd& residual_diag,
        double initial_genetic_var,
        int maxiter,
        double cc_par,
        double cc_gra,
        double cc_logL,
        std::int64_t num_id,
        const std::string& log_prefix) {
    BinaryVarianceFitState result;
    result.genetic_var = std::max(initial_genetic_var, 1e-6);

    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < maxiter; ++iter) {
        spdlog::info("{} variance iteration {}", log_prefix, iter + 1);

        result.projection_state = evaluate_binary_projection_state(
            grm_blocks,
            block_offsets,
            grm_mat,
            xmat,
            working_response,
            residual_diag,
            result.genetic_var,
            num_id);
        result.neg2_reml = compute_binary_working_neg2_reml(working_response, result.projection_state);

        const double quad_pgp = result.projection_state.projected_response.dot(
            result.projection_state.genetic_times_projected_response);
        // The note writes the score for the negative objective as
        // 0.5 * [tr(PG) - y*'PGPy*]. Here we update by maximizing the working
        // quasi-likelihood, so the equivalent sign-flipped derivative is used.
        result.score = 0.5 * (quad_pgp - result.projection_state.trace_pg);
        result.ai = 0.5 * result.projection_state.genetic_times_projected_response.dot(
            result.projection_state.projected_genetic_response);
        result.em = num_id / (2.0 * result.genetic_var * result.genetic_var);

        const ScalarVarianceUpdateState update_state = find_positive_scalar_variance_update(
            result.score, result.ai, result.em, result.genetic_var);
        const double updated_variance = update_state.updated_variance;
        const double var_change =
            std::abs(updated_variance - result.genetic_var) /
            std::max(1.0, std::abs(result.genetic_var));
        const double logL_change = std::isfinite(previous_neg2_reml)
            ? std::abs(result.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();

        spdlog::info("{} -2 working REML: {}", log_prefix, result.neg2_reml);
        spdlog::info("{} Score: {}, AI: {}, EM: {}, gamma: {}",
            log_prefix, result.score, result.ai, result.em, update_state.gamma);
        spdlog::info("{} Updated genetic variance: {}", log_prefix, updated_variance);

        result.gamma = update_state.gamma;
        result.converged = (var_change < cc_par) ||
                           (std::abs(result.score) < cc_gra) ||
                           (logL_change < cc_logL);
        result.genetic_var = updated_variance;
        previous_neg2_reml = result.neg2_reml;

        if (result.converged) {
            break;
        }
    }

    result.projection_state = evaluate_binary_projection_state(
        grm_blocks,
        block_offsets,
        grm_mat,
        xmat,
        working_response,
        residual_diag,
        result.genetic_var,
        num_id);
    result.neg2_reml = compute_binary_working_neg2_reml(working_response, result.projection_state);
    const double quad_pgp = result.projection_state.projected_response.dot(
        result.projection_state.genetic_times_projected_response);
    result.score = 0.5 * (quad_pgp - result.projection_state.trace_pg);
    result.ai = 0.5 * result.projection_state.genetic_times_projected_response.dot(
        result.projection_state.projected_genetic_response);
    result.em = num_id / (2.0 * result.genetic_var * result.genetic_var);
    return result;
}

JointMainWorkingState build_joint_main_working_state(const MatrixXd& y, const VectorXd& eta_binary) {
    require_binary_continuous_trait_matrix(y);

    JointMainWorkingState state;
    state.binary_mu = clamp_probability_vector(logistic_mean(eta_binary));
    state.binary_weight = (state.binary_mu.array() * (1.0 - state.binary_mu.array())).max(1e-6).matrix();
    state.binary_residual_diag = (1.0 / state.binary_weight.array()).matrix();
    state.binary_working_response =
        (eta_binary.array() +
         (y.col(0).array() - state.binary_mu.array()) / state.binary_weight.array()).matrix();

    state.stacked_response.resize(y.rows() * 2);
    state.stacked_response.head(y.rows()) = state.binary_working_response;
    state.stacked_response.tail(y.rows()) = y.col(1);
    return state;
}

VectorXd apply_joint_variance_component(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& rhs,
        int component,
        std::int64_t num_id) {
    VectorXd result = VectorXd::Zero(rhs.size());

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const VectorXd rhs_binary = rhs.segment(start_index, block_size);
        const VectorXd rhs_continuous = rhs.segment(num_id + start_index, block_size);
        VectorXd gx_binary(block_size);
        VectorXd gx_continuous(block_size);

        if (block_idx == 0) {
            const VectorXd gdiag = grm_block.col(0);
            gx_binary = gdiag.cwiseProduct(rhs_binary);
            gx_continuous = gdiag.cwiseProduct(rhs_continuous);
        } else {
            gx_binary = grm_block * rhs_binary;
            gx_continuous = grm_block * rhs_continuous;
        }

        if (component == 0) {
            result.segment(start_index, block_size) = gx_binary;
        } else if (component == 1) {
            result.segment(start_index, block_size) = gx_continuous;
            result.segment(num_id + start_index, block_size) = gx_binary;
        } else if (component == 2) {
            result.segment(num_id + start_index, block_size) = gx_continuous;
        } else if (component == 3) {
            result.segment(num_id + start_index, block_size) = rhs_continuous;
        }
    }

    return result;
}

MatrixXd apply_joint_variance_component(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& rhs,
        int component,
        std::int64_t num_id) {
    MatrixXd result = MatrixXd::Zero(rhs.rows(), rhs.cols());

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const MatrixXd rhs_binary = rhs.middleRows(start_index, block_size);
        const MatrixXd rhs_continuous = rhs.middleRows(num_id + start_index, block_size);
        MatrixXd gx_binary(block_size, rhs.cols());
        MatrixXd gx_continuous(block_size, rhs.cols());

        if (block_idx == 0) {
            const VectorXd gdiag = grm_block.col(0);
            gx_binary = (rhs_binary.array().colwise() * gdiag.array()).matrix();
            gx_continuous = (rhs_continuous.array().colwise() * gdiag.array()).matrix();
        } else {
            gx_binary = grm_block * rhs_binary;
            gx_continuous = grm_block * rhs_continuous;
        }

        if (component == 0) {
            result.middleRows(start_index, block_size) = gx_binary;
        } else if (component == 1) {
            result.middleRows(start_index, block_size) = gx_continuous;
            result.middleRows(num_id + start_index, block_size) = gx_binary;
        } else if (component == 2) {
            result.middleRows(num_id + start_index, block_size) = gx_continuous;
        } else if (component == 3) {
            result.middleRows(num_id + start_index, block_size) = rhs_continuous;
        }
    }

    return result;
}

VectorXd apply_multitrait_continuous_variance_component(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& rhs,
        int component,
        std::int64_t num_id,
        std::int64_t num_traits) {
    VectorXd result = VectorXd::Zero(rhs.size());
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const bool is_residual_component = component >= num_covariance_terms;
    const auto [trait_row, trait_col] = lower_triangle_pair_from_index(
        component % num_covariance_terms, num_traits);

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const VectorXd rhs_row =
            rhs.segment(trait_row * num_id + start_index, block_size);
        const VectorXd rhs_col =
            rhs.segment(trait_col * num_id + start_index, block_size);

        VectorXd transformed_row(block_size);
        VectorXd transformed_col(block_size);
        if (is_residual_component) {
            transformed_row = rhs_row;
            transformed_col = rhs_col;
        } else if (block_idx == 0) {
            const VectorXd gdiag = grm_block.col(0);
            transformed_row = gdiag.cwiseProduct(rhs_row);
            transformed_col = gdiag.cwiseProduct(rhs_col);
        } else {
            transformed_row = grm_block * rhs_row;
            transformed_col = grm_block * rhs_col;
        }

        if (trait_row == trait_col) {
            result.segment(trait_row * num_id + start_index, block_size) = transformed_row;
        } else {
            result.segment(trait_row * num_id + start_index, block_size) = transformed_col;
            result.segment(trait_col * num_id + start_index, block_size) = transformed_row;
        }
    }

    return result;
}

MatrixXd apply_multitrait_continuous_variance_component(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& rhs,
        int component,
        std::int64_t num_id,
        std::int64_t num_traits) {
    MatrixXd result = MatrixXd::Zero(rhs.rows(), rhs.cols());
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const bool is_residual_component = component >= num_covariance_terms;
    const auto [trait_row, trait_col] = lower_triangle_pair_from_index(
        component % num_covariance_terms, num_traits);

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const MatrixXd rhs_row =
            rhs.middleRows(trait_row * num_id + start_index, block_size);
        const MatrixXd rhs_col =
            rhs.middleRows(trait_col * num_id + start_index, block_size);

        MatrixXd transformed_row(block_size, rhs.cols());
        MatrixXd transformed_col(block_size, rhs.cols());
        if (is_residual_component) {
            transformed_row = rhs_row;
            transformed_col = rhs_col;
        } else if (block_idx == 0) {
            const VectorXd gdiag = grm_block.col(0);
            transformed_row = (rhs_row.array().colwise() * gdiag.array()).matrix();
            transformed_col = (rhs_col.array().colwise() * gdiag.array()).matrix();
        } else {
            transformed_row = grm_block * rhs_row;
            transformed_col = grm_block * rhs_col;
        }

        if (trait_row == trait_col) {
            result.middleRows(trait_row * num_id + start_index, block_size) = transformed_row;
        } else {
            result.middleRows(trait_row * num_id + start_index, block_size) = transformed_col;
            result.middleRows(trait_col * num_id + start_index, block_size) = transformed_row;
        }
    }

    return result;
}

VectorXd apply_joint_genetic_random_effect(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& rhs,
        const VectorXd& varcom,
        std::int64_t num_id) {
    VectorXd result = VectorXd::Zero(rhs.size());

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const VectorXd rhs_binary = rhs.segment(start_index, block_size);
        const VectorXd rhs_continuous = rhs.segment(num_id + start_index, block_size);
        VectorXd gx_binary(block_size);
        VectorXd gx_continuous(block_size);

        if (block_idx == 0) {
            const VectorXd gdiag = grm_block.col(0);
            gx_binary = gdiag.cwiseProduct(rhs_binary);
            gx_continuous = gdiag.cwiseProduct(rhs_continuous);
        } else {
            gx_binary = grm_block * rhs_binary;
            gx_continuous = grm_block * rhs_continuous;
        }

        result.segment(start_index, block_size) =
            varcom(0) * gx_binary + varcom(1) * gx_continuous;
        result.segment(num_id + start_index, block_size) =
            varcom(1) * gx_binary + varcom(2) * gx_continuous;
    }

    return result;
}

double append_joint_inverse_singleton_block(
        const MatrixXd& grm_block,
        const VectorXd& residual_diag,
        const VectorXd& varcom,
        std::int64_t start_index,
        std::int64_t num_id,
        std::vector<GrmTriplet>& triplets,
        VectorXd& trace_vi_components) {
    if (grm_block.size() == 0) {
        return 0.0;
    }

    double logdet = 0.0;
    for (std::int64_t row_idx = 0; row_idx < grm_block.rows(); ++row_idx) {
        const double g_value = grm_block(row_idx, 0);
        const double variance_binary = varcom(0) * g_value + residual_diag(start_index + row_idx);
        const double covariance = varcom(1) * g_value;
        const double variance_continuous = varcom(2) * g_value + varcom(3);
        const double determinant = variance_binary * variance_continuous - covariance * covariance;
        if (!(determinant > 0.0)) {
            fatal_error("The joint covariance block for a singleton sample is not positive definite.");
        }

        const double inverse_binary = variance_continuous / determinant;
        const double inverse_cross = -covariance / determinant;
        const double inverse_continuous = variance_binary / determinant;
        const std::int64_t binary_index = start_index + row_idx;
        const std::int64_t continuous_index = num_id + binary_index;

        triplets.emplace_back(binary_index, binary_index, inverse_binary);
        triplets.emplace_back(binary_index, continuous_index, inverse_cross);
        triplets.emplace_back(continuous_index, binary_index, inverse_cross);
        triplets.emplace_back(continuous_index, continuous_index, inverse_continuous);

        trace_vi_components(0) += inverse_binary * g_value;
        trace_vi_components(1) += 2.0 * inverse_cross * g_value;
        trace_vi_components(2) += inverse_continuous * g_value;
        trace_vi_components(3) += inverse_continuous;
        logdet += std::log(determinant);
    }

    return logdet;
}

double append_joint_inverse_dense_block(
        const MatrixXd& grm_block,
        const VectorXd& residual_diag,
        const VectorXd& varcom,
        std::int64_t start_index,
        std::int64_t num_id,
        std::vector<GrmTriplet>& triplets,
        VectorXd& trace_vi_components) {
    const std::int64_t block_size = grm_block.rows();
    MatrixXd covariance_block = MatrixXd::Zero(block_size * 2, block_size * 2);
    covariance_block.topLeftCorner(block_size, block_size) = grm_block * varcom(0);
    covariance_block.topLeftCorner(block_size, block_size).diagonal().array() += residual_diag.array();
    covariance_block.topRightCorner(block_size, block_size) = grm_block * varcom(1);
    covariance_block.bottomLeftCorner(block_size, block_size) = grm_block * varcom(1);
    covariance_block.bottomRightCorner(block_size, block_size) = grm_block * varcom(2);
    covariance_block.bottomRightCorner(block_size, block_size).diagonal().array() += varcom(3);

    Eigen::LDLT<MatrixXd> ldlt(covariance_block);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0) {
        fatal_error("Failed to factorize binary-continuous covariance block starting at index {}.", start_index);
    }

    const MatrixXd inverse_block =
        ldlt.solve(MatrixXd::Identity(covariance_block.rows(), covariance_block.cols()));
    const MatrixXd inverse_binary =
        inverse_block.topLeftCorner(block_size, block_size);
    const MatrixXd inverse_cross =
        inverse_block.topRightCorner(block_size, block_size);
    const MatrixXd inverse_continuous =
        inverse_block.bottomRightCorner(block_size, block_size);

    trace_vi_components(0) += (inverse_binary.cwiseProduct(grm_block)).sum();
    trace_vi_components(1) += 2.0 * (inverse_cross.cwiseProduct(grm_block)).sum();
    trace_vi_components(2) += (inverse_continuous.cwiseProduct(grm_block)).sum();
    trace_vi_components(3) += inverse_continuous.diagonal().sum();

    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        for (std::int64_t col_idx = 0; col_idx < block_size; ++col_idx) {
            const std::int64_t binary_row = start_index + row_idx;
            const std::int64_t binary_col = start_index + col_idx;
            const std::int64_t continuous_row = num_id + binary_row;
            const std::int64_t continuous_col = num_id + binary_col;

            triplets.emplace_back(binary_row, binary_col, inverse_binary(row_idx, col_idx));
            triplets.emplace_back(binary_row, continuous_col, inverse_cross(row_idx, col_idx));
            triplets.emplace_back(continuous_row, binary_col,
                inverse_block(block_size + row_idx, col_idx));
            triplets.emplace_back(continuous_row, continuous_col, inverse_continuous(row_idx, col_idx));
        }
    }

    return ldlt.vectorD().array().log().sum();
}

JointInverseCovarianceResult build_joint_main_inverse_covariance_matrix(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& residual_diag,
        const VectorXd& varcom,
        std::int64_t num_id) {
    JointInverseCovarianceResult result;
    result.trace_vi_components = VectorXd::Zero(4);

    std::vector<GrmTriplet> triplets;
    triplets.reserve(4 * count_grm_triplets(grm_blocks));

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        const VectorXd residual_block = residual_diag.segment(start_index, block_size);
        if (block_idx == 0) {
            result.logdet += append_joint_inverse_singleton_block(
                grm_block, residual_block, varcom, start_index, num_id, triplets, result.trace_vi_components);
        } else {
            result.logdet += append_joint_inverse_dense_block(
                grm_block, residual_block, varcom, start_index, num_id, triplets, result.trace_vi_components);
        }
    }

    result.inverse_matrix.resize(num_id * 2, num_id * 2);
    result.inverse_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

MatrixXd build_bivariate_em_information_block(
        double var_trait1,
        double covariance,
        double var_trait2,
        std::int64_t num_id) {
    MatrixXd covariance_block = MatrixXd::Zero(3, 3);
    covariance_block <<
        2.0 * std::pow(var_trait1, 2), 2.0 * var_trait1 * covariance, 2.0 * std::pow(covariance, 2),
        2.0 * var_trait1 * covariance, var_trait1 * var_trait2 + std::pow(covariance, 2), 2.0 * covariance * var_trait2,
        2.0 * std::pow(covariance, 2), 2.0 * covariance * var_trait2, 2.0 * std::pow(var_trait2, 2);
    covariance_block /= (2.0 * static_cast<double>(num_id));
    covariance_block.diagonal().array() += kEmStabilityNudge;

    Eigen::LDLT<MatrixXd> ldlt(covariance_block);
    if (ldlt.info() != Eigen::Success) {
        fatal_error("Failed to build the bivariate EM information block.");
    }

    return 0.5 * ldlt.solve(MatrixXd::Identity(3, 3));
}

MatrixXd build_joint_em_information(const VectorXd& varcom, std::int64_t num_id) {
    MatrixXd em = MatrixXd::Zero(4, 4);
    em.topLeftCorner(3, 3) = build_bivariate_em_information_block(
        varcom(0), varcom(1), varcom(2), num_id);
    em(3, 3) = num_id / (2.0 * std::pow(varcom(3), 2));
    return em;
}

MatrixXd build_lower_triangle_basis_matrix(std::int64_t matrix_size, std::int64_t packed_index) {
    MatrixXd basis_matrix = MatrixXd::Zero(matrix_size, matrix_size);
    const auto [row, col] = lower_triangle_pair_from_index(packed_index, matrix_size);
    basis_matrix(row, col) = 1.0;
    basis_matrix(col, row) = 1.0;
    return basis_matrix;
}

MatrixXd build_covariance_matrix_em_information(
        const MatrixXd& covariance_matrix,
        std::int64_t num_id) {
    const std::int64_t matrix_size = covariance_matrix.rows();
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(matrix_size);
    MatrixXd em_information = MatrixXd::Zero(num_covariance_terms, num_covariance_terms);

    Eigen::LDLT<MatrixXd> ldlt(covariance_matrix);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0) {
        fatal_error("Failed to build the multitrait continuous EM information block.");
    }

    const MatrixXd inverse_covariance =
        ldlt.solve(MatrixXd::Identity(matrix_size, matrix_size));
    std::vector<MatrixXd> inverse_basis_products(num_covariance_terms);
    for (std::int64_t component = 0; component < num_covariance_terms; ++component) {
        inverse_basis_products[component] =
            inverse_covariance * build_lower_triangle_basis_matrix(matrix_size, component);
    }

    for (std::int64_t row = 0; row < num_covariance_terms; ++row) {
        for (std::int64_t col = 0; col <= row; ++col) {
            const double information =
                0.5 * static_cast<double>(num_id) *
                (inverse_basis_products[row] * inverse_basis_products[col]).trace();
            em_information(row, col) = em_information(col, row) = information;
        }
    }

    em_information.diagonal().array() += 1e-10;
    return em_information;
}

MatrixXd build_multitrait_continuous_em_information(
        const VectorXd& varcom,
        std::int64_t num_id,
        std::int64_t num_traits) {
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    MatrixXd em_information = MatrixXd::Zero(2 * num_covariance_terms, 2 * num_covariance_terms);
    em_information.topLeftCorner(num_covariance_terms, num_covariance_terms) =
        build_covariance_matrix_em_information(
            unpack_lower_triangle_to_symmetric_matrix(varcom, num_traits, 0),
            num_id);
    em_information.bottomRightCorner(num_covariance_terms, num_covariance_terms) =
        build_covariance_matrix_em_information(
            unpack_lower_triangle_to_symmetric_matrix(varcom, num_traits, num_covariance_terms),
            num_id);
    return em_information;
}

VarianceUpdateResult find_valid_joint_variance_update(
        const MatrixXd& ai_mat,
        const MatrixXd& em_mat,
        const VectorXd& score,
        const VectorXd& varcom) {
    VarianceUpdateResult update_result;

    for (int step = 0; step <= 100; ++step) {
        const double gamma = step * 0.01;
        const MatrixXd blended_information = (1.0 - gamma) * ai_mat + gamma * em_mat;
        Eigen::LDLT<MatrixXd> ldlt(blended_information);
        if (ldlt.info() != Eigen::Success) {
            continue;
        }

        update_result.delta = ldlt.solve(score);
        update_result.updated_varcom = varcom + update_result.delta;
        if (is_valid_joint_main_variance_components(update_result.updated_varcom)) {
            spdlog::info("Binary-continuous EM weight value: {}", gamma);
            return update_result;
        }
    }

    Eigen::LDLT<MatrixXd> em_ldlt(em_mat);
    if (em_ldlt.info() != Eigen::Success) {
        fatal_error("Failed to find a valid AI/EM update for the binary-continuous null model.");
    }

    const VectorXd full_step = em_ldlt.solve(score);
    for (int shrink_step = 0; shrink_step <= 60; ++shrink_step) {
        const double shrink = std::pow(0.5, shrink_step);
        update_result.delta = shrink * full_step;
        update_result.updated_varcom = varcom + update_result.delta;
        if (is_valid_joint_main_variance_components(update_result.updated_varcom)) {
            spdlog::info("Binary-continuous EM fallback shrink factor: {}", shrink);
            return update_result;
        }
    }

    fatal_error("Unable to find a valid variance-component update for the binary-continuous null model.");
}

JointProjectionState evaluate_joint_main_projection_state(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& joint_design,
        const VectorXd& working_response,
        const VectorXd& residual_diag,
        const VectorXd& varcom,
        std::int64_t num_id) {
    JointProjectionState state;
    state.inverse_result = build_joint_main_inverse_covariance_matrix(
        grm_blocks, block_offsets, residual_diag, varcom, num_id);
    state.fixed_effect_projection = build_projection_cache(
        state.inverse_result.inverse_matrix, joint_design, true);
    state.projected_response = apply_projection(
        state.inverse_result.inverse_matrix,
        state.fixed_effect_projection.inverse_times_design,
        state.fixed_effect_projection.design_cross_inverse,
        working_response);
    state.projected_random_response = apply_joint_genetic_random_effect(
        grm_blocks, block_offsets, state.projected_response, varcom, num_id);
    state.neg2_reml = state.inverse_result.logdet +
                      state.fixed_effect_projection.design_logdet +
                      working_response.dot(state.projected_response);
    return state;
}

JointVarianceFitState fit_joint_main_working_variance(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& joint_design,
        const VectorXd& working_response,
        const VectorXd& residual_diag,
        const VectorXd& initial_varcom,
        int maxiter,
        double cc_par,
        double cc_gra,
        double cc_logL,
        std::int64_t num_id,
        const std::string& log_prefix) {
    JointVarianceFitState result;
    result.varcom = initial_varcom;

    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < maxiter; ++iter) {
        spdlog::info("{} variance iteration {}", log_prefix, iter + 1);

        result.projection_state = evaluate_joint_main_projection_state(
            grm_blocks, block_offsets, joint_design, working_response, residual_diag, result.varcom, num_id);
        result.neg2_reml = result.projection_state.neg2_reml;

        std::array<VectorXd, 4> kpy_vec;
        std::array<VectorXd, 4> pkpy_vec;
        for (int component = 0; component < 4; ++component) {
            kpy_vec[component] = apply_joint_variance_component(
                grm_blocks, block_offsets, result.projection_state.projected_response, component, num_id);
            pkpy_vec[component] = apply_projection(
                result.projection_state.inverse_result.inverse_matrix,
                result.projection_state.fixed_effect_projection.inverse_times_design,
                result.projection_state.fixed_effect_projection.design_cross_inverse,
                kpy_vec[component]);

            const MatrixXd design_cross_component =
                result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
                apply_joint_variance_component(
                    grm_blocks,
                    block_offsets,
                    result.projection_state.fixed_effect_projection.inverse_times_design,
                    component,
                    num_id);
            const double trace_projection_component =
                result.projection_state.inverse_result.trace_vi_components(component) -
                (result.projection_state.fixed_effect_projection.design_cross_inverse
                    .cwiseProduct(design_cross_component)).sum();
            result.score(component) = 0.5 * (
                result.projection_state.projected_response.dot(kpy_vec[component]) -
                trace_projection_component);
        }

        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col <= row; ++col) {
                result.ai(row, col) = result.ai(col, row) =
                    0.5 * kpy_vec[row].dot(pkpy_vec[col]);
            }
        }
        result.em = build_joint_em_information(result.varcom, num_id);

        const VarianceUpdateResult update_result = find_valid_joint_variance_update(
            result.ai, result.em, result.score, result.varcom);
        const double cc_par_val =
            std::sqrt(update_result.delta.squaredNorm() /
                      std::max(1.0, result.varcom.squaredNorm()));
        const double cc_gra_val = result.score.norm();
        const double cc_logL_val = std::isfinite(previous_neg2_reml)
            ? std::abs(result.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();
        const bool parameter_stable = cc_par_val < cc_par;
        const bool gradient_small = cc_gra_val < cc_gra;
        const bool objective_stable = cc_logL_val < cc_logL;

        spdlog::info("{} -2 working REML: {}", log_prefix, result.neg2_reml);
        spdlog::info("{} Score: {}", log_prefix, format_vector_for_log(result.score));
        spdlog::info("{} Updated variances: {}", log_prefix,
            format_vector_for_log(update_result.updated_varcom));
        spdlog::info(
            "{} Convergence metrics: cc_par={}, cc_gra={}, cc_logL={}",
            log_prefix, cc_par_val, cc_gra_val, cc_logL_val);

        result.varcom = update_result.updated_varcom;
        // For the joint model, a tiny step or a flat objective alone is not
        // enough: the score must also be small before we accept convergence.
        result.converged = gradient_small && (parameter_stable || objective_stable);
        previous_neg2_reml = result.neg2_reml;
        if (result.converged) {
            break;
        }
    }

    result.projection_state = evaluate_joint_main_projection_state(
        grm_blocks, block_offsets, joint_design, working_response, residual_diag, result.varcom, num_id);
    result.neg2_reml = result.projection_state.neg2_reml;
    result.score.setZero();
    result.ai.setZero();
    result.em = build_joint_em_information(result.varcom, num_id);

    std::array<VectorXd, 4> final_kpy_vec;
    std::array<VectorXd, 4> final_pkpy_vec;
    for (int component = 0; component < 4; ++component) {
        final_kpy_vec[component] = apply_joint_variance_component(
            grm_blocks, block_offsets, result.projection_state.projected_response, component, num_id);
        final_pkpy_vec[component] = apply_projection(
            result.projection_state.inverse_result.inverse_matrix,
            result.projection_state.fixed_effect_projection.inverse_times_design,
            result.projection_state.fixed_effect_projection.design_cross_inverse,
            final_kpy_vec[component]);

        const MatrixXd design_cross_component =
            result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
            apply_joint_variance_component(
                grm_blocks,
                block_offsets,
                result.projection_state.fixed_effect_projection.inverse_times_design,
                component,
                num_id);
        const double trace_projection_component =
            result.projection_state.inverse_result.trace_vi_components(component) -
            (result.projection_state.fixed_effect_projection.design_cross_inverse
                .cwiseProduct(design_cross_component)).sum();
        result.score(component) = 0.5 * (
            result.projection_state.projected_response.dot(final_kpy_vec[component]) -
            trace_projection_component);
    }

    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col <= row; ++col) {
            result.ai(row, col) = result.ai(col, row) =
                0.5 * final_kpy_vec[row].dot(final_pkpy_vec[col]);
        }
    }

    return result;
}

JointCovarianceFitState fit_joint_main_working_genetic_covariance(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& joint_design,
        const VectorXd& working_response,
        const VectorXd& residual_diag,
        const VectorXd& fixed_varcom,
        int maxiter,
        double cc_par,
        double cc_gra,
        double cc_logL,
        std::int64_t num_id,
        const std::string& log_prefix) {
    JointCovarianceFitState result;
    result.varcom = fixed_varcom;

    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < maxiter; ++iter) {
        spdlog::info("{} covariance iteration {}", log_prefix, iter + 1);

        result.projection_state = evaluate_joint_main_projection_state(
            grm_blocks, block_offsets, joint_design, working_response, residual_diag, result.varcom, num_id);
        result.neg2_reml = result.projection_state.neg2_reml;

        const VectorXd kpy = apply_joint_variance_component(
            grm_blocks, block_offsets, result.projection_state.projected_response, 1, num_id);
        const VectorXd pkpy = apply_projection(
            result.projection_state.inverse_result.inverse_matrix,
            result.projection_state.fixed_effect_projection.inverse_times_design,
            result.projection_state.fixed_effect_projection.design_cross_inverse,
            kpy);

        const MatrixXd design_cross_component =
            result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
            apply_joint_variance_component(
                grm_blocks,
                block_offsets,
                result.projection_state.fixed_effect_projection.inverse_times_design,
                1,
                num_id);
        const double trace_projection_component =
            result.projection_state.inverse_result.trace_vi_components(1) -
            (result.projection_state.fixed_effect_projection.design_cross_inverse
                .cwiseProduct(design_cross_component)).sum();
        result.score = 0.5 * (
            result.projection_state.projected_response.dot(kpy) -
            trace_projection_component);
        result.ai = 0.5 * kpy.dot(pkpy);
        result.em = build_joint_em_information(result.varcom, num_id)(1, 1);

        const double covariance_bound = joint_genetic_covariance_bound(result.varcom);
        const ScalarVarianceUpdateState update_state = find_bounded_scalar_update(
            result.score,
            result.ai,
            result.em,
            result.varcom(1),
            -covariance_bound,
            covariance_bound);
        const double updated_covariance = update_state.updated_variance;
        const double cc_par_val =
            std::abs(updated_covariance - result.varcom(1)) /
            std::max(1.0, std::abs(result.varcom(1)));
        const double cc_gra_val = std::abs(result.score);
        const double cc_logL_val = std::isfinite(previous_neg2_reml)
            ? std::abs(result.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();
        const bool parameter_stable = cc_par_val < cc_par;
        const bool gradient_small = cc_gra_val < cc_gra;
        const bool objective_stable = cc_logL_val < cc_logL;

        spdlog::info("{} -2 working REML: {}", log_prefix, result.neg2_reml);
        spdlog::info("{} Covariance score: {}, AI: {}, EM: {}, gamma: {}",
            log_prefix, result.score, result.ai, result.em, update_state.gamma);
        spdlog::info("{} Updated genetic covariance: {}", log_prefix, updated_covariance);
        spdlog::info(
            "{} Convergence metrics: cc_par={}, cc_gra={}, cc_logL={}",
            log_prefix, cc_par_val, cc_gra_val, cc_logL_val);

        result.gamma = update_state.gamma;
        result.varcom(1) = updated_covariance;
        result.converged = gradient_small && (parameter_stable || objective_stable);
        previous_neg2_reml = result.neg2_reml;
        if (result.converged) {
            break;
        }
    }

    result.projection_state = evaluate_joint_main_projection_state(
        grm_blocks, block_offsets, joint_design, working_response, residual_diag, result.varcom, num_id);
    result.neg2_reml = result.projection_state.neg2_reml;

    const VectorXd final_kpy = apply_joint_variance_component(
        grm_blocks, block_offsets, result.projection_state.projected_response, 1, num_id);
    const VectorXd final_pkpy = apply_projection(
        result.projection_state.inverse_result.inverse_matrix,
        result.projection_state.fixed_effect_projection.inverse_times_design,
        result.projection_state.fixed_effect_projection.design_cross_inverse,
        final_kpy);
    const MatrixXd final_design_cross_component =
        result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
        apply_joint_variance_component(
            grm_blocks,
            block_offsets,
            result.projection_state.fixed_effect_projection.inverse_times_design,
            1,
            num_id);
    const double final_trace_projection_component =
        result.projection_state.inverse_result.trace_vi_components(1) -
        (result.projection_state.fixed_effect_projection.design_cross_inverse
            .cwiseProduct(final_design_cross_component)).sum();
    result.score = 0.5 * (
        result.projection_state.projected_response.dot(final_kpy) -
        final_trace_projection_component);
    result.ai = 0.5 * final_kpy.dot(final_pkpy);
    result.em = build_joint_em_information(result.varcom, num_id)(1, 1);
    return result;
}

double append_multitrait_continuous_inverse_singleton_block(
        const MatrixXd& grm_block,
        const MatrixXd& genetic_covariance,
        const MatrixXd& residual_covariance,
        std::int64_t start_index,
        std::int64_t num_id,
        std::int64_t num_traits,
        std::vector<GrmTriplet>& triplets,
        VectorXd& trace_vi_components) {
    if (grm_block.size() == 0) {
        return 0.0;
    }

    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    double logdet = 0.0;
    for (std::int64_t row_idx = 0; row_idx < grm_block.rows(); ++row_idx) {
        const double g_value = grm_block(row_idx, 0);
        const MatrixXd covariance_block =
            g_value * genetic_covariance + residual_covariance;
        Eigen::LDLT<MatrixXd> ldlt(covariance_block);
        if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0) {
            fatal_error("The multitrait continuous covariance block for a singleton sample is not positive definite.");
        }

        const MatrixXd inverse_block =
            ldlt.solve(MatrixXd::Identity(num_traits, num_traits));
        for (std::int64_t trait_row = 0; trait_row < num_traits; ++trait_row) {
            for (std::int64_t trait_col = 0; trait_col < num_traits; ++trait_col) {
                triplets.emplace_back(
                    trait_row * num_id + start_index + row_idx,
                    trait_col * num_id + start_index + row_idx,
                    inverse_block(trait_row, trait_col));
            }
        }

        for (std::int64_t component = 0; component < num_covariance_terms; ++component) {
            const auto [trait_row, trait_col] =
                lower_triangle_pair_from_index(component, num_traits);
            const double trace_contribution = (trait_row == trait_col)
                ? inverse_block(trait_row, trait_col)
                : (inverse_block(trait_row, trait_col) + inverse_block(trait_col, trait_row));
            trace_vi_components(component) += trace_contribution * g_value;
            trace_vi_components(num_covariance_terms + component) += trace_contribution;
        }
        logdet += ldlt.vectorD().array().log().sum();
    }

    return logdet;
}

double append_multitrait_continuous_inverse_dense_block(
        const MatrixXd& grm_block,
        const MatrixXd& genetic_covariance,
        const MatrixXd& residual_covariance,
        std::int64_t start_index,
        std::int64_t num_id,
        std::int64_t num_traits,
        std::vector<GrmTriplet>& triplets,
        VectorXd& trace_vi_components) {
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const std::int64_t block_size = grm_block.rows();
    MatrixXd covariance_block = MatrixXd::Zero(block_size * num_traits, block_size * num_traits);

    for (std::int64_t trait_row = 0; trait_row < num_traits; ++trait_row) {
        for (std::int64_t trait_col = 0; trait_col < num_traits; ++trait_col) {
            MatrixXd trait_covariance_block = grm_block * genetic_covariance(trait_row, trait_col);
            trait_covariance_block.diagonal().array() += residual_covariance(trait_row, trait_col);
            covariance_block.block(
                trait_row * block_size,
                trait_col * block_size,
                block_size,
                block_size) = trait_covariance_block;
        }
    }

    Eigen::LDLT<MatrixXd> ldlt(covariance_block);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0.0) {
        fatal_error("Failed to factorize a multitrait continuous covariance block starting at index {}.",
            start_index);
    }

    const MatrixXd inverse_block =
        ldlt.solve(MatrixXd::Identity(covariance_block.rows(), covariance_block.cols()));

    for (std::int64_t component = 0; component < num_covariance_terms; ++component) {
        const auto [trait_row, trait_col] =
            lower_triangle_pair_from_index(component, num_traits);
        const MatrixXd inverse_row_col = inverse_block.block(
            trait_row * block_size,
            trait_col * block_size,
            block_size,
            block_size);
        if (trait_row == trait_col) {
            trace_vi_components(component) +=
                (inverse_row_col.cwiseProduct(grm_block)).sum();
            trace_vi_components(num_covariance_terms + component) +=
                inverse_row_col.diagonal().sum();
        } else {
            const MatrixXd inverse_col_row = inverse_block.block(
                trait_col * block_size,
                trait_row * block_size,
                block_size,
                block_size);
            trace_vi_components(component) +=
                (inverse_row_col.cwiseProduct(grm_block)).sum() +
                (inverse_col_row.cwiseProduct(grm_block)).sum();
            trace_vi_components(num_covariance_terms + component) +=
                inverse_row_col.diagonal().sum() +
                inverse_col_row.diagonal().sum();
        }
    }

    for (std::int64_t trait_row = 0; trait_row < num_traits; ++trait_row) {
        for (std::int64_t trait_col = 0; trait_col < num_traits; ++trait_col) {
            const MatrixXd inverse_trait_block = inverse_block.block(
                trait_row * block_size,
                trait_col * block_size,
                block_size,
                block_size);
            for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
                for (std::int64_t col_idx = 0; col_idx < block_size; ++col_idx) {
                    triplets.emplace_back(
                        trait_row * num_id + start_index + row_idx,
                        trait_col * num_id + start_index + col_idx,
                        inverse_trait_block(row_idx, col_idx));
                }
            }
        }
    }

    return ldlt.vectorD().array().log().sum();
}

MultitraitContinuousInverseCovarianceResult build_multitrait_continuous_inverse_covariance_matrix(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const VectorXd& varcom,
        std::int64_t num_id,
        std::int64_t num_traits) {
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const MatrixXd genetic_covariance =
        unpack_lower_triangle_to_symmetric_matrix(varcom, num_traits, 0);
    const MatrixXd residual_covariance =
        unpack_lower_triangle_to_symmetric_matrix(varcom, num_traits, num_covariance_terms);

    MultitraitContinuousInverseCovarianceResult result;
    result.trace_vi_components = VectorXd::Zero(2 * num_covariance_terms);

    std::vector<GrmTriplet> triplets;
    triplets.reserve(
        static_cast<std::size_t>(num_traits * num_traits) *
        count_grm_triplets(grm_blocks));

    for (std::size_t block_idx = 0; block_idx < grm_blocks.size(); ++block_idx) {
        const MatrixXd& grm_block = grm_blocks[block_idx];
        const std::int64_t start_index = block_offsets[block_idx];
        const std::int64_t block_size = grm_block.rows();
        if (block_size == 0) {
            continue;
        }

        if (block_idx == 0) {
            result.logdet += append_multitrait_continuous_inverse_singleton_block(
                grm_block,
                genetic_covariance,
                residual_covariance,
                start_index,
                num_id,
                num_traits,
                triplets,
                result.trace_vi_components);
        } else {
            result.logdet += append_multitrait_continuous_inverse_dense_block(
                grm_block,
                genetic_covariance,
                residual_covariance,
                start_index,
                num_id,
                num_traits,
                triplets,
                result.trace_vi_components);
        }
    }

    result.inverse_matrix.resize(num_id * num_traits, num_id * num_traits);
    result.inverse_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

VarianceUpdateResult find_valid_multitrait_continuous_variance_update(
        const MatrixXd& ai_mat,
        const MatrixXd& em_mat,
        const VectorXd& score,
        const VectorXd& varcom,
        std::int64_t num_traits) {
    VarianceUpdateResult update_result;

    for (int step = 0; step <= 100; ++step) {
        const double gamma = step * 0.01;
        const MatrixXd blended_information = (1.0 - gamma) * ai_mat + gamma * em_mat;
        Eigen::LDLT<MatrixXd> ldlt(blended_information);
        if (ldlt.info() != Eigen::Success) {
            continue;
        }

        update_result.delta = ldlt.solve(score);
        update_result.updated_varcom = varcom + update_result.delta;
        if (is_valid_multitrait_continuous_variance_components(
                update_result.updated_varcom, num_traits)) {
            spdlog::info("Multitrait continuous EM weight value: {}", gamma);
            return update_result;
        }
    }

    Eigen::LDLT<MatrixXd> em_ldlt(em_mat);
    if (em_ldlt.info() != Eigen::Success) {
        fatal_error("Failed to find a valid AI/EM update for the multitrait continuous null model.");
    }

    const VectorXd full_step = em_ldlt.solve(score);
    for (int shrink_step = 0; shrink_step <= 60; ++shrink_step) {
        const double shrink = std::pow(0.5, shrink_step);
        update_result.delta = shrink * full_step;
        update_result.updated_varcom = varcom + update_result.delta;
        if (is_valid_multitrait_continuous_variance_components(
                update_result.updated_varcom, num_traits)) {
            spdlog::info("Multitrait continuous EM fallback shrink factor: {}", shrink);
            return update_result;
        }
    }

    fatal_error("Unable to find a valid variance-component update for the multitrait continuous null model.");
}

MultitraitContinuousProjectionState evaluate_multitrait_continuous_projection_state(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& joint_design,
        const VectorXd& stacked_response,
        const VectorXd& varcom,
        std::int64_t num_id,
        std::int64_t num_traits) {
    MultitraitContinuousProjectionState state;
    state.inverse_result = build_multitrait_continuous_inverse_covariance_matrix(
        grm_blocks, block_offsets, varcom, num_id, num_traits);
    state.fixed_effect_projection = build_projection_cache(
        state.inverse_result.inverse_matrix, joint_design, true);
    state.projected_response = apply_projection(
        state.inverse_result.inverse_matrix,
        state.fixed_effect_projection.inverse_times_design,
        state.fixed_effect_projection.design_cross_inverse,
        stacked_response);
    state.neg2_reml = state.inverse_result.logdet +
                      state.fixed_effect_projection.design_logdet +
                      stacked_response.dot(state.projected_response);
    return state;
}

MultitraitContinuousVarianceFitState fit_multitrait_continuous_variance(
        const std::vector<MatrixXd>& grm_blocks,
        const std::vector<std::int64_t>& block_offsets,
        const MatrixXd& joint_design,
        const VectorXd& stacked_response,
        const VectorXd& initial_varcom,
        int maxiter,
        double cc_par,
        double cc_gra,
        double cc_logL,
        std::int64_t num_id,
        std::int64_t num_traits,
        const std::string& log_prefix) {
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const std::int64_t num_variance_terms = 2 * num_covariance_terms;

    MultitraitContinuousVarianceFitState result;
    result.varcom = initial_varcom;
    result.score = VectorXd::Zero(num_variance_terms);
    result.ai = MatrixXd::Zero(num_variance_terms, num_variance_terms);
    result.em = MatrixXd::Zero(num_variance_terms, num_variance_terms);

    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < maxiter; ++iter) {
        spdlog::info("{} variance iteration {}", log_prefix, iter + 1);

        result.projection_state = evaluate_multitrait_continuous_projection_state(
            grm_blocks,
            block_offsets,
            joint_design,
            stacked_response,
            result.varcom,
            num_id,
            num_traits);
        result.neg2_reml = result.projection_state.neg2_reml;

        std::vector<VectorXd> kpy_vec(num_variance_terms);
        std::vector<VectorXd> pkpy_vec(num_variance_terms);
        for (int component = 0; component < num_variance_terms; ++component) {
            kpy_vec[component] = apply_multitrait_continuous_variance_component(
                grm_blocks,
                block_offsets,
                result.projection_state.projected_response,
                component,
                num_id,
                num_traits);
            pkpy_vec[component] = apply_projection(
                result.projection_state.inverse_result.inverse_matrix,
                result.projection_state.fixed_effect_projection.inverse_times_design,
                result.projection_state.fixed_effect_projection.design_cross_inverse,
                kpy_vec[component]);

            const MatrixXd design_cross_component =
                result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
                apply_multitrait_continuous_variance_component(
                    grm_blocks,
                    block_offsets,
                    result.projection_state.fixed_effect_projection.inverse_times_design,
                    component,
                    num_id,
                    num_traits);
            const double trace_projection_component =
                result.projection_state.inverse_result.trace_vi_components(component) -
                (result.projection_state.fixed_effect_projection.design_cross_inverse
                    .cwiseProduct(design_cross_component)).sum();
            result.score(component) = 0.5 * (
                result.projection_state.projected_response.dot(kpy_vec[component]) -
                trace_projection_component);
        }

        for (int row = 0; row < num_variance_terms; ++row) {
            for (int col = 0; col <= row; ++col) {
                result.ai(row, col) = result.ai(col, row) =
                    0.5 * kpy_vec[row].dot(pkpy_vec[col]);
            }
        }
        result.em = build_multitrait_continuous_em_information(
            result.varcom, num_id, num_traits);

        const VarianceUpdateResult update_result =
            find_valid_multitrait_continuous_variance_update(
                result.ai, result.em, result.score, result.varcom, num_traits);
        const double cc_par_val =
            std::sqrt(update_result.delta.squaredNorm() /
                      std::max(1.0, result.varcom.squaredNorm()));
        const double cc_gra_val = result.score.norm();
        const double cc_logL_val = std::isfinite(previous_neg2_reml)
            ? std::abs(result.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();
        const bool parameter_stable = cc_par_val < cc_par;
        const bool gradient_small = cc_gra_val < cc_gra;
        const bool objective_stable = cc_logL_val < cc_logL;

        spdlog::info("{} -2 REML: {}", log_prefix, result.neg2_reml);
        spdlog::info("{} Score: {}", log_prefix, format_vector_for_log(result.score));
        spdlog::info("{} Updated variances: {}", log_prefix,
            format_vector_for_log(update_result.updated_varcom));
        spdlog::info(
            "{} Convergence metrics: cc_par={}, cc_gra={}, cc_logL={}",
            log_prefix, cc_par_val, cc_gra_val, cc_logL_val);

        result.varcom = update_result.updated_varcom;
        result.converged = gradient_small && (parameter_stable || objective_stable);
        previous_neg2_reml = result.neg2_reml;
        if (result.converged) {
            break;
        }
    }

    result.projection_state = evaluate_multitrait_continuous_projection_state(
        grm_blocks,
        block_offsets,
        joint_design,
        stacked_response,
        result.varcom,
        num_id,
        num_traits);
    result.neg2_reml = result.projection_state.neg2_reml;
    result.score.setZero();
    result.ai.setZero();
    result.em = build_multitrait_continuous_em_information(
        result.varcom, num_id, num_traits);

    std::vector<VectorXd> final_kpy_vec(num_variance_terms);
    std::vector<VectorXd> final_pkpy_vec(num_variance_terms);
    for (int component = 0; component < num_variance_terms; ++component) {
        final_kpy_vec[component] = apply_multitrait_continuous_variance_component(
            grm_blocks,
            block_offsets,
            result.projection_state.projected_response,
            component,
            num_id,
            num_traits);
        final_pkpy_vec[component] = apply_projection(
            result.projection_state.inverse_result.inverse_matrix,
            result.projection_state.fixed_effect_projection.inverse_times_design,
            result.projection_state.fixed_effect_projection.design_cross_inverse,
            final_kpy_vec[component]);

        const MatrixXd design_cross_component =
            result.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
            apply_multitrait_continuous_variance_component(
                grm_blocks,
                block_offsets,
                result.projection_state.fixed_effect_projection.inverse_times_design,
                component,
                num_id,
                num_traits);
        const double trace_projection_component =
            result.projection_state.inverse_result.trace_vi_components(component) -
            (result.projection_state.fixed_effect_projection.design_cross_inverse
                .cwiseProduct(design_cross_component)).sum();
        result.score(component) = 0.5 * (
            result.projection_state.projected_response.dot(final_kpy_vec[component]) -
            trace_projection_component);
    }

    for (int row = 0; row < num_variance_terms; ++row) {
        for (int col = 0; col <= row; ++col) {
            result.ai(row, col) = result.ai(col, row) =
                0.5 * final_kpy_vec[row].dot(final_pkpy_vec[col]);
        }
    }

    return result;
}

VarianceUpdateResult find_positive_variance_update(const MatrixXd& ai_mat, const MatrixXd& em_mat,
        const VectorXd& fd_mat, const VectorXd& varcom) {
    VarianceUpdateResult update_result;

    for (int step = 0; step <= 100; ++step) {
        const double gamma = step * 0.01;
        const MatrixXd blended_information = (1 - gamma) * ai_mat + gamma * em_mat;
        update_result.delta = blended_information.inverse() * fd_mat;
        update_result.updated_varcom = varcom + update_result.delta;
        if (update_result.updated_varcom.minCoeff() > 0) {
            spdlog::info("EM weight value: {}", gamma);
            break;
        }
    }

    return update_result;
}

void write_variance_components_with_se(const std::string& out_file, const VectorXd& varcom, const VectorXd& se_varcom) {
    ofstream fout(out_file + ".var");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.var", out_file);
    }

    for(std::int64_t i = 0; i < varcom.size(); ++i){
        fout << varcom(i) << " " << se_varcom(i) << std::endl;
    }
}

void write_variance_components(const std::string& out_file, const VectorXd& varcom) {
    ofstream fout(out_file + ".var");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.var", out_file);
    }

    fout << varcom << std::endl;
}

MainModelIterationState evaluate_main_model_iteration(const VectorXd& grm_eigenvals,
        const MatrixXd& xmat_trans, const MatrixXd& y_trans, const VectorXd& varcom) {
    MainModelIterationState iteration_state;

    const VectorXd inverse_diagonal =
        (1.0 / (grm_eigenvals.array() * varcom(0) + varcom(1))).matrix();
    iteration_state.logL = -(inverse_diagonal.array().log()).sum();

    ProjectionCache fixed_effect_projection =
        build_projection_cache(inverse_diagonal, xmat_trans, true);
    iteration_state.logL += fixed_effect_projection.design_logdet;

    const VectorXd py = apply_projection(
        inverse_diagonal,
        fixed_effect_projection.inverse_times_design,
        fixed_effect_projection.design_cross_inverse,
        y_trans).col(0);
    iteration_state.logL += (y_trans.transpose() * py).sum();

    const VectorXd dpy = grm_eigenvals.cwiseProduct(py);
    const VectorXd pdpy = apply_projection(
        inverse_diagonal,
        fixed_effect_projection.inverse_times_design,
        fixed_effect_projection.design_cross_inverse,
        dpy);
    const VectorXd ppy = apply_projection(
        inverse_diagonal,
        fixed_effect_projection.inverse_times_design,
        fixed_effect_projection.design_cross_inverse,
        py);

    iteration_state.fd_mat(0) = (inverse_diagonal.cwiseProduct(grm_eigenvals)).sum();
    MatrixXd inverse_times_genetic_design =
        fixed_effect_projection.inverse_times_design.array().colwise() * grm_eigenvals.array();
    iteration_state.fd_mat(0) -= (
        fixed_effect_projection.design_cross_inverse.cwiseProduct(
            inverse_times_genetic_design.transpose() * fixed_effect_projection.inverse_times_design)).sum();
    iteration_state.fd_mat(0) -= (py.transpose() * dpy).sum();
    iteration_state.fd_mat(0) *= -0.5;

    iteration_state.fd_mat(1) = inverse_diagonal.sum() - (
        fixed_effect_projection.design_cross_inverse.cwiseProduct(
            fixed_effect_projection.inverse_times_design.transpose() *
            fixed_effect_projection.inverse_times_design)).sum();
    iteration_state.fd_mat(1) -= (py.transpose() * py).sum();
    iteration_state.fd_mat(1) *= -0.5;

    iteration_state.ai_mat(0, 0) = 0.5 * dpy.dot(pdpy);
    iteration_state.ai_mat(1, 1) = 0.5 * py.dot(ppy);
    iteration_state.ai_mat(0, 1) = 0.5 * dpy.dot(ppy);
    iteration_state.ai_mat(1, 0) = iteration_state.ai_mat(0, 1);

    return iteration_state;
}

BedScanContext prepare_bed_scan_context(GENO& geno, const std::vector<string>& sample_ids) {
    BedScanContext context;
    context.id_index_in_bed_vec = geno.find_fam_index(sample_ids);
    context.num_id_used = context.id_index_in_bed_vec.size();
    context.num_snp = geno.get_num_sid();
    context.num_id_in_bed = geno.get_num_iid();
    spdlog::info("The number of individuals and SNPs in the plink file are: {} {}",
        context.num_id_in_bed, context.num_snp);
    geno.validate_bed_size();
    return context;
}

std::vector<std::int64_t> sample_random_snp_indices(std::int64_t num_snp, std::int64_t num_random_snp) {
    return sample_unique_integers(0, num_snp, num_random_snp);
}

void load_centered_random_snp_panel(GENO& geno, const string& bed_file,
        const std::vector<std::int64_t>& snp_indices, const std::vector<std::int64_t>& id_index_in_bed_vec,
        MatrixXd& snp_mat, VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr) {
    geno.read_bed_by_snp_indices(
        bed_file + ".bed", snp_mat, freq_arr, missing_rate_arr, nobs_geno_arr, snp_indices, id_index_in_bed_vec);
    snp_mat.rowwise() -= 2 * freq_arr.transpose();
}

SnpPartition get_snp_partition(std::int64_t start_pos, std::int64_t end_pos, int npart_snp, int partition_index) {
    const std::int64_t base_partition_size = (end_pos - start_pos) / npart_snp;
    SnpPartition partition;
    partition.start_snp = partition_index * base_partition_size + start_pos;
    partition.num_snp_read = base_partition_size;
    if (partition_index == npart_snp - 1) {
        partition.num_snp_read += (end_pos - start_pos) % npart_snp;
    }
    return partition;
}

void require_random_snp_fraction(std::int64_t num_used, std::int64_t num_random_snp,
        double min_fraction, const char* error_template) {
    if (num_used * 1.0 / num_random_snp < min_fraction) {
        fatal_error("{}", fmt::format(fmt::runtime(error_template), num_used));
    }
}

ofstream open_output_stream(const string& file_path) {
    ofstream fout(file_path);
    if(!fout.is_open()){
        fatal_error("Fail to open output file: {}", file_path);
    }
    return fout;
}

void configure_run_cli(CLI::App& app, RunOptions& options) {
    app.description(R"(
    Quick Start:

      Test GxE:
        fastgxe --test-gxe --grm grm_file --bfile bed_file --data data_file --trait BMI --env-int age smoking:alcohol --out test_gxe

      Test SNP Main Effects:
        fastgxe --test-main --grm grm_file --bfile bed_file --data data_file --trait BMI --out test_main

      Test Binary-Trait SNP Main Effects:
        fastgxe --test-main-binary --grm grm_file --bfile bed_file --data data_file --trait disease --out test_main_binary

      Test Joint Binary + Continuous SNP Main Effects:
        fastgxe --test-main-binary-continuous --grm grm_file --bfile bed_file --data data_file --trait disease BMI --out test_main_binary_continuous

      Test Multitrait Continuous SNP Main Effects:
        fastgxe --test-main-multitrait-continuous --grm grm_file --bfile bed_file --data data_file --trait trait1 trait2 trait3 --out test_main_multitrait_continuous
    )");

    app.add_option("-p,--threads", options.threads,
        "Number of threads to use (default: 10).")
        ->default_val(10);

    auto* test_main_flag = app.add_flag("--test-main", options.test_main,
        "Enable testing of SNP main effects for continuous traits (mutually exclusive with --test-main-binary, --test-main-binary-continuous and --test-gxe).");
    auto* test_main_binary_flag = app.add_flag("--test-main-binary", options.test_main_binary,
        "Enable testing of SNP main effects for binary traits using the logistic mixed model (mutually exclusive with --test-main, --test-main-binary-continuous and --test-gxe).");
    auto* test_main_binary_continuous_flag = app.add_flag(
        "--test-main-binary-continuous",
        options.test_main_binary_continuous,
        "Enable the joint binary-plus-continuous SNP main-effect model. The first --trait is binary and the second --trait is continuous (mutually exclusive with --test-main, --test-main-binary and --test-gxe).");
    auto* test_main_multitrait_continuous_flag = app.add_flag(
        "--test-main-multitrait-continuous",
        options.test_main_multitrait_continuous,
        "Enable the multitrait continuous SNP main-effect model. Use two or more continuous --trait values (mutually exclusive with --test-main, --test-main-binary, --test-main-binary-continuous and --test-gxe).");
    auto* test_gxe_flag = app.add_flag("--test-gxe", options.test_gxe,
        "Enable testing of SNP-environment interactions (mutually exclusive with --test-main, --test-main-binary, --test-main-binary-continuous and --test-main-multitrait-continuous).");
    test_main_flag->excludes(test_gxe_flag);
    test_main_flag->excludes(test_main_binary_flag);
    test_main_flag->excludes(test_main_binary_continuous_flag);
    test_main_flag->excludes(test_main_multitrait_continuous_flag);
    test_main_binary_flag->excludes(test_main_flag);
    test_main_binary_flag->excludes(test_gxe_flag);
    test_main_binary_flag->excludes(test_main_binary_continuous_flag);
    test_main_binary_flag->excludes(test_main_multitrait_continuous_flag);
    test_main_binary_continuous_flag->excludes(test_main_flag);
    test_main_binary_continuous_flag->excludes(test_main_binary_flag);
    test_main_binary_continuous_flag->excludes(test_gxe_flag);
    test_main_binary_continuous_flag->excludes(test_main_multitrait_continuous_flag);
    test_main_multitrait_continuous_flag->excludes(test_main_flag);
    test_main_multitrait_continuous_flag->excludes(test_main_binary_flag);
    test_main_multitrait_continuous_flag->excludes(test_main_binary_continuous_flag);
    test_main_multitrait_continuous_flag->excludes(test_gxe_flag);
    test_gxe_flag->excludes(test_main_flag);
    test_gxe_flag->excludes(test_main_binary_flag);
    test_gxe_flag->excludes(test_main_binary_continuous_flag);
    test_gxe_flag->excludes(test_main_multitrait_continuous_flag);

    app.add_flag("--no-noisebye", [&options](int count) {
        if (count > 0) options.no_noisebye = true;
    }, "Disable noise-by-environment interaction terms.\n"
       "  - Noise-by-environment interactions are ENABLED by default.\n"
       "  - Use --no-noisebye to turn it OFF.");

    app.add_option("--data", options.data_file, "Path to input data file (required).")->required();
    app.add_option("--trait", options.trait_vec,
        "Trait(s) to analyze.\n"
        "  - Use 1 trait for --test-main, --test-main-binary and --test-gxe.\n"
        "  - Use 2 traits for --test-main-binary-continuous: first binary, second continuous.\n"
        "  - Use 2 or more traits for --test-main-multitrait-continuous: all continuous.")
        ->expected(-1)->required();
    app.add_option("--env-int", options.bye_vec,
        "List of interacting environmental covariates.\n"
        "  - Supports multiple values (e.g., --env-int age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --env-int age:BMI).\n"
        "  - Expands to include all covariates in the specified range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);
    app.add_option("--out", options.out_file, "Path to output file (required).")->required();
    app.add_option("--grm", options.agrm_file, "Path to genetic relationship matrix (GRM) file (required).")->required();
    app.add_option("--bfile", options.bed_file, "Path to binary PLINK BED file (optional).");

    app.add_option("--covar", options.covariate_vec,
        "List of covariates for the analysis (required).\n"
        "  - Supports multiple values (e.g., --covar age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --covar age:BMI).\n"
        "  - Expands to include all covariates in the range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);
    app.add_option("--class", options.class_vec,
        "List of categorical class variables (required).\n"
        "  - Supports multiple values (e.g., --class gender ethnicity region).\n"
        "  - Use ':' to specify a range (e.g., --class gender:region).\n"
        "  - Expands to include all variables in the range.\n"
        "  - Example order: {gender, ethnicity, region, education, income}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);

    app.add_flag("--no-standardize-env", [&options](int count) {
        if (count > 0) options.standardize_env = false;
    }, "Disable standardization of interacting environmental covariates.\n"
       "  - Standardization is ENABLED by default (mean 0, std 1).\n"
       "  - Use --no-standardize-env to turn it OFF.");

    app.add_option("--missing-data", options.missing_in_data_vec,
        "List of missing value indicators for phenotype/covariates.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-data . -999 \"?\").")
        ->expected(-1);

    app.add_option("--maxiter", options.maxiter, "Maximum number of optimization iterations (default: 200).")
        ->default_val(200);
    app.add_option("--cc-par", options.cc_par, "Convergence threshold for parameter updates (default: 1e-7).")
        ->default_val(1e-7);
    app.add_option("--cc-gra", options.cc_gra, "Convergence threshold for gradient norm (default: 1e-6).")
        ->default_val(1e-6);
    app.add_option("--cc-logL", options.cc_logL, "Convergence threshold for log-likelihood change (default: 5e-5).")
        ->default_val(5e-5);
    app.add_option("--snp-range", options.snp_range_vec, "Specify start and end SNPs (expects 2 values).")->expected(2);
    app.add_option("--npart-SNP", options.npart_snp, "Number of SNP partitions (default: 1000).")->default_val(1000);
    app.add_option("--speed", options.speed, "Computation speed level (default: 2, higher is faster).")->default_val(2);
    app.add_option("--p-approx-cut", options.p_approx_cut, "P-value approximation threshold (default: 1e-3).")->default_val(1e-3);
    app.add_option("--p-cut", options.p_cut, "P-value significance cutoff for output (default: 2).")->default_val(2);
    app.add_option("--maf", options.maf_cut, "Minor allele frequency (MAF) cutoff (default: 0.01).")->default_val(0.01);
    app.add_option("--missing-rate", options.missing_rate_cut, "Missing genotype rate cutoff (default: 0.05).")->default_val(0.05);
    app.add_option("--num-random-snp", options.num_random_snp, "Number of randomly selected SNPs (default: 2000).")->default_val(2000);
    app.add_option("--split-task", options.split_task_vec,
        "Partition the task for parallel execution (expects 2 values: total parts, current part).")
        ->expected(2);
}

void validate_run_options(RunOptions& options) {
    const int num_analysis_modes =
        static_cast<int>(options.test_main) +
        static_cast<int>(options.test_main_binary) +
        static_cast<int>(options.test_main_binary_continuous) +
        static_cast<int>(options.test_main_multitrait_continuous) +
        static_cast<int>(options.test_gxe);
    if (num_analysis_modes != 1) {
        fatal_error("Exactly one of --test-main, --test-main-binary, --test-main-binary-continuous, --test-main-multitrait-continuous or --test-gxe must be specified.");
    }
    if (options.test_gxe && options.bye_vec.empty()) {
        fatal_error("--env-int is required when --test-gxe is specified.");
    }
    const bool requires_exactly_two_traits =
        options.test_main_binary_continuous;
    const bool requires_at_least_two_traits =
        options.test_main_multitrait_continuous;
    const std::size_t expected_num_traits =
        requires_exactly_two_traits ? 2 : 1;
    if ((requires_exactly_two_traits && options.trait_vec.size() != expected_num_traits) ||
            (requires_at_least_two_traits && options.trait_vec.size() < 2) ||
            (!requires_exactly_two_traits && !requires_at_least_two_traits && options.trait_vec.size() != 1)) {
        if (options.test_main_binary_continuous) {
            fatal_error("--test-main-binary-continuous requires exactly two traits: first binary, second continuous.");
        } else if (options.test_main_multitrait_continuous) {
            fatal_error("--test-main-multitrait-continuous requires at least two continuous traits.");
        } else {
            fatal_error("This analysis mode requires exactly one trait.");
        }
    }
    if (options.threads <= 0) {
        options.threads = 10;
    }
}

void log_run_options(const RunOptions& options) {
    spdlog::info("=== Parsed Arguments ===");
    spdlog::info("Threads: {}", options.threads);
    if (options.test_main) {
        spdlog::info("Analysis Mode: Continuous SNP main effects");
    } else if (options.test_main_binary) {
        spdlog::info("Analysis Mode: Binary SNP main effects");
    } else if (options.test_main_binary_continuous) {
        spdlog::info("Analysis Mode: Joint binary + continuous SNP main effects");
    } else if (options.test_main_multitrait_continuous) {
        spdlog::info("Analysis Mode: Multitrait continuous SNP main effects");
    } else if (options.test_gxe) {
        spdlog::info("Analysis Mode: SNP-environment interactions");
    }
    spdlog::info("Input Data File: {}", options.data_file);
    spdlog::info("Output File: {}", options.out_file);
    spdlog::info("GRM File: {}", options.agrm_file);
    spdlog::info("BED File: {}", options.bed_file.empty() ? "Not Provided" : options.bed_file);
    spdlog::info("Trait(s): {}", join_string(options.trait_vec, ", "));
    spdlog::info("Interacting environmental covariates: {}", join_or_none(options.bye_vec));
    spdlog::info("Standardization of Interacting environments: {}", options.standardize_env ? "ENABLED" : "DISABLED");
    spdlog::info("Covariates: {}", join_or_none(options.covariate_vec));
    spdlog::info("Class Variables: {}", join_or_none(options.class_vec));
    spdlog::info("Missing Data Indicators: {}", join_string(options.missing_in_data_vec, ", "));
    spdlog::info("Max Iterations: {}", options.maxiter);
    spdlog::info("CC Par: {}, CC Gra: {}, CC LogL: {}", options.cc_par, options.cc_gra, options.cc_logL);
    if (options.snp_range_vec.size() == 2) {
        spdlog::info("SNP Range: {} to {}", options.snp_range_vec[0], options.snp_range_vec[1]);
    }
    spdlog::info("SNP Partitions: {}", options.npart_snp);
    spdlog::info("P-value Approx Cutoff: {}", options.p_approx_cut);
    spdlog::info("P-value Output Cutoff: {}", options.p_cut);
    spdlog::info("Minor Allele Frequency Cutoff: {}", options.maf_cut);
    spdlog::info("Missing Rate Cutoff: {}", options.missing_rate_cut);
    spdlog::info("Number of Random SNPs: {}", options.num_random_snp);
    if(!options.split_task_vec.empty()) {
        spdlog::info("Task Partitioning: {} parts, running part {}", options.split_task_vec[0], options.split_task_vec[1]);
    }
    spdlog::info("========================");
}

ResolvedInputColumns resolve_input_columns(RunOptions& options) {
    const std::vector<std::string> head_vec = read_header_columns(options.data_file);
    std::vector<std::string> missing_names;
    ResolvedInputColumns resolved;

    resolved.trait_index_vec = find_index(head_vec, options.trait_vec, missing_names);
    if (!missing_names.empty()) {
        fatal_error("Trait names not found in the header: {}", join_string(missing_names));
    }
    missing_names.clear();

    options.covariate_vec = expand_variable_ranges(options.covariate_vec, head_vec);
    if(!options.covariate_vec.empty()){
        spdlog::info("Number of covariates: {}, Covariates: {}",
             options.covariate_vec.size(), join_string(options.covariate_vec, ", "));
    }
    resolved.covariate_index_vec = find_index(head_vec, options.covariate_vec, missing_names);
    missing_names.clear();

    options.class_vec = expand_variable_ranges(options.class_vec, head_vec);
    if(!options.class_vec.empty()){
        spdlog::info("Number of class variables: {}, class variables: {}",
             options.class_vec.size(), join_string(options.class_vec, ", "));
    }
    resolved.class_index_vec = find_index(head_vec, options.class_vec, missing_names);
    missing_names.clear();

    options.bye_vec = expand_variable_ranges(options.bye_vec, head_vec);
    spdlog::info("Number of interacting environmental covariates: {}, Interacting environmental covariates: {}",
         options.bye_vec.size(), join_string(options.bye_vec, ", "));
    resolved.bye_index_vec = find_index(head_vec, options.bye_vec, missing_names);
    if(!missing_names.empty()){
        fatal_error("Interacting environments not found in the header: {}", join_string(missing_names));
    }
    missing_names.clear();

    if (!resolved.bye_index_vec.empty()) {
        const std::vector<std::int64_t> covar_overlap = find_index(options.covariate_vec, options.bye_vec, missing_names);
        const std::vector<std::int64_t> class_overlap = find_index(options.class_vec, options.bye_vec, missing_names);
        if (!covar_overlap.empty() || !class_overlap.empty()) {
            fatal_error("Interacting environments should not be included in --covar or --class. They are treated as covariates by default.");
        }
    }

    return resolved;
}

}  // namespace

vector<string> fastGxE::extract_retained_sample_ids(const PHEN& phenoA) const {
    std::vector<string> sample_ids;
    phenoA.get_given_column(0, sample_ids);

    std::unordered_set<string> unique_ids;
    unique_ids.reserve(sample_ids.size());
    for (const auto& sample_id : sample_ids) {
        if (!unique_ids.insert(sample_id).second) {
            fatal_error("Duplicated sample IDs exist in the data file!");
        }
    }

    if (sample_ids.empty()) {
        fatal_error("No samples remain after aligning the phenotype data with the GRM IDs.");
    }

    return sample_ids;
}

vector<fastGxE::GroupRecord> fastGxE::load_and_sort_group_records(const string& agrm_file,
        const vector<string>& retained_sample_ids, std::int64_t expected_num_samples) const {
    spdlog::info("Read the group files: {}.grm.group", agrm_file);
    std::unordered_set<string> retained_id_set(retained_sample_ids.begin(), retained_sample_ids.end());

    ifstream fin(agrm_file + ".grm.group");
    if(!fin.is_open()){
        fatal_error("Fail to open the {}.grm.group", agrm_file);
    }

    string line;
    std::vector<GroupRecord> group_records;
    group_records.reserve(expected_num_samples);
    std::unordered_map<std::int64_t, std::int64_t> group_size_map;
    group_size_map.reserve(expected_num_samples);

    while(std::getline(fin, line)){
        process_line(line);
        if(line.empty()){
            continue;
        }

        vector<string> vec = split_string(line);
        if(vec.size() < 4){
            fatal_error("Malformed record in {}.grm.group: {}", agrm_file, line);
        }

        if(retained_id_set.find(vec[1]) == retained_id_set.end()){
            continue;
        }

        GroupRecord record;
        record.grm_index = std::atoll(vec[0].c_str());
        record.group_id = std::atoll(vec[2].c_str());
        group_records.push_back(record);
        ++group_size_map[record.group_id];
    }

    if(group_records.size() != expected_num_samples){
        fatal_error("The number of retained records in {}.grm.group ({}) does not match the retained phenotype sample count ({}).",
                      agrm_file, group_records.size(), expected_num_samples);
    }

    spdlog::info("Re-calculate the group size");
    for(auto& record : group_records){
        record.group_size = group_size_map[record.group_id];
    }

    spdlog::info("Sort by group size and group index");
    sort(group_records.begin(), group_records.end(), [](const GroupRecord& lhs, const GroupRecord& rhs){
        if(lhs.group_size != rhs.group_size){
            return lhs.group_size < rhs.group_size;
        }
        if(lhs.group_id != rhs.group_id){
            return lhs.group_id < rhs.group_id;
        }
        return lhs.grm_index < rhs.grm_index;
    });

    return group_records;
}

vector<string> fastGxE::build_reordered_sample_ids(const vector<GroupRecord>& group_records,
        const vector<string>& grm_sample_ids) const {
    vector<string> reordered_sample_ids;
    reordered_sample_ids.reserve(group_records.size());
    for(const auto& record : group_records){
        if(record.grm_index <= 0 || record.grm_index > static_cast<std::int64_t>(grm_sample_ids.size())){
            fatal_error("GRM index {} in group file is out of range.", record.grm_index);
        }
        reordered_sample_ids.push_back(grm_sample_ids[record.grm_index - 1]);
    }
    return reordered_sample_ids;
}

void fastGxE::log_removed_design_columns(const vector<std::int64_t>& removed_column_indices) const {
    if (removed_column_indices.empty()) {
        return;
    }

    spdlog::info("Dependent columns of xmat are removed: ");
    std::string removed_columns;
    for (const auto& col : removed_column_indices) {
        removed_columns += std::to_string(col) + " ";
    }
    spdlog::info("{}", removed_columns);
}

void fastGxE::build_fixed_and_environment_matrices(PHEN& phenoA,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait) {
    spdlog::info("Design matrix for fixed effects");
    vector<std::int64_t> index_fixed_effect_vec;
    vector<string> fixed_effect_name_vec;
    m_xmat = phenoA.build_fixed_effect_design_matrix(
        covariate_arr, class_arr, index_fixed_effect_vec, fixed_effect_name_vec, trait, m_y);

    vector<std::int64_t> removed_column_indices;
    phenoA.remove_dependent_design_columns(
        m_xmat, index_fixed_effect_vec, fixed_effect_name_vec, removed_column_indices);
    log_removed_design_columns(removed_column_indices);

    spdlog::info("Get interaction environment covariates");
    if(!bye_arr.empty()){
        phenoA.get_columns_as_matrix(bye_arr, m_bye_mat);
    }
}

fastGxE::GroupLayout fastGxE::build_group_layout(const vector<GroupRecord>& group_records) const {
    spdlog::info("Store sample index within subgroup");

    std::vector<std::int64_t> singleton_indices;
    std::vector<vector<std::int64_t>> dense_group_indices;
    singleton_indices.reserve(group_records.size());
    dense_group_indices.reserve(group_records.size());

    std::int64_t current_group_id = std::numeric_limits<std::int64_t>::min();
    for(const auto& record : group_records){
        if(record.group_size == 1){
            singleton_indices.push_back(record.grm_index);
            continue;
        }

        if(dense_group_indices.empty() || record.group_id != current_group_id){
            dense_group_indices.push_back({record.grm_index});
            current_group_id = record.group_id;
        }else{
            dense_group_indices.back().push_back(record.grm_index);
        }
    }

    GroupLayout layout;
    layout.sample_location_map.reserve(group_records.size());
    layout.grm_blocks.reserve(dense_group_indices.size() + 1);
    layout.grm_blocks.emplace_back(MatrixXd::Zero(singleton_indices.size(), 1));

    for(std::int64_t i = 0; i < static_cast<std::int64_t>(singleton_indices.size()); ++i){
        layout.sample_location_map.emplace(singleton_indices[i], SampleBlockLocation{-(i + 1), i});
    }

    for(std::int64_t block_idx = 0; block_idx < static_cast<std::int64_t>(dense_group_indices.size()); ++block_idx){
        const auto& group_indices = dense_group_indices[block_idx];
        layout.grm_blocks.emplace_back(MatrixXd::Zero(group_indices.size(), group_indices.size()));
        const std::int64_t block_id = block_idx + 1;
        for(std::int64_t row_idx = 0; row_idx < static_cast<std::int64_t>(group_indices.size()); ++row_idx){
            layout.sample_location_map.emplace(group_indices[row_idx], SampleBlockLocation{block_id, row_idx});
        }
    }

    return layout;
}

void fastGxE::load_grouped_grm_blocks(const string& agrm_file,
        const std::unordered_map<std::int64_t, SampleBlockLocation>& sample_location_map,
        vector<MatrixXd>& grm_blocks) const {
    spdlog::info("Read the GRMs");

    std::ifstream fin(agrm_file + ".grm.sp.bin", std::ios::binary);
    if(!fin.is_open()){
        fatal_error("Fail to open: {}.grm.sp.bin", agrm_file);
    }

    double value;
    std::int64_t index0, index1;

    while(fin.read((char*)&index0, sizeof(std::int64_t)) &&
            fin.read((char*)&index1, sizeof(std::int64_t)) &&
            fin.read((char*)&value, sizeof(double))){
        const auto loc0_it = sample_location_map.find(index0 + 1);
        if(loc0_it == sample_location_map.end()){
            continue;
        }

        const auto loc1_it = sample_location_map.find(index1 + 1);
        if(loc1_it == sample_location_map.end()){
            continue;
        }

        const auto& loc0 = loc0_it->second;
        const auto& loc1 = loc1_it->second;
        if(loc0.block_id != loc1.block_id){
            continue;
        }

        if(loc0.block_id > 0){
            grm_blocks[loc0.block_id](loc0.row_index, loc1.row_index) = value;
            grm_blocks[loc0.block_id](loc1.row_index, loc0.row_index) = value;
        }else if(loc0.row_index == loc1.row_index){
            grm_blocks[0](loc0.row_index, 0) = value;
        }
    }
}

fastGxE::fastGxE(){

}

fastGxE::~fastGxE() {
}

void fastGxE::initialize_analysis_metadata(const string& out_file,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait) {
    spdlog::info("Start analysis...");
    m_out_file = out_file;

    m_num_trait = trait.size();
    const std::int64_t num_covariate = covariate_arr.size();
    const std::int64_t num_class = class_arr.size();
    m_num_bye = bye_arr.size();
    spdlog::info("The number of analyzed traits: {}", m_num_trait);
    spdlog::info("The number of covariates: {}", num_covariate);
    spdlog::info("The number of classes: {}", num_class);
    spdlog::info("The number of interaction environments: {}", m_num_bye);
}

void fastGxE::load_filtered_phenotype(PHEN& phenoA, const string& data_file,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait,
        const vector<string>& grm_sample_ids, const vector<string>& missing_in_data_vec) {
    spdlog::info("Read data file");
    const vector<std::int64_t> col_used_vec = merge_unique_vectors({covariate_arr, class_arr, bye_arr, trait});
    phenoA.read_and_filter(data_file, {0}, {grm_sample_ids}, col_used_vec, missing_in_data_vec);

    m_id_in_data_vec = extract_retained_sample_ids(phenoA);
    m_num_id = m_id_in_data_vec.size();
    spdlog::info("The number of used iids in data file: {}", m_num_id);
}

vector<fastGxE::GroupRecord> fastGxE::reorder_phenotype_by_group(PHEN& phenoA, const string& agrm_file,
        const vector<string>& grm_sample_ids) {
    const vector<GroupRecord> group_records =
        load_and_sort_group_records(agrm_file, m_id_in_data_vec, m_num_id);

    spdlog::info("Update the dataframe using the current sample id orders");
    phenoA.reorder_records_by_sample_id_order(build_reordered_sample_ids(group_records, grm_sample_ids));
    m_id_in_data_vec = extract_retained_sample_ids(phenoA);
    return group_records;
}

/**
 * @Description: Read the data file and GRM file
 */
int fastGxE::pre_data(const string& out_file, const string& data_file, const string& agrm_file,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr,  const vector<std::int64_t>& trait,
        const vector<string>& missing_in_data_vec){
    // Stage 1: collect metadata and load the phenotype view restricted to the
    // samples that also exist in the GRM.
    initialize_analysis_metadata(out_file, covariate_arr, class_arr, bye_arr, trait);

    // id in grm
    spdlog::info("Read iids in GRM");
    ProcessGRM ProcessGRMA;
    const vector<string> id_in_gmat_vec = ProcessGRMA.read_grm_id(agrm_file + ".grm");
    const std::int64_t num_id_in_grm = id_in_gmat_vec.size();
    spdlog::info("The sample size in the GRM is: {}", num_id_in_grm);

    PHEN phenoA;
    load_filtered_phenotype(
        phenoA, data_file, covariate_arr, class_arr, bye_arr, trait, id_in_gmat_vec, missing_in_data_vec);

    // Stage 2: align sample order between the phenotype table and the grouped
    // GRM representation before building any matrices.
    //
    // Keep the phenotype view and GRM view in the same sample order before any
    // downstream matrix construction, otherwise every later multiplication will
    // silently mix incompatible row/column orders.
    const vector<GroupRecord> group_records = reorder_phenotype_by_group(phenoA, agrm_file, id_in_gmat_vec);

    // Stage 3: build the fixed-effect design and the grouped GRM blocks that
    // every downstream variance-component and SNP-scan routine reuses.
    build_fixed_and_environment_matrices(phenoA, covariate_arr, class_arr, bye_arr, trait);

    GroupLayout group_layout = build_group_layout(group_records);
    m_grm_mat_group_vec = std::move(group_layout.grm_blocks);
    load_grouped_grm_blocks(agrm_file, group_layout.sample_location_map, m_grm_mat_group_vec);


    return 0;
}








/**
 * @Description: Sparse GRM or perform Eigen Decompostion of GRM
 * @param {bool} use_eigen
 */
void fastGxE::process_grm(bool use_eigen){
    // Block 0 stores singleton samples as an N x 1 vector, while later blocks are
    // dense within-group GRMs. Keep their global row offsets cached once here.
    m_grm_index_vec = build_grm_block_offsets(m_grm_mat_group_vec);

    std::vector<GrmTriplet> triplet_list;
    triplet_list.reserve(count_grm_triplets(m_grm_mat_group_vec));

    if(!use_eigen){
        spdlog::info("Prepare sparse GRMs");
        // The main-effect workflow keeps the sparse GRM directly because later
        // scans use block-wise sparse products rather than a rotated basis.
        for(std::size_t block_idx = 0; block_idx < m_grm_mat_group_vec.size(); ++block_idx){
            const MatrixXd& block = m_grm_mat_group_vec[block_idx];
            const std::int64_t start_index = m_grm_index_vec[block_idx];
            if(block_idx == 0){
                append_singleton_block_triplets(block, start_index, triplet_list, false);
            }else{
                append_dense_block_triplets(block, start_index, triplet_list);
            }
        }
        m_grm_mat.resize(m_num_id, m_num_id);
        m_grm_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
    }else{
        spdlog::info("Eigen Decompostion of GRM");
        // The variance-component solver for SNP main effects works in the GRM
        // eigen basis, so precompute both eigenvalues and the sparse rotation.
        m_grm_eigenvals.resize(m_num_id);
        for(std::size_t block_idx = 0; block_idx < m_grm_mat_group_vec.size(); ++block_idx){
            const MatrixXd& block = m_grm_mat_group_vec[block_idx];
            const std::int64_t start_index = m_grm_index_vec[block_idx];
            if(block_idx == 0){
                load_singleton_block_eigenvalues(block, start_index, m_grm_eigenvals);
                append_singleton_block_triplets(block, start_index, triplet_list, true);
            }else{
                append_dense_block_eigendecomposition(block, start_index, m_grm_eigenvals, triplet_list);
            }
        }

        m_grm_mat_group_vec.shrink_to_fit();
        m_grm_eigenvecs.resize(m_num_id, m_num_id);
        m_grm_eigenvecs.setFromTriplets(triplet_list.begin(), triplet_list.end());

        m_y_trans = m_grm_eigenvecs.transpose() * m_y;
        m_xmat_trans = m_grm_eigenvecs.transpose() * m_xmat;
    }
}



/**
 * @Description: Scale the interaction environment covariates and add to xmat
 */
void fastGxE::pre_data_GxE(bool standardize_env, bool phen_correct){
    spdlog::info("Scale the interaction environment covariates");
    require_full_column_rank(
        m_bye_mat,
        m_num_bye,
        "There are dependent columns in the interaction environment covariates!");

    if(standardize_env){
        standardize_environment_matrix(m_bye_mat);
    }

    const std::int64_t nrows = m_xmat.rows();
    append_columns(m_xmat, m_bye_mat);
    const ColPivHouseholderQR<MatrixXd> design_qr = require_full_column_rank(
        m_xmat,
        m_xmat.cols(),
        "There are dependent columns between interaction environment covariates and other covariates!");

    if(phen_correct){
        // After projecting out fixed effects and environments, downstream GxE fitting
        // only needs an intercept in the fixed-effect design.
        m_y = m_y - m_xmat * design_qr.solve(m_y);
        m_xmat = MatrixXd::Ones(nrows, 1);
    }
}


/**
 * @Description: Improved structLMM model; Include the GRM and GxE rm 
 * @param {VectorXd&} init_varcom
 * @param {int} maxiter0
 * @param {double} cc_par0
 * @param {double} cc_gra0
 * @param {double} cc_logL0
 */
VectorXd fastGxE::varcom_GxE(VectorXd& init_varcom, bool no_noisebye, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0){
    spdlog::info("Estimate variances...");
    (void)cc_logL0;

    const std::int64_t num_cov = no_noisebye ? 3 : 4;
    VectorXd varcom = initialize_variance_components(init_varcom, num_cov);
    VectorXd ai_mat_inv_diag = VectorXd::Zero(num_cov);

    spdlog::info("Initial variances are: {}", format_vector_for_log(varcom));

    spdlog::info("Calculate GxE relationship matrix");
    // Precompute the environment-weighted GRM once, then reuse it across every
    // AI/EM iteration instead of rebuilding it from scratch each time.
    const GxeRelationshipData gxe_relationship = build_gxe_relationship_data(
        m_grm_mat_group_vec, m_grm_index_vec, m_bye_mat, m_num_bye, m_num_id);

    spdlog::info("Start the iteration");

    bool isCC = true;
    const VectorXd nxe_vec = build_noise_by_environment_vector(m_bye_mat);
    const SparseMatrix<double> nxe_diag = build_sparse_diagonal_matrix(nxe_vec);

    for(int iter_count = 0; iter_count < maxiter0; iter_count++){
        spdlog::info("Iteration {}", iter_count + 1);

        // Rebuild V^{-1} under the current variance components, then evaluate
        // the REML objective, score vector and AI matrix for the next update.
        const InverseCovarianceResult inverse_result = build_inverse_covariance_matrix(
            m_grm_mat_group_vec,
            gxe_relationship.group_blocks,
            m_grm_index_vec,
            nxe_vec,
            varcom,
            no_noisebye,
            m_num_id);
        const SparseMatrix<double>& vmat = inverse_result.inverse_matrix;

        // Evaluate the REML objective, gradient and AI matrix under the current V^{-1}.
        vector <MatrixXd> KPy_vec(num_cov);
        VectorXd fd_mat(num_cov);
        
        const ProjectionCache fixed_effect_projection =
            build_projection_cache(vmat, m_xmat, true);
        MatrixXd design_cross_genetic_design =
            fixed_effect_projection.inverse_times_design.transpose() * m_grm_mat *
            fixed_effect_projection.inverse_times_design;
        MatrixXd projected_response = apply_projection(
            vmat,
            fixed_effect_projection.inverse_times_design,
            fixed_effect_projection.design_cross_inverse,
            m_y);
        const double logL = inverse_result.logdet + fixed_effect_projection.design_logdet +
            (m_y.transpose() * projected_response).sum();
        spdlog::info("-2logL: {}", logL);
        MatrixXd KPy = m_grm_mat * projected_response;
        double tr_ViK = (vmat.cwiseProduct(m_grm_mat)).sum();
        double tr_design_cross_genetic = (
            fixed_effect_projection.design_cross_inverse.cwiseProduct(design_cross_genetic_design)).sum();
        fd_mat(0) = tr_ViK - tr_design_cross_genetic - (projected_response.transpose() * KPy).sum();
        KPy_vec[0] = KPy;

        double tr_ViKe = (vmat.cwiseProduct(gxe_relationship.sparse_matrix)).sum();
        MatrixXd design_cross_gxe_design =
            fixed_effect_projection.inverse_times_design.transpose() *
            gxe_relationship.sparse_matrix *
            fixed_effect_projection.inverse_times_design;
        double tr_design_cross_gxe = (
            fixed_effect_projection.design_cross_inverse.cwiseProduct(design_cross_gxe_design)).sum();
        MatrixXd KePy = gxe_relationship.sparse_matrix * projected_response;
        fd_mat(1) = tr_ViKe - tr_design_cross_gxe - (projected_response.transpose() * KePy).sum();
        KPy_vec[1] = KePy;

        if(!no_noisebye){
            double tr_ViNxe = (vmat.cwiseProduct(nxe_diag)).sum();
            MatrixXd design_cross_noise_design =
                fixed_effect_projection.inverse_times_design.transpose() * nxe_diag *
                fixed_effect_projection.inverse_times_design;
            double tr_design_cross_noise = (
                fixed_effect_projection.design_cross_inverse.cwiseProduct(design_cross_noise_design)).sum();
            MatrixXd NxePy = nxe_diag * projected_response;
            fd_mat(2) = tr_ViNxe - tr_design_cross_noise - (projected_response.transpose() * NxePy).sum();
            KPy_vec[2] = NxePy;
        }



        fd_mat(num_cov - 1) = vmat.diagonal().sum() - (
            fixed_effect_projection.design_cross_inverse.cwiseProduct(
                fixed_effect_projection.inverse_times_design.transpose() *
                fixed_effect_projection.inverse_times_design)).sum() -
            (projected_response.transpose() * projected_response).sum();
        fd_mat *= -0.5;
        KPy_vec[num_cov - 1] = projected_response;

        vector <MatrixXd> PKPy_vec(num_cov);
        for(std::int64_t m = 0; m < num_cov; m++){
            const MatrixXd& KPy_tmp = KPy_vec[m];
            PKPy_vec[m] = apply_projection(
                vmat,
                fixed_effect_projection.inverse_times_design,
                fixed_effect_projection.design_cross_inverse,
                KPy_tmp);
        }

        // AI
        MatrixXd ai_mat(num_cov, num_cov);
        for(std::int64_t m = 0; m < num_cov; m++){
            ai_mat(m, m) = (KPy_vec[m].transpose() * PKPy_vec[m]).sum();
            for(std::int64_t n = 0; n < m; n++){
                ai_mat(m, n) = ai_mat(n, m) = (KPy_vec[m].transpose() * PKPy_vec[n]).sum();
            }
        }
        ai_mat *= 0.5;
        ai_mat_inv_diag = ai_mat.inverse().diagonal();

        // EM information diagonal: n/(2σ⁴) per variance parameter; a valid lower bound when AI is ill-conditioned.
        MatrixXd em_mat = MatrixXd::Zero(num_cov, num_cov);
        em_mat.diagonal() = m_num_id / (2 * varcom.array() * varcom.array());
        
        // Update
        const VarianceUpdateResult update_result = find_positive_variance_update(ai_mat, em_mat, fd_mat, varcom);
        spdlog::info("Updated variances: {}", format_vector_for_log(update_result.updated_varcom));

        varcom = update_result.updated_varcom;
        
        double cc_par_val = (update_result.delta.cwiseProduct(update_result.delta)).sum() / (varcom.cwiseProduct(varcom)).sum();
        cc_par_val = sqrt(cc_par_val);
        double cc_gra_val = (fd_mat.cwiseProduct(fd_mat)).sum();
        cc_gra_val = sqrt(cc_gra_val);

        if (cc_par_val < cc_par0 || cc_gra_val < cc_gra0){
            isCC = true;
            break;
        }else{
            isCC = false;
        }
    }

   if (isCC)
        spdlog::info("Variances Converged");
    else
        spdlog::info("Variances Not Converged");
    
    const InverseCovarianceResult final_inverse = build_inverse_covariance_matrix(
        m_grm_mat_group_vec,
        gxe_relationship.group_blocks,
        m_grm_index_vec,
        nxe_vec,
        varcom,
        no_noisebye,
        m_num_id);
    this->m_V0_logdet = final_inverse.logdet;
    this->m_Vi0 = final_inverse.inverse_matrix;

    // logL
    const ProjectionCache fixed_effect_projection =
        build_projection_cache(m_Vi0, m_xmat, true);
    double logL = m_V0_logdet + fixed_effect_projection.design_logdet;
    MatrixXd projected_response = apply_projection(
        this->m_Vi0,
        fixed_effect_projection.inverse_times_design,
        fixed_effect_projection.design_cross_inverse,
        m_y);
    logL += (m_y.transpose() * projected_response).sum();
    spdlog::info("-2logL in null model: {}", logL);
    m_logL_null = logL;
    
    m_varcom_null = varcom;

    VectorXd se_varcom = ai_mat_inv_diag.cwiseSqrt();
    write_variance_components_with_se(m_out_file, varcom, se_varcom);

    return m_varcom_null;
}




/**
 * @Description: Estimate variances for model with only one random effect
 * @param {VectorXd&} init
 * @param {int} maxiter0
 * @param {double} cc_par0
 * @param {double} cc_gra0
 * @param {double} cc_logL0
 */
VectorXd fastGxE::varcom_main(VectorXd& init, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0){

    spdlog::info("Estimate variances using eigen decompostion...");

    // Preserve the original behavior: the main-effect model always starts from
    // a unit initialization rather than reusing any caller-supplied values.
    init = VectorXd::Ones(2);
    VectorXd varcom = init;

    spdlog::info("Initial variances are: {}", format_vector_for_log(varcom));

    double logL = -1e100;
    bool isCC = true;

    for (int i = 0; i < maxiter0; i++) {
        spdlog::info("Iteration {}", i + 1);

        // In the eigen basis the null model reduces to diagonal algebra, so one
        // helper can produce the log-likelihood, score and AI terms together.
        const MainModelIterationState iteration_state = evaluate_main_model_iteration(
            m_grm_eigenvals, m_xmat_trans, m_y_trans, varcom);
        spdlog::info("-2logL: {}", iteration_state.logL);

        MatrixXd em_mat = MatrixXd::Zero(2, 2);
        em_mat(0, 0) = m_num_id / (2 * varcom(0) * varcom(0));
        em_mat(1, 1) = m_num_id / (2 * varcom(1) * varcom(1));

        const VarianceUpdateResult update_result =
            find_positive_variance_update(iteration_state.ai_mat, em_mat, iteration_state.fd_mat, varcom);
        spdlog::info("Updated variances: {}", format_vector_for_log(update_result.updated_varcom));

        varcom = update_result.updated_varcom;

        double cc_par_val = 1000.0, cc_gra_val = 1000.0;
        isCC = convergence_criteria(iteration_state.logL, logL, cc_logL0,
                    update_result.delta, cc_par_val, cc_par0,
                    iteration_state.fd_mat, cc_gra_val, cc_gra0);
        logL = iteration_state.logL;
        if (isCC)
            break;
    }

    if (isCC)
        spdlog::info("Variances Converged");
    else
        spdlog::info("Variances Not Converged");
    m_varcom_null = varcom;
    write_variance_components(m_out_file, varcom);

    return m_varcom_null;
}

VectorXd fastGxE::varcom_main_binary(VectorXd& init, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0){
    spdlog::info("Estimate binary-trait variance using logistic mixed model...");
    require_binary_trait_matrix(m_y);

    const VectorXd y = m_y.col(0);
    if (y.minCoeff() < 0.0 || y.maxCoeff() > 1.0) {
        fatal_error("Binary-trait logistic mixed model requires phenotype values in {{0, 1}}.");
    }
    if (y.mean() <= 0.0 || y.mean() >= 1.0) {
        fatal_error("Binary-trait logistic mixed model requires both cases and controls.");
    }

    init = initialize_variance_components(init, 1);
    double genetic_var = std::max(init(0), 1e-6);
    spdlog::info("Initial genetic variance: {}", genetic_var);

    // Step 1 in the working-LMM algorithm: start from a standard logistic
    // regression and set the random effect to zero.
    VectorXd fixed_effect = fit_logistic_fixed_effects(m_xmat, y);
    VectorXd eta = m_xmat * fixed_effect;
    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    bool isCC = false;

    for (int outer_iter = 0; outer_iter < maxiter0; ++outer_iter) {
        spdlog::info("Binary outer iteration {}", outer_iter + 1);

        // Outer loop: update mu, W and y* from the current eta.
        const BinaryWorkingLmmState working_state = build_binary_working_lmm_state(y, eta);

        // Inner loop: with y* and W fixed, iterate the genetic variance until
        // the working-LMM variance component converges.
        const BinaryVarianceFitState variance_fit = fit_binary_working_genetic_variance(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            m_grm_mat,
            m_xmat,
            working_state.working_response,
            working_state.residual_diag,
            genetic_var,
            maxiter0,
            cc_par0,
            cc_gra0,
            cc_logL0,
            m_num_id,
            "Inner");

        const double genetic_var_new = variance_fit.genetic_var;
        const BinaryProjectionState& projection_state = variance_fit.projection_state;
        const VectorXd fixed_effect_new =
            projection_state.fixed_effect_projection.design_cross_inverse *
            (projection_state.fixed_effect_projection.inverse_times_design.transpose() *
             working_state.working_response);
        const VectorXd random_effect_new =
            genetic_var_new * projection_state.genetic_times_projected_response;
        const VectorXd eta_new = m_xmat * fixed_effect_new + random_effect_new;

        const double beta_change =
            (fixed_effect_new - fixed_effect).norm() /
            std::max(1.0, fixed_effect.norm());
        const double eta_change =
            (eta_new - eta).norm() /
            std::sqrt(static_cast<double>(m_num_id));
        const double var_change =
            std::abs(genetic_var_new - genetic_var) / std::max(1.0, std::abs(genetic_var));
        const double logL_change = std::isfinite(previous_neg2_reml)
            ? std::abs(variance_fit.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();

        spdlog::info("Outer {} fitted genetic variance: {}", outer_iter + 1, genetic_var_new);
        if (!variance_fit.converged) {
            spdlog::info("Outer {} inner variance loop hit the maximum number of iterations.", outer_iter + 1);
        }

        fixed_effect = fixed_effect_new;
        eta = eta_new;
        genetic_var = genetic_var_new;
        previous_neg2_reml = variance_fit.neg2_reml;

        isCC = ((var_change < cc_par0) && (beta_change < cc_par0) && (eta_change < cc_par0)) ||
               (logL_change < cc_logL0);
        if (isCC) {
            break;
        }
    }

    if (isCC) {
        spdlog::info("Binary null-model variances converged");
    } else {
        spdlog::info("Binary null-model variances not converged");
    }

    // Finalize with a couple of outer-style refinements so the cached null
    // model is consistent with both the working response and the converged
    // inner variance loop.
    for (int refine_iter = 0; refine_iter < 2; ++refine_iter) {
        const BinaryWorkingLmmState refine_working_state = build_binary_working_lmm_state(y, eta);
        const BinaryVarianceFitState refine_variance_fit = fit_binary_working_genetic_variance(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            m_grm_mat,
            m_xmat,
            refine_working_state.working_response,
            refine_working_state.residual_diag,
            genetic_var,
            maxiter0,
            cc_par0,
            cc_gra0,
            cc_logL0,
            m_num_id,
            "Final inner");
        genetic_var = refine_variance_fit.genetic_var;
        const BinaryProjectionState& refine_projection_state = refine_variance_fit.projection_state;
        fixed_effect =
            refine_projection_state.fixed_effect_projection.design_cross_inverse *
            (refine_projection_state.fixed_effect_projection.inverse_times_design.transpose() *
             refine_working_state.working_response);
        const VectorXd random_effect =
            genetic_var * refine_projection_state.genetic_times_projected_response;
        eta = m_xmat * fixed_effect + random_effect;
    }

    m_binary_eta = eta;
    const BinaryWorkingLmmState final_working_state = build_binary_working_lmm_state(y, m_binary_eta);
    const BinaryVarianceFitState final_variance_fit = fit_binary_working_genetic_variance(
        m_grm_mat_group_vec,
        m_grm_index_vec,
        m_grm_mat,
        m_xmat,
        final_working_state.working_response,
        final_working_state.residual_diag,
        genetic_var,
        maxiter0,
        cc_par0,
        cc_gra0,
        cc_logL0,
        m_num_id,
        "Cached inner");
    genetic_var = final_variance_fit.genetic_var;
    const BinaryProjectionState& final_projection_state = final_variance_fit.projection_state;

    m_binary_mu = final_working_state.mu;
    m_binary_working_response = final_working_state.working_response;
    m_V0_logdet = final_projection_state.inverse_result.logdet;
    m_Vi0 = final_projection_state.inverse_result.inverse_matrix;
    m_varcom_null = VectorXd::Constant(1, genetic_var);
    m_logL_null = m_V0_logdet +
                  final_projection_state.fixed_effect_projection.design_logdet +
                  m_binary_working_response.dot(final_projection_state.projected_response);

    spdlog::info("-2 working REML in binary null model: {}", m_logL_null);
    spdlog::info("Bernoulli log-likelihood at fitted eta: {}",
        bernoulli_loglikelihood(y, m_binary_mu));

    write_variance_components(m_out_file, m_varcom_null);
    return m_varcom_null;
}

VectorXd fastGxE::varcom_main_binary_continuous(
        VectorXd& init, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0) {
    spdlog::info("Estimate joint binary-continuous variances using the bivariate working mixed model...");
    require_binary_continuous_trait_matrix(m_y);

    const MatrixXd original_y = m_y;
    const VectorXd y_binary = m_y.col(0);
    const VectorXd y_continuous = m_y.col(1);
    const MatrixXd joint_design = build_repeated_fixed_effect_design(m_xmat, 2);

    MatrixXd binary_trait = MatrixXd::Zero(m_num_id, 1);
    binary_trait.col(0) = y_binary;
    m_y = binary_trait;
    VectorXd binary_init;
    const VectorXd binary_varcom = varcom_main_binary(
        binary_init, maxiter0, cc_par0, cc_gra0, cc_logL0);

    MatrixXd continuous_trait = MatrixXd::Zero(m_num_id, 1);
    continuous_trait.col(0) = y_continuous;
    m_y = continuous_trait;
    if (m_grm_eigenvals.size() != m_num_id || m_grm_eigenvecs.rows() != m_num_id) {
        process_grm(true);
    } else {
        m_y_trans = m_grm_eigenvecs.transpose() * m_y;
        m_xmat_trans = m_grm_eigenvecs.transpose() * m_xmat;
    }
    VectorXd continuous_init;
    const VectorXd continuous_varcom = varcom_main(
        continuous_init, maxiter0, cc_par0, cc_gra0, cc_logL0);

    m_y = original_y;

    VectorXd beta_binary = fit_logistic_fixed_effects(m_xmat, y_binary);
    VectorXd beta_continuous = fit_linear_fixed_effects(m_xmat, y_continuous);
    VectorXd joint_beta(m_xmat.cols() * 2);
    joint_beta.head(m_xmat.cols()) = beta_binary;
    joint_beta.tail(m_xmat.cols()) = beta_continuous;

    VectorXd eta_binary = m_binary_eta.size() == m_num_id
        ? m_binary_eta
        : VectorXd(m_xmat * beta_binary);

    VectorXd varcom = init;
    if (!is_valid_joint_main_variance_components(varcom)) {
        varcom.resize(4);
        varcom << binary_varcom(0), 0.0, continuous_varcom(0), continuous_varcom(1);
    }
    spdlog::info("Initial binary-continuous variances are: {}", format_vector_for_log(varcom));

    double previous_neg2_reml = std::numeric_limits<double>::infinity();
    bool isCC = false;

    for (int outer_iter = 0; outer_iter < maxiter0; ++outer_iter) {
        spdlog::info("Binary-continuous outer iteration {}", outer_iter + 1);

        const JointMainWorkingState working_state =
            build_joint_main_working_state(m_y, eta_binary);
        const JointVarianceFitState variance_fit = fit_joint_main_working_variance(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            joint_design,
            working_state.stacked_response,
            working_state.binary_residual_diag,
            varcom,
            maxiter0,
            cc_par0,
            cc_gra0,
            cc_logL0,
            m_num_id,
            "Joint inner");

        const VectorXd joint_beta_new =
            variance_fit.projection_state.fixed_effect_projection.design_cross_inverse *
            (variance_fit.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
             working_state.stacked_response);
        const VectorXd joint_random_new = apply_joint_genetic_random_effect(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            variance_fit.projection_state.projected_response,
            variance_fit.varcom,
            m_num_id);
        const VectorXd eta_binary_new =
            m_xmat * joint_beta_new.head(m_xmat.cols()) +
            joint_random_new.head(m_num_id);

        const double beta_change =
            (joint_beta_new - joint_beta).norm() / std::max(1.0, joint_beta.norm());
        const double eta_change =
            (eta_binary_new - eta_binary).norm() / std::sqrt(static_cast<double>(m_num_id));
        const double var_change =
            (variance_fit.varcom - varcom).norm() / std::max(1.0, varcom.norm());
        const double logL_change = std::isfinite(previous_neg2_reml)
            ? std::abs(variance_fit.neg2_reml - previous_neg2_reml)
            : std::numeric_limits<double>::infinity();
        const double gradient_norm = variance_fit.score.norm();
        const bool fixed_effect_stable =
            (var_change < cc_par0) && (beta_change < cc_par0) && (eta_change < cc_par0);
        const bool objective_stable = logL_change < cc_logL0;
        const bool gradient_small = gradient_norm < cc_gra0;

        spdlog::info("Outer {} fitted variances: {}",
            outer_iter + 1, format_vector_for_log(variance_fit.varcom));
        spdlog::info(
            "Outer {} convergence metrics: var_change={}, beta_change={}, eta_change={}, grad_norm={}, logL_change={}",
            outer_iter + 1, var_change, beta_change, eta_change, gradient_norm, logL_change);

        joint_beta = joint_beta_new;
        eta_binary = eta_binary_new;
        varcom = variance_fit.varcom;
        previous_neg2_reml = variance_fit.neg2_reml;

        // Outer-loop convergence also requires the inner score to be small.
        isCC = variance_fit.converged && gradient_small &&
               (fixed_effect_stable || objective_stable);
        if (isCC) {
            break;
        }
    }

    if (isCC) {
        spdlog::info("Binary-continuous null-model variances converged");
    } else {
        spdlog::info("Binary-continuous null-model variances not converged");
    }

    for (int refine_iter = 0; refine_iter < 2; ++refine_iter) {
        const JointMainWorkingState refine_working_state =
            build_joint_main_working_state(m_y, eta_binary);
        const JointVarianceFitState refine_variance_fit = fit_joint_main_working_variance(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            joint_design,
            refine_working_state.stacked_response,
            refine_working_state.binary_residual_diag,
            varcom,
            maxiter0,
            cc_par0,
            cc_gra0,
            cc_logL0,
            m_num_id,
            "Joint final inner");
        joint_beta =
            refine_variance_fit.projection_state.fixed_effect_projection.design_cross_inverse *
            (refine_variance_fit.projection_state.fixed_effect_projection.inverse_times_design.transpose() *
             refine_working_state.stacked_response);
        const VectorXd joint_random = apply_joint_genetic_random_effect(
            m_grm_mat_group_vec,
            m_grm_index_vec,
            refine_variance_fit.projection_state.projected_response,
            refine_variance_fit.varcom,
            m_num_id);
        eta_binary = m_xmat * joint_beta.head(m_xmat.cols()) +
                     joint_random.head(m_num_id);
        varcom = refine_variance_fit.varcom;
    }

    const JointMainWorkingState final_working_state =
        build_joint_main_working_state(m_y, eta_binary);
    const JointVarianceFitState final_variance_fit = fit_joint_main_working_variance(
        m_grm_mat_group_vec,
        m_grm_index_vec,
        joint_design,
        final_working_state.stacked_response,
        final_working_state.binary_residual_diag,
        varcom,
        maxiter0,
        cc_par0,
        cc_gra0,
        cc_logL0,
        m_num_id,
        "Joint cached inner");

    m_joint_binary_eta = eta_binary;
    m_joint_binary_mu = final_working_state.binary_mu;
    m_joint_working_response = final_working_state.stacked_response;
    m_V0_logdet = final_variance_fit.projection_state.inverse_result.logdet;
    m_Vi0 = final_variance_fit.projection_state.inverse_result.inverse_matrix;
    m_varcom_null = final_variance_fit.varcom;
    m_logL_null = final_variance_fit.neg2_reml;

    spdlog::info("-2 working REML in binary-continuous null model: {}", m_logL_null);
    spdlog::info("Bernoulli log-likelihood at fitted eta: {}",
        bernoulli_loglikelihood(y_binary, m_joint_binary_mu));

    write_variance_components(m_out_file, m_varcom_null);
    return m_varcom_null;
}

VectorXd fastGxE::varcom_main_multitrait_continuous(
        VectorXd& init, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0) {
    spdlog::info("Estimate multitrait continuous variances using the multivariate mixed model...");
    require_multitrait_continuous_trait_matrix(m_y);

    const MatrixXd original_y = m_y;
    const std::int64_t num_traits = original_y.cols();
    const std::int64_t num_covariance_terms = num_lower_triangle_elements(num_traits);
    const MatrixXd joint_design = build_repeated_fixed_effect_design(m_xmat, num_traits);

    MatrixXd genetic_covariance = MatrixXd::Zero(num_traits, num_traits);
    MatrixXd residual_covariance = MatrixXd::Zero(num_traits, num_traits);
    for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
        MatrixXd single_trait_matrix = MatrixXd::Zero(m_num_id, 1);
        single_trait_matrix.col(0) = original_y.col(trait_index);
        m_y = single_trait_matrix;
        if (m_grm_eigenvals.size() != m_num_id || m_grm_eigenvecs.rows() != m_num_id) {
            process_grm(true);
        } else {
            m_y_trans = m_grm_eigenvecs.transpose() * m_y;
            m_xmat_trans = m_grm_eigenvecs.transpose() * m_xmat;
        }

        VectorXd single_trait_init;
        const VectorXd single_trait_varcom = varcom_main(
            single_trait_init, maxiter0, cc_par0, cc_gra0, cc_logL0);
        genetic_covariance(trait_index, trait_index) = single_trait_varcom(0);
        residual_covariance(trait_index, trait_index) = single_trait_varcom(1);
    }

    m_y = original_y;
    const VectorXd stacked_response = stack_trait_matrix_rows(original_y);

    VectorXd varcom = init;
    if (!is_valid_multitrait_continuous_variance_components(varcom, num_traits)) {
        varcom.resize(2 * num_covariance_terms);
        varcom.head(num_covariance_terms) = pack_symmetric_matrix_lower_triangle(genetic_covariance);
        varcom.tail(num_covariance_terms) = pack_symmetric_matrix_lower_triangle(residual_covariance);
    }
    spdlog::info("Initial multitrait continuous variances are: {}", format_vector_for_log(varcom));

    const MultitraitContinuousVarianceFitState variance_fit = fit_multitrait_continuous_variance(
        m_grm_mat_group_vec,
        m_grm_index_vec,
        joint_design,
        stacked_response,
        varcom,
        maxiter0,
        cc_par0,
        cc_gra0,
        cc_logL0,
        m_num_id,
        num_traits,
        "Multitrait continuous");

    if (variance_fit.converged) {
        spdlog::info("Multitrait continuous null-model variances converged");
    } else {
        spdlog::info("Multitrait continuous null-model variances not converged");
    }

    m_joint_working_response = stacked_response;
    m_V0_logdet = variance_fit.projection_state.inverse_result.logdet;
    m_Vi0 = variance_fit.projection_state.inverse_result.inverse_matrix;
    m_varcom_null = variance_fit.varcom;
    m_logL_null = variance_fit.neg2_reml;

    spdlog::info("-2 REML in multitrait continuous null model: {}", m_logL_null);
    write_variance_components(m_out_file, m_varcom_null);
    return m_varcom_null;
}



/**
 * @Description: Test SNPs with one random-effect model
 */
void fastGxE::test_main(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut, 
                int maxiter, double cc_par, double cc_gra){
    (void)speed;
    (void)p_approx_cut;
    (void)maxiter;
    (void)cc_par;
    (void)cc_gra;

    spdlog::info("Assocation...");
    GENO GenoA(bed_file);
    const BedScanContext bed_context = prepare_bed_scan_context(GenoA, m_id_in_data_vec);

    // Project the fixed-effect design into the GRM eigen basis once, then keep
    // reusing the resulting null-model quantities across all SNP partitions.
    m_xmat_trans = m_grm_eigenvecs.transpose() * m_xmat;

    // V = Q diag(λ_i σ²_g + σ²_e) Q' in GRM eigenbasis; diagonal form avoids an O(n²) dense inverse.
    const VectorXd inverse_diagonal =
        1.0 / (m_grm_eigenvals.array() * m_varcom_null(0) + m_varcom_null(1));
    const ProjectionCache null_projection =
        build_projection_cache(inverse_diagonal, m_xmat_trans);
    
    spdlog::info("Randomly select SNPs and calculate the gamma");
    const vector<std::int64_t> random_snp_index_vec = sample_random_snp_indices(bed_context.num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
    load_centered_random_snp_panel(
        GenoA, bed_file, random_snp_index_vec, bed_context.id_index_in_bed_vec,
        snp_mat_part_random, freq_arr, missing_rate_arr, nobs_geno_arr);
    snp_mat_part_random = m_grm_eigenvecs.transpose() * snp_mat_part_random;

    VectorXd gamma_correct_Vec(num_random_snp);
    VectorXd p_Vec(num_random_snp);
    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t k = 0; k < num_random_snp; k++){
        VectorXd isnp_arr = snp_mat_part_random.col(k);
        // !score test
        VectorXd psnp = apply_projection(
            inverse_diagonal,
            null_projection.inverse_times_design,
            null_projection.design_cross_inverse,
            isnp_arr);
        double p_score = gsl_cdf_chisq_Q(std::pow((psnp.transpose() * m_y_trans).sum(), 2) / (isnp_arr.transpose() * psnp).sum(), 1);

        double geno_var_isnp = isnp_arr.transpose() * isnp_arr;
        gamma_correct_Vec(k) = (isnp_arr.transpose() * psnp).sum() / geno_var_isnp;
        p_Vec(k) = p_score;
    }
    snp_mat_part_random.resize(0, 0);

    std::int64_t num_random_snp_used = 0;
    double gamma_correct = 0;
    ofstream fout = open_output_stream(m_out_file + ".random.res");
    fout << "order af missing gamma p" << std::endl;
    for(std::int64_t k = 0; k < num_random_snp; k++){
        if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
            fout << random_snp_index_vec[k] << " " << freq_arr(k) << " " << 
                    missing_rate_arr(k) << " " << gamma_correct_Vec(k) << " " << p_Vec(k) << std::endl;
            num_random_snp_used++;
            gamma_correct += gamma_correct_Vec[k];
        }
    }
    fout.close();

    require_random_snp_fraction(
        num_random_snp_used, num_random_snp, 0.1,
        "The number of used random SNPs to calculate gamma is: {}, which is less than 10% of random SNPs");
    gamma_correct /= num_random_snp_used;
    spdlog::info("Gamma factor: {}", gamma_correct);
    spdlog::info("The number of used random SNPs is: {}", num_random_snp_used);

    // Move back to the sample space for fast genome-wide SNP scanning once the
    // random-SNP calibration factor has been estimated.
    spdlog::info("Start testing SNPs...");
    fout = open_output_stream(m_out_file + ".res");
    fout << "order chrom SNP cm base allele1 allele2 af missing beta se p" << std::endl;
    
    const VectorXd pyarr0 = m_grm_eigenvecs * apply_projection(
        inverse_diagonal,
        null_projection.inverse_times_design,
        null_projection.design_cross_inverse,
        m_y_trans).col(0);
    const double* pyarr0_pt = pyarr0.data();

    SnpScalarScanBuffers bufs;
    bufs.resize(get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read, bed_context.num_id_used);

    for(std::int64_t i = 0; i < npart_snp; i++){
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, i);
        if(bufs.needs_resize(partition.num_snp_read)){
            bufs.resize(partition.num_snp_read, bed_context.num_id_used);
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}",
             i + 1, npart_snp, partition.start_snp, partition.num_snp_read);

        // Each partition reuses the same work buffers; only the logical SNP span
        // changes as we stream the genotype matrix block by block.
        GenoA.read_bed_centered_to_buffer_omp(bed_file + ".bed", bufs.snp_mat.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
                partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);

        // Test
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            bufs.geno_var[j] = cblas_ddot(
                bed_context.num_id_used, bufs.snp_mat.data() + j * bed_context.num_id_used, 1,
                bufs.snp_mat.data() + j * bed_context.num_id_used, 1);
        }

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            bufs.se[j] = 1.0 / (gamma_correct * bufs.geno_var[j]);
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, partition.num_snp_read, bed_context.num_id_used, 1.0, bufs.snp_mat.data(),
                bed_context.num_id_used, pyarr0_pt, 1, 0.0, bufs.prod.data(), 1);

        vdMul(partition.num_snp_read, bufs.se.data(), bufs.prod.data(), bufs.eff.data());

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            bufs.se[j] = sqrt(bufs.se[j]);
        }

        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t k = 0; k < partition.num_snp_read; k++) {
            if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                double eff = bufs.eff[k];
                double se = bufs.se[k];
                bufs.p(k) = gsl_cdf_chisq_Q(eff*eff/(se*se), 1);
            }else{
                bufs.eff[k] = NAN;
                bufs.se[k] = NAN;
                bufs.p(k) = NAN;
            }
        }
        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        fastGxE::output(fout, snp_info_vec, freq_arr, missing_rate_arr,
            bufs.eff.data(), bufs.se.data(), bufs.p, partition.start_snp, partition.num_snp_read, p_cut);
    }
    fout.close();
}

void fastGxE::test_main_binary(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut,
                int maxiter, double cc_par, double cc_gra){
    (void)speed;
    (void)p_approx_cut;
    (void)maxiter;
    (void)cc_par;
    (void)cc_gra;

    spdlog::info("Association for binary trait using logistic mixed model...");
    require_binary_trait_matrix(m_y);
    if (m_binary_working_response.size() != m_num_id) {
        fatal_error("Binary null-model working response is unavailable. Run varcom_main_binary first.");
    }

    GENO GenoA(bed_file);
    const BedScanContext bed_context = prepare_bed_scan_context(GenoA, m_id_in_data_vec);

    const ProjectionCache null_projection = build_projection_cache(m_Vi0, m_xmat);
    const VectorXd projected_working_response = apply_projection(
        m_Vi0,
        null_projection.inverse_times_design,
        null_projection.design_cross_inverse,
        m_binary_working_response);

    spdlog::info("Randomly select SNPs and calibrate the binary score variance");
    const vector<std::int64_t> random_snp_index_vec = sample_random_snp_indices(bed_context.num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
    load_centered_random_snp_panel(
        GenoA, bed_file, random_snp_index_vec, bed_context.id_index_in_bed_vec,
        snp_mat_part_random, freq_arr, missing_rate_arr, nobs_geno_arr);

    VectorXd gamma_correct_Vec = VectorXd::Constant(num_random_snp, std::numeric_limits<double>::quiet_NaN());
    VectorXd p_Vec = VectorXd::Constant(num_random_snp, std::numeric_limits<double>::quiet_NaN());
    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t k = 0; k < num_random_snp; ++k){
        const VectorXd snp_vec = snp_mat_part_random.col(k);
        const VectorXd projected_snp = apply_projection(
            m_Vi0,
            null_projection.inverse_times_design,
            null_projection.design_cross_inverse,
            snp_vec);
        const double geno_var = snp_vec.squaredNorm();
        const double score = snp_vec.dot(projected_working_response);
        const double projected_var = snp_vec.dot(projected_snp);
        gamma_correct_Vec(k) = projected_var / geno_var;
        if (projected_var > 0.0) {
            p_Vec(k) = gsl_cdf_chisq_Q(score * score / projected_var, 1);
        }
    }
    snp_mat_part_random.resize(0, 0);

    std::int64_t num_random_snp_used = 0;
    double gamma_correct = 0.0;
    ofstream fout = open_output_stream(m_out_file + ".random.res");
    fout << "order af missing gamma p" << std::endl;
    for(std::int64_t k = 0; k < num_random_snp; ++k){
        if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut) &&
                std::isfinite(gamma_correct_Vec(k)) && gamma_correct_Vec(k) > 0.0){
            fout << random_snp_index_vec[k] << " " << freq_arr(k) << " "
                 << missing_rate_arr(k) << " " << gamma_correct_Vec(k) << " "
                 << p_Vec(k) << std::endl;
            ++num_random_snp_used;
            gamma_correct += gamma_correct_Vec(k);
        }
    }
    fout.close();

    require_random_snp_fraction(
        num_random_snp_used, num_random_snp, 0.1,
        "The number of used random SNPs to calibrate the binary score is: {}, which is less than 10% of random SNPs");
    gamma_correct /= num_random_snp_used;
    if (!(gamma_correct > 0.0)) {
        fatal_error("The calibrated binary score variance factor must be positive.");
    }

    spdlog::info("Binary gamma factor: {}", gamma_correct);
    spdlog::info("The number of used random SNPs is: {}", num_random_snp_used);

    spdlog::info("Start testing SNPs for the binary trait...");
    fout = open_output_stream(m_out_file + ".res");
    fout << "order chrom SNP cm base allele1 allele2 af missing beta se p" << std::endl;

    SnpScalarScanBuffers bufs;
    bufs.resize(get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read, bed_context.num_id_used);

    for(std::int64_t i = 0; i < npart_snp; ++i){
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, i);
        if(bufs.needs_resize(partition.num_snp_read)){
            bufs.resize(partition.num_snp_read, bed_context.num_id_used);
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}",
             i + 1, npart_snp, partition.start_snp, partition.num_snp_read);

        GenoA.read_bed_centered_to_buffer_omp(
            bed_file + ".bed", bufs.snp_mat.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
            partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            bufs.geno_var[j] = cblas_ddot(
                bed_context.num_id_used, bufs.snp_mat.data() + j * bed_context.num_id_used, 1,
                bufs.snp_mat.data() + j * bed_context.num_id_used, 1);
            bufs.se[j] = 1.0 / (gamma_correct * bufs.geno_var[j]);
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, partition.num_snp_read, bed_context.num_id_used, 1.0,
            bufs.snp_mat.data(), bed_context.num_id_used, projected_working_response.data(), 1, 0.0, bufs.prod.data(), 1);
        vdMul(partition.num_snp_read, bufs.se.data(), bufs.prod.data(), bufs.eff.data());

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            const double variance = bufs.se[j];
            if (std::isfinite(variance) && variance > 0.0) {
                bufs.se[j] = std::sqrt(variance);
            } else {
                bufs.se[j] = NAN;
            }
        }

        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t k = 0; k < partition.num_snp_read; ++k) {
            if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut) &&
                    std::isfinite(bufs.eff[k]) && std::isfinite(bufs.se[k]) && bufs.se[k] > 0.0) {
                const double beta = bufs.eff[k];
                const double se = bufs.se[k];
                bufs.p(k) = gsl_cdf_chisq_Q(beta * beta / (se * se), 1);
            } else {
                bufs.eff[k] = NAN;
                bufs.se[k] = NAN;
                bufs.p(k) = NAN;
            }
        }

        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        fastGxE::output(
            fout, snp_info_vec, freq_arr, missing_rate_arr, bufs.eff.data(), bufs.se.data(),
            bufs.p, partition.start_snp, partition.num_snp_read, p_cut);
    }
    fout.close();
}

void fastGxE::test_main_binary_continuous(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut,
                int maxiter, double cc_par, double cc_gra) {
    (void)speed;
    (void)num_random_snp;
    (void)p_approx_cut;
    (void)p_cut;
    (void)maxiter;
    (void)cc_par;
    (void)cc_gra;

    spdlog::info("Association for the joint binary-continuous model...");
    require_binary_continuous_trait_matrix(m_y);
    if (m_joint_working_response.size() != 2 * m_num_id) {
        fatal_error("Joint binary-continuous working response is unavailable. Run a joint binary-continuous null model first.");
    }

    GENO GenoA(bed_file);
    const BedScanContext bed_context = prepare_bed_scan_context(GenoA, m_id_in_data_vec);
    const MatrixXd joint_design = build_repeated_fixed_effect_design(m_xmat, 2);
    const ProjectionCache null_projection = build_projection_cache(m_Vi0, joint_design);
    const VectorXd projected_joint_response = apply_projection(
        m_Vi0,
        null_projection.inverse_times_design,
        null_projection.design_cross_inverse,
        m_joint_working_response);

    ofstream fout = open_output_stream(m_out_file + ".res");
    fout << "order chrom SNP cm base allele1 allele2 af missing"
         << " beta_binary se_binary p_binary"
         << " beta_continuous se_continuous p_continuous"
         << " p_joint" << std::endl;

    SnpMatrixScanBuffers bufs;
    bufs.resize(get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read,
                bed_context.num_id_used, 2, 3);
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;

    for (std::int64_t ipart = 0; ipart < npart_snp; ++ipart) {
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, ipart);
        if (bufs.needs_resize(partition.num_snp_read)) {
            bufs.resize(partition.num_snp_read, bed_context.num_id_used, 2, 3);
        } else {
            bufs.reset_values();
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}",
            ipart + 1, npart_snp, partition.start_snp, partition.num_snp_read);

        GenoA.read_bed_centered_to_buffer_omp(
            bed_file + ".bed", bufs.snp_mat.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
            partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t k = 0; k < partition.num_snp_read; ++k) {
            if (!isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)) {
                continue;
            }

            Eigen::Map<Eigen::VectorXd> snp_vec(
                bufs.snp_mat.data() + k * bed_context.num_id_used, bed_context.num_id_used);
            MatrixXd snp_design = MatrixXd::Zero(2 * m_num_id, 2);
            snp_design.block(0, 0, m_num_id, 1) = snp_vec;
            snp_design.block(m_num_id, 1, m_num_id, 1) = snp_vec;

            const MatrixXd projected_snp_design = apply_projection(
                m_Vi0,
                null_projection.inverse_times_design,
                null_projection.design_cross_inverse,
                snp_design);
            const MatrixXd information = snp_design.transpose() * projected_snp_design;
            Eigen::LDLT<MatrixXd> ldlt(information);
            if (ldlt.info() != Eigen::Success) {
                continue;
            }

            const VectorXd score = snp_design.transpose() * projected_joint_response;
            const VectorXd beta = ldlt.solve(score);
            const MatrixXd covariance = ldlt.solve(MatrixXd::Identity(2, 2));
            if (!beta.allFinite() || !covariance.allFinite() ||
                    covariance(0, 0) <= 0.0 || covariance(1, 1) <= 0.0) {
                continue;
            }

            bufs.beta_mat.row(k) = beta.transpose();
            bufs.se_mat(k, 0) = std::sqrt(covariance(0, 0));
            bufs.se_mat(k, 1) = std::sqrt(covariance(1, 1));
            bufs.p_mat(k, 0) = gsl_cdf_chisq_Q(beta(0) * beta(0) / covariance(0, 0), 1);
            bufs.p_mat(k, 1) = gsl_cdf_chisq_Q(beta(1) * beta(1) / covariance(1, 1), 1);
            const double joint_chisq = std::max(score.dot(beta), 0.0);
            bufs.p_mat(k, 2) = gsl_cdf_chisq_Q(joint_chisq, 2);
        }

        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        for (std::int64_t k = 0; k < partition.num_snp_read; ++k) {
            fout << partition.start_snp + k << " " << snp_info_vec[k] << " "
                 << freq_arr(k) << " " << missing_rate_arr(k) << " "
                 << bufs.beta_mat(k, 0) << " " << bufs.se_mat(k, 0) << " " << bufs.p_mat(k, 0) << " "
                 << bufs.beta_mat(k, 1) << " " << bufs.se_mat(k, 1) << " " << bufs.p_mat(k, 1) << " "
                 << bufs.p_mat(k, 2) << std::endl;
        }
    }

    fout.close();
}

void fastGxE::test_main_multitrait_continuous(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut,
                int maxiter, double cc_par, double cc_gra) {
    (void)speed;
    (void)num_random_snp;
    (void)p_approx_cut;
    (void)p_cut;
    (void)maxiter;
    (void)cc_par;
    (void)cc_gra;

    spdlog::info("Association for the multitrait continuous model...");
    require_multitrait_continuous_trait_matrix(m_y);
    const std::int64_t num_traits = m_y.cols();
    if (m_joint_working_response.size() != num_traits * m_num_id) {
        fatal_error("Multitrait continuous response is unavailable. Run varcom_main_multitrait_continuous first.");
    }

    GENO GenoA(bed_file);
    const BedScanContext bed_context = prepare_bed_scan_context(GenoA, m_id_in_data_vec);
    const MatrixXd joint_design = build_repeated_fixed_effect_design(m_xmat, num_traits);
    const ProjectionCache null_projection = build_projection_cache(m_Vi0, joint_design);
    const VectorXd projected_joint_response = apply_projection(
        m_Vi0,
        null_projection.inverse_times_design,
        null_projection.design_cross_inverse,
        m_joint_working_response);

    ofstream fout = open_output_stream(m_out_file + ".res");
    fout << "order chrom SNP cm base allele1 allele2 af missing";
    for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
        const std::string raw_label =
            trait_index < static_cast<std::int64_t>(m_trait_name_vec.size())
                ? m_trait_name_vec[trait_index]
                : ("trait" + std::to_string(trait_index + 1));
        const std::string trait_label = sanitize_output_label(raw_label);
        fout << " beta_" << trait_label
             << " se_" << trait_label
             << " p_" << trait_label;
    }
    fout << " p_joint" << std::endl;

    SnpMatrixScanBuffers bufs;
    bufs.resize(get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read,
                bed_context.num_id_used, num_traits, num_traits + 1);
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;

    for (std::int64_t ipart = 0; ipart < npart_snp; ++ipart) {
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, ipart);
        if (bufs.needs_resize(partition.num_snp_read)) {
            bufs.resize(partition.num_snp_read, bed_context.num_id_used,
                        num_traits, num_traits + 1);
        } else {
            bufs.reset_values();
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}",
            ipart + 1, npart_snp, partition.start_snp, partition.num_snp_read);

        GenoA.read_bed_centered_to_buffer_omp(
            bed_file + ".bed", bufs.snp_mat.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
            partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t k = 0; k < partition.num_snp_read; ++k) {
            if (!isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)) {
                continue;
            }

            Eigen::Map<Eigen::VectorXd> snp_vec(
                bufs.snp_mat.data() + k * bed_context.num_id_used, bed_context.num_id_used);
            MatrixXd snp_design = MatrixXd::Zero(num_traits * m_num_id, num_traits);
            for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
                snp_design.block(trait_index * m_num_id, trait_index, m_num_id, 1) = snp_vec;
            }

            const MatrixXd projected_snp_design = apply_projection(
                m_Vi0,
                null_projection.inverse_times_design,
                null_projection.design_cross_inverse,
                snp_design);
            const MatrixXd information = snp_design.transpose() * projected_snp_design;
            Eigen::LDLT<MatrixXd> ldlt(information);
            if (ldlt.info() != Eigen::Success) {
                continue;
            }

            const VectorXd score = snp_design.transpose() * projected_joint_response;
            const VectorXd beta = ldlt.solve(score);
            const MatrixXd covariance = ldlt.solve(MatrixXd::Identity(num_traits, num_traits));
            if (!beta.allFinite() || !covariance.allFinite()) {
                continue;
            }

            bufs.beta_mat.row(k) = beta.transpose();
            bool has_valid_standard_errors = true;
            for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
                if (covariance(trait_index, trait_index) <= 0.0) {
                    has_valid_standard_errors = false;
                    break;
                }
                bufs.se_mat(k, trait_index) = std::sqrt(covariance(trait_index, trait_index));
                bufs.p_mat(k, trait_index) = gsl_cdf_chisq_Q(
                    beta(trait_index) * beta(trait_index) / covariance(trait_index, trait_index), 1);
            }
            if (!has_valid_standard_errors) {
                bufs.beta_mat.row(k).setConstant(std::numeric_limits<double>::quiet_NaN());
                bufs.se_mat.row(k).setConstant(std::numeric_limits<double>::quiet_NaN());
                bufs.p_mat.row(k).setConstant(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            const double joint_chisq = std::max(score.dot(beta), 0.0);
            bufs.p_mat(k, num_traits) = gsl_cdf_chisq_Q(joint_chisq, num_traits);
        }

        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        for (std::int64_t k = 0; k < partition.num_snp_read; ++k) {
            fout << partition.start_snp + k << " " << snp_info_vec[k] << " "
                 << freq_arr(k) << " " << missing_rate_arr(k);
            for (std::int64_t trait_index = 0; trait_index < num_traits; ++trait_index) {
                fout << " " << bufs.beta_mat(k, trait_index)
                     << " " << bufs.se_mat(k, trait_index)
                     << " " << bufs.p_mat(k, trait_index);
            }
            fout << " " << bufs.p_mat(k, num_traits) << std::endl;
        }
    }

    fout.close();
}


void fastGxE::output(ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
                const double* eff_arr, const double* se_arr, const VectorXd& p_arr,
                std::int64_t start_snp, std::int64_t num_snp_read, double p_cut){
    for(std::int64_t i = 0; i < num_snp_read; i++){
        fout << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " "  << missing_rate_arr(i) << " "
                << " " << eff_arr[i] << " " << se_arr[i] 
                << " "<< p_arr(i) << std::endl;
    }
}

void fastGxE::reset_output_prefix(const string& out_file, std::int64_t first, std::int64_t second){
    m_out_file = out_file + "." + std::to_string(first) + "_"  + std::to_string(second);
}

bool fastGxE::isValidSnp(double af, double missing_rate, double maf_cut, double missing_rate_cut) {
    return af > maf_cut && af < 1 - maf_cut && missing_rate < missing_rate_cut;
}



/**
 * @Description: Assocation of GxE using an improved structLMM model
 */
void fastGxE::test_GxE_multi(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut){
    spdlog::info("Assocation...");
    GENO GenoA(bed_file);
    vector<std::int64_t> id_index_in_bed_vec = GenoA.find_fam_index(m_id_in_data_vec);
    std::int64_t num_id_used = id_index_in_bed_vec.size();

    std::int64_t num_snp = GenoA.get_num_sid(); // the number of SNP in the plink file.
    std::int64_t num_id_in_bed = GenoA.get_num_iid(); // the number of IDs in the plink file.
    spdlog::info("The number of individuals and SNPs in the plink file are: {} {}", num_id_in_bed, num_snp);
    GenoA.validate_bed_size();

    VectorXd Vi0y = this->m_Vi0 * m_y;
    double *bye_mat_pt = this->m_bye_mat.data();
    double *Vi0y_pt = Vi0y.data();

    spdlog::info("Randomly select SNPs and calculate the gamma");
    vector <std::int64_t> random_snp_index_vec = sample_unique_integers(0, num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
    GenoA.read_bed_by_snp_indices(bed_file + ".bed", snp_mat_part_random, freq_arr, missing_rate_arr, nobs_geno_arr, random_snp_index_vec, id_index_in_bed_vec);
    snp_mat_part_random.rowwise() -= 2*freq_arr.transpose();

    // main 
    VectorXd beta_main_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd se_main_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd p_main_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd gamma_main_Vec= VectorXd::Zero(num_random_snp);

    // GxE without main
    VectorXd score_GxEnoMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd p_GxEnoMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
    vector<MatrixXd> gamma_GxEnoMain_vec(num_random_snp);

    // GxE with main
    VectorXd score_GxEwithMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd p_GxEwithMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
    vector<MatrixXd> gamma_GxEwithMain_vec(num_random_snp);

    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t i = 0; i < num_random_snp; i++){
        if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
            VectorXd snpk = snp_mat_part_random.col(i);
            double snp_norm2 = snpk.transpose() * snpk;

            // main
            double beta_main_var = 1 / (snpk.transpose() * this->m_Vi0 * snpk);
            se_main_Vec(i) = sqrt(beta_main_var);
            double beta_main = beta_main_var * snpk.transpose() * Vi0y;
            beta_main_Vec(i) = beta_main;
            double p_main = gsl_cdf_chisq_Q(beta_main * beta_main / beta_main_var, 1);
            p_main_Vec(i) = p_main;
            gamma_main_Vec(i) = (1 / beta_main_var) / snp_norm2;

            // GxE without main effect
            MatrixXd EoG = this->m_bye_mat.array().colwise() * snpk.array();
            VectorXd EoGtViy = EoG.transpose() * Vi0y;
            double score = EoGtViy.transpose() * EoGtViy;
            score_GxEnoMain_Vec(i) = score;

            MatrixXd EoGtPEoG = EoG.transpose() * this->m_Vi0 * EoG;
            gamma_GxEnoMain_vec[i] = EoGtPEoG / snp_norm2;
            Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
            eigensolver.compute(EoGtPEoG);
            VectorXd a = eigensolver.eigenvalues().real();
            p_GxEnoMain_Vec(i) = saddle(score, a);

            // GxE with main effect
            ScalarProjectionCache snp_main_projection = build_scalar_projection_cache(this->m_Vi0, snpk);
            VectorXd projected_null_response = apply_projection(
                this->m_Vi0,
                snp_main_projection.inverse_times_design,
                snp_main_projection.design_cross_inverse,
                m_y).col(0);

            VectorXd EoGtPy = EoG.transpose() * projected_null_response;
            score = EoGtPy.transpose() * EoGtPy;
            score_GxEwithMain_Vec(i) = score;
            EoGtPEoG = EoG.transpose() * apply_projection(
                this->m_Vi0,
                snp_main_projection.inverse_times_design,
                snp_main_projection.design_cross_inverse,
                EoG);
            gamma_GxEwithMain_vec[i] = EoGtPEoG / snp_norm2;
            eigensolver.compute(EoGtPEoG);
            a = eigensolver.eigenvalues().real();
            p_GxEwithMain_Vec(i) = saddle(score, a);
        }
    }
    snp_mat_part_random.resize(0, 0);

    std::int64_t num_random_snp_used = 0;
    double gamma_main = 0;
    MatrixXd gamma_GxEnoMain = MatrixXd::Zero(m_num_bye, m_num_bye);
    MatrixXd gamma_GxEwithMain = MatrixXd::Zero(m_num_bye, m_num_bye);
    ofstream fout;
    fout.open(m_out_file + ".random.res");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.random.res", m_out_file);
    }

    fout << "order af missing beta se p_main score_noMain p_noMain score_withMain p_withMain" << std::endl;
    for(std::int64_t i = 0; i < num_random_snp; i++){
        if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
            fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) 
            << " " << beta_main_Vec(i) << " " << se_main_Vec(i) << " " << p_main_Vec(i) 
            << " " << score_GxEnoMain_Vec(i) << " " << p_GxEnoMain_Vec(i)
            << " " << score_GxEwithMain_Vec(i) << " " << p_GxEwithMain_Vec(i) << std::endl;
            gamma_main += gamma_main_Vec(i);
            gamma_GxEnoMain += gamma_GxEnoMain_vec[i];
            gamma_GxEwithMain += gamma_GxEwithMain_vec[i];
            num_random_snp_used++;
        }
    }
    fout.close();
    
    if(num_random_snp_used * 1.0 / num_random_snp < 0.5){
        fatal_error("The number of used random SNPs to calculate gamma is: {}, which is less than 50% of random SNPs",
                num_random_snp_used);
    }
    gamma_main /= num_random_snp_used;
    gamma_GxEnoMain /= num_random_snp_used;
    gamma_GxEwithMain /= num_random_snp_used;
    spdlog::info("The number of used random SNPs is: {}", num_random_snp_used);
    
    Eigen::VectorXd chi2_coef_GxEnoMain_Vec; // coefficient weights for chi2
    Eigen::VectorXd chi2_coef_GxEwithMain_Vec;
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
    eigensolver.compute(gamma_GxEnoMain);
    chi2_coef_GxEnoMain_Vec = eigensolver.eigenvalues().real();
    eigensolver.compute(gamma_GxEwithMain);
    chi2_coef_GxEwithMain_Vec = eigensolver.eigenvalues().real();
    
    fout.open(m_out_file + ".res");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.res", m_out_file);
    }
    fout << "order chrom SNP cm base allele1 allele2 af missing beta se p_main score p_gxe" << std::endl;

    const int _omp_max_threads = omp_get_max_threads();
    std::vector<std::vector<double>> EP0y_part(_omp_max_threads, std::vector<double>(num_id_used));
    std::vector<std::vector<double>> WEP0y_part(_omp_max_threads, std::vector<double>(m_num_bye));

    std::int64_t num_snp_read = (std::int64_t)((end_pos - start_pos)/npart_snp);
    std::vector<double> snp_mat_part(num_snp_read * num_id_used);
    beta_main_Vec.resize(num_snp_read);
    se_main_Vec.resize(num_snp_read);
    p_main_Vec.resize(num_snp_read);
    VectorXd score_Vec(num_snp_read);
    VectorXd p_gxe_Vec(num_snp_read);
    
    for(std::int64_t ipart = 0; ipart < npart_snp; ipart++){
        std::int64_t num_snp_read = (std::int64_t)((end_pos - start_pos)/npart_snp);
        std::int64_t start_snp = ipart*num_snp_read + start_pos;

        if(ipart == npart_snp - 1 && (end_pos - start_pos) % npart_snp != 0){
            num_snp_read = num_snp_read + (end_pos - start_pos) % npart_snp;
            snp_mat_part.resize(num_snp_read * num_id_used);
            beta_main_Vec.resize(num_snp_read);
            se_main_Vec.resize(num_snp_read);
            p_main_Vec.resize(num_snp_read);
            score_Vec.resize(num_snp_read);
            p_gxe_Vec.resize(num_snp_read);
        }
        
        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}", 
             ipart + 1, npart_snp, start_snp, num_snp_read);
        
        GenoA.read_bed_centered_to_buffer_omp(bed_file + ".bed", snp_mat_part.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
                start_snp, num_snp_read, id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for(std::int64_t k = 0; k < num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                std::int64_t thread_id = omp_get_thread_num();

                // snpk.T * snpk
                double snpk_norm2 = cblas_ddot(num_id_used, snp_mat_part.data() + k*num_id_used, 1,
                                    snp_mat_part.data() + k*num_id_used, 1);

                // main
                double beta_main_var = 1 / (gamma_main * snpk_norm2);
                se_main_Vec(k) = sqrt(beta_main_var);
                double beta_main = cblas_ddot(num_id_used, snp_mat_part.data() + k*num_id_used, 1,
                                    Vi0y_pt, 1) * beta_main_var;
                beta_main_Vec(k) = beta_main;
                double p_main = gsl_cdf_chisq_Q(beta_main * beta_main / beta_main_var, 1);
                p_main_Vec(k) = p_main;

                // GxE without main effect
                vdMul(num_id_used, Vi0y_pt, snp_mat_part.data() + k*num_id_used, EP0y_part[thread_id].data());
                cblas_dgemv(CblasColMajor, CblasTrans, num_id_used, m_num_bye, 1.0, bye_mat_pt, num_id_used,
                            EP0y_part[thread_id].data(), 1, 0.0, WEP0y_part[thread_id].data(), 1);
                double score = cblas_ddot(m_num_bye, WEP0y_part[thread_id].data(), 1, WEP0y_part[thread_id].data(), 1);
                VectorXd chi_coef_Vec_kth = chi2_coef_GxEnoMain_Vec * snpk_norm2;
                double p_gxe = saddle(score, chi_coef_Vec_kth);

                if(p_main < p_approx_cut || p_gxe < p_approx_cut){
                    // GxE with main effect
                    Eigen::Map<Eigen::VectorXd> snpk(snp_mat_part.data() + k * num_id_used, num_id_used);
                    ScalarProjectionCache snp_main_projection = build_scalar_projection_cache(this->m_Vi0, snpk);
                    VectorXd projected_null_response = apply_projection(
                        this->m_Vi0,
                        snp_main_projection.inverse_times_design,
                        snp_main_projection.design_cross_inverse,
                        m_y).col(0);
                    
                    if(speed == 0){
                        // exact test
                        MatrixXd EoG = this->m_bye_mat.array().colwise() * snpk.array();
                        VectorXd EoGtPy = EoG.transpose() * projected_null_response;
                        score = EoGtPy.transpose() * EoGtPy;
                        MatrixXd EoGtPEoG = EoG.transpose() * apply_projection(
                            this->m_Vi0,
                            snp_main_projection.inverse_times_design,
                            snp_main_projection.design_cross_inverse,
                            EoG);
                        
                        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
                        eigensolver.compute(EoGtPEoG);
                        chi_coef_Vec_kth = eigensolver.eigenvalues().real();
                        p_gxe = saddle(score, chi_coef_Vec_kth);
                    }else{
                        // approximate test
                        VectorXd EtGoPy = this->m_bye_mat.transpose() * (snpk.cwiseProduct(projected_null_response));
                        score = EtGoPy.transpose() * EtGoPy;
                        chi_coef_Vec_kth = chi2_coef_GxEwithMain_Vec * snpk_norm2;
                        p_gxe = saddle(score, chi_coef_Vec_kth);
                    }
                    
                }
                score_Vec(k) = score;
                p_gxe_Vec(k) = p_gxe;
            }else{
                score_Vec(k) = NAN;
                p_gxe_Vec(k) = NAN;
            }
        }

        // Output
        vector <std::string> snp_info_vec = GenoA.snp_anno(start_snp, num_snp_read);
        fastGxE::output_GxE_multi(fout, snp_info_vec, freq_arr, missing_rate_arr,
            beta_main_Vec, se_main_Vec, p_main_Vec, score_Vec, p_gxe_Vec, start_snp, num_snp_read);
    }

    fout.close();
}

void fastGxE::output_GxE_multi(std::ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
            const VectorXd& beta_main_Vec, const VectorXd& se_main_Vec, const VectorXd& p_main_Vec,
            const VectorXd& score_Vec, const VectorXd& p_Vec, std::int64_t start_snp, std::int64_t num_snp_read){
    for(std::int64_t i = 0; i < num_snp_read; i++){
        fout << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) << " "
            << beta_main_Vec(i) << " " << se_main_Vec(i) << " " << p_main_Vec(i) << " "
            << score_Vec(i) << " " << p_Vec(i) << std::endl;
    }
}

void fastGxE::test_GxE(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                int speed, std::int64_t num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut){
    (void)p_cut;

    spdlog::info("GxE analysis for individual environments");
    GENO GenoA(bed_file);
    const BedScanContext bed_context = prepare_bed_scan_context(GenoA, m_id_in_data_vec);
    spdlog::info("The number of used iids: {}", bed_context.num_id_used);
    const VectorXd Vi0y = m_Vi0 * m_y;
    
    spdlog::info("Randomly select SNPs and calculate the gamma");
    const vector<std::int64_t> random_snp_index_vec = sample_random_snp_indices(bed_context.num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
    load_centered_random_snp_panel(
        GenoA, bed_file, random_snp_index_vec, bed_context.id_index_in_bed_vec,
        snp_mat_part_random, freq_arr, missing_rate_arr, nobs_geno_arr);
    const VectorXd snp_squaredNorm_Vec = snp_mat_part_random.colwise().squaredNorm();

    // 1) Random-SNP calibration for the SNP main effect.
    spdlog::info("main");
    VectorXd p_main_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd beta_main_var_Vec = 1 / ((snp_mat_part_random).cwiseProduct(m_Vi0 * snp_mat_part_random)).colwise().sum().array();
    VectorXd se_main_Vec = beta_main_var_Vec.array().sqrt();
    VectorXd beta_main_Vec = beta_main_var_Vec.cwiseProduct(snp_mat_part_random.transpose() * Vi0y);
    VectorXd gamma_main_Vec = 1 / (beta_main_var_Vec.array() * snp_squaredNorm_Vec.array());  // gamma = x'V(-1)x / x'x

    std::int64_t num_random_snp_used = 0;
    double gamma_main = 0;
    for(std::int64_t i = 0; i < num_random_snp; i++){
        if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
            num_random_snp_used++;
            gamma_main += gamma_main_Vec[i];
            double beta = beta_main_Vec[i];
            double se = se_main_Vec[i];
            p_main_Vec[i] = gsl_cdf_chisq_Q(beta * beta / (se * se), 1);
        }else{
            p_main_Vec[i] = NAN;
            se_main_Vec[i] = NAN;
            beta_main_Vec[i] = NAN;
            gamma_main_Vec[i] = NAN;
        }
    }
    require_random_snp_fraction(num_random_snp_used, num_random_snp, 0.5, "Used random SNPs: {} (<50%)");
    gamma_main /= num_random_snp_used;

    ofstream fout = open_output_stream(m_out_file + ".main.random.res");
    
    fout << "order af missing gamma beta se p" << std::endl;
    for(std::int64_t i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) << " "
              << gamma_main_Vec[i] << " " << beta_main_Vec[i] << " " << se_main_Vec[i] << " " << p_main_Vec[i]
               << std::endl;
    }
    fout.close();
    gamma_main_Vec.resize(0); beta_main_Vec.resize(0); se_main_Vec.resize(0); p_main_Vec.resize(0); beta_main_var_Vec.resize(0);

    // 2) Random-SNP calibration for marginal GxE terms without the SNP main effect.
    spdlog::info("GxE without main");
    VectorXd gamma_GxEnoMain_Vec = VectorXd::Zero(m_num_bye);
    MatrixXd beta_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    MatrixXd se_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    MatrixXd p_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    MatrixXd gamma_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    for(std::int64_t j = 0; j < m_num_bye; j++){
        MatrixXd snp_mat_part_random_xE = snp_mat_part_random.array().colwise() * m_bye_mat.col(j).array();
        VectorXd beta_GxEnoMain_var_Vec = 1 / ((snp_mat_part_random_xE).cwiseProduct(m_Vi0 * snp_mat_part_random_xE)).colwise().sum().array();
        VectorXd se_GxEnoMain_Vec = beta_GxEnoMain_var_Vec.array().sqrt();
        VectorXd beta_GxEnoMain_Vec = beta_GxEnoMain_var_Vec.cwiseProduct(snp_mat_part_random_xE.transpose() * Vi0y);
        VectorXd gamma_GxEnoMain_j = 1 / (beta_GxEnoMain_var_Vec.array() * snp_squaredNorm_Vec.array());
        VectorXd p_GxEnoMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
        for(std::int64_t i = 0; i < num_random_snp; i++){
            if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
                gamma_GxEnoMain_Vec[j] += gamma_GxEnoMain_j[i];
                double beta = beta_GxEnoMain_Vec[i];
                double se = se_GxEnoMain_Vec[i];
                p_GxEnoMain_Vec[i] = gsl_cdf_chisq_Q(beta * beta / (se * se), 1);
            }else{
                p_GxEnoMain_Vec[i] = NAN;
                se_GxEnoMain_Vec[i] = NAN;
                beta_GxEnoMain_Vec[i] = NAN;
                gamma_GxEnoMain_j[i] = NAN;
            }
        }
        beta_GxEnoMain_Mat.col(j) = beta_GxEnoMain_Vec;
        se_GxEnoMain_Mat.col(j) = se_GxEnoMain_Vec;
        p_GxEnoMain_Mat.col(j) = p_GxEnoMain_Vec;
        gamma_GxEnoMain_Mat.col(j) = gamma_GxEnoMain_j;
    }
    gamma_GxEnoMain_Vec = gamma_GxEnoMain_Vec / num_random_snp_used;

    fout = open_output_stream(m_out_file + ".GxEnoMain.random.res");
    fout << "order af missing";
    for(std::int64_t i = 0; i < m_num_bye; i++){
        fout << " gamma" << i + 1 << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << std::endl;

    for(std::int64_t i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        for(std::int64_t j = 0; j < m_num_bye; j++){
            fout <<  " " << gamma_GxEnoMain_Mat(i, j) << " " << beta_GxEnoMain_Mat(i, j) << " " 
            << se_GxEnoMain_Mat(i, j) << " " << p_GxEnoMain_Mat(i, j);
        }
        fout << std::endl;
    }
    fout.close();
    gamma_GxEnoMain_Mat.resize(0, 0); beta_GxEnoMain_Mat.resize(0, 0); se_GxEnoMain_Mat.resize(0, 0); p_GxEnoMain_Mat.resize(0, 0);

    // 3) Random-SNP calibration for joint SNP + single-environment fits.
    spdlog::info("GxE with main");
    MatrixXd beta_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    MatrixXd se_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    MatrixXd p_Mat = MatrixXd::Ones(num_random_snp, m_num_bye) * NAN;
    vector<vector<MatrixXd>> gamma_Mat_2Dvec;
    gamma_Mat_2Dvec.reserve(m_num_bye);
    for(std::int64_t j = 0; j < m_num_bye; j++){
        vector<MatrixXd> gamma_Mat_vec;
        gamma_Mat_vec.reserve(num_random_snp);
        for(std::int64_t i = 0; i < num_random_snp; i++){
            gamma_Mat_vec.push_back(MatrixXd::Zero(2, 2));
        }
        gamma_Mat_2Dvec.push_back(gamma_Mat_vec);
    }

    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t i = 0; i < num_random_snp; i++){
        if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
            VectorXd snpi = snp_mat_part_random.col(i);
            double snp_norm2 = snp_squaredNorm_Vec(i);
            MatrixXd snpiME = m_bye_mat.array().colwise() * snpi.array();
            MatrixXd xmati = MatrixXd::Zero(bed_context.num_id_used, 2);
            xmati.col(0) = snpi;
            for(std::int64_t j = 0; j < m_num_bye; j++){
                xmati.col(1) = snpiME.col(j);
                MatrixXd design_cross_product = xmati.transpose() * m_Vi0 * xmati;
                MatrixXd design_cross_inverse = design_cross_product.inverse();
                VectorXd beta_Vec = design_cross_inverse * (xmati.transpose() * Vi0y);
                double beta_var = design_cross_inverse(1, 1);
                double beta = beta_Vec(1);
                double p = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);
                beta_Mat(i, j) = beta;
                se_Mat(i, j) = sqrt(beta_var);
                p_Mat(i, j) = p;
                gamma_Mat_2Dvec[j][i] = design_cross_product / snp_norm2;
            }
        }
    }

    vector<MatrixXd> gamma_Mat_vec;
    gamma_Mat_vec.reserve(m_num_bye);
    for(std::int64_t j = 0; j < m_num_bye; j++){
        MatrixXd gamma_Mat = MatrixXd::Zero(2, 2);
        for(std::int64_t i = 0; i < num_random_snp; i++){
            gamma_Mat += gamma_Mat_2Dvec[j][i];
        }
        gamma_Mat_vec.push_back(gamma_Mat / num_random_snp_used);
    }

    
    fout = open_output_stream(m_out_file + ".GxE.random.res");
    fout << "order af missing";
    for(std::int64_t i = 0; i < m_num_bye; i++){
        fout << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << std::endl;
    for(std::int64_t i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        for(std::int64_t j = 0; j < m_num_bye; j++){
            fout << " " << beta_Mat(i, j) << " " << se_Mat(i, j) << " " << p_Mat(i, j);
        }
        fout << std::endl;
    }
    fout.close();

    // The environment correlation is used as a lightweight approximation for
    // the multi-environment combined score test in the genome-wide scan.
    MatrixXd corrE = computeCorrelationMatrix(m_bye_mat);
    corrE.diagonal() += Eigen::VectorXd::Constant(corrE.rows(), kCorrMatrixNudge);

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(corrE);
    if (eigensolver.info() != Eigen::Success) {
        fatal_error("Fail to eigen decomposition of GRM");
    }
    
    VectorXd corrE_eigenvalues = eigensolver.eigenvalues();
    
    // 4) Genome-wide scan using the random-SNP calibration factors above.
    spdlog::info("Genome-wide association");
    fout = open_output_stream(m_out_file + ".res");
    fout << "order chrom SNP cm base allele1 allele2 af missing beta_main se_main p_main";
    for(std::int64_t i = 0; i < m_num_bye; i++){
        fout << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << " p_single p_multi p_gxe" << std::endl;

    const std::int64_t initial_partition_size = get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read;
    std::vector<double> snp_mat_part(initial_partition_size * bed_context.num_id_used);
    VectorXd geno_var_Vec = VectorXd::Zero(initial_partition_size);
    beta_Mat.resize(initial_partition_size, m_num_bye + 1);
    se_Mat.resize(initial_partition_size, m_num_bye + 1);
    p_Mat.resize(initial_partition_size, m_num_bye + 1);
    MatrixXd p_combined_Mat = MatrixXd::Zero(initial_partition_size, 3);
    VectorXd EVy = VectorXd::Zero(bed_context.num_id_used);
    VectorXd prod_Vec = VectorXd::Zero(initial_partition_size);
    
    for(std::int64_t ipart = 0; ipart < npart_snp; ipart++){
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, ipart);
        if(partition.num_snp_read != geno_var_Vec.size()){
            snp_mat_part.resize(partition.num_snp_read * bed_context.num_id_used);
            geno_var_Vec.resize(partition.num_snp_read);
            beta_Mat.resize(partition.num_snp_read, m_num_bye + 1);
            se_Mat.resize(partition.num_snp_read, m_num_bye + 1);
            p_Mat.resize(partition.num_snp_read, m_num_bye + 1);
            p_combined_Mat.resize(partition.num_snp_read, 3);
            prod_Vec.resize(partition.num_snp_read);
        }
        
        spdlog::info("Part {}/{}: Start SNP {}, read SNPs {}", ipart + 1, npart_snp, partition.start_snp, partition.num_snp_read);
        GenoA.read_bed_centered_to_buffer_omp(
            bed_file + ".bed", snp_mat_part.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
            partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);

        // First get cheap marginal approximations for the main effect and each
        // environment-specific interaction; fall back to exact per-SNP refits
        // only when one of those approximations looks promising.
        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t j = 0; j < partition.num_snp_read; ++j) {
            geno_var_Vec(j) = cblas_ddot(
                bed_context.num_id_used, snp_mat_part.data() + j * bed_context.num_id_used, 1,
                snp_mat_part.data() + j * bed_context.num_id_used, 1);
        }

        // main
        se_Mat.col(0) = 1.0 / (gamma_main * geno_var_Vec.array());

        cblas_dgemv(CblasRowMajor, CblasNoTrans, partition.num_snp_read, bed_context.num_id_used, 1.0, snp_mat_part.data(), 
                bed_context.num_id_used, Vi0y.data(), 1, 0.0, prod_Vec.data(), 1);
        vdMul(partition.num_snp_read, se_Mat.col(0).data(), prod_Vec.data(), beta_Mat.col(0).data());
        
        // GxE
        for(std::int64_t k = 0; k < m_num_bye; k++){
            se_Mat.col(k+1) = 1.0 / (gamma_GxEnoMain_Vec(k) * geno_var_Vec.array());
            vdMul(bed_context.num_id_used, m_bye_mat.col(k).data(), Vi0y.data(), EVy.data());
            cblas_dgemv(CblasRowMajor, CblasNoTrans, partition.num_snp_read, bed_context.num_id_used, 1.0, snp_mat_part.data(), 
                bed_context.num_id_used, EVy.data(), 1, 0.0, prod_Vec.data(), 1);
            vdMul(partition.num_snp_read, se_Mat.col(k+1).data(), prod_Vec.data(), beta_Mat.col(k+1).data());
        }
        se_Mat = se_Mat.array().sqrt();
        MatrixXd z_mat = beta_Mat.array() / se_Mat.array();

        #pragma omp parallel for schedule(dynamic)
        for(std::int64_t k = 0; k < partition.num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                VectorXd p_Vec = VectorXd::Ones(m_num_bye + 1) * NAN;

                for(std::int64_t m = 0; m < m_num_bye + 1; m++){
                    double z = z_mat(k, m);
                    double chi2_val = z * z;
                    p_Vec(m) = gsl_cdf_chisq_Q(chi2_val, 1);
                }

                if(p_Vec.minCoeff() < p_approx_cut){
                    Eigen::Map<Eigen::VectorXd> snpk(
                        snp_mat_part.data() + k * bed_context.num_id_used, bed_context.num_id_used);
                    double snpk_norm2 = geno_var_Vec(k);
                    MatrixXd snpkME = m_bye_mat.array().colwise() * snpk.array();
                    MatrixXd xmatk = MatrixXd::Zero(bed_context.num_id_used, 2);
                    xmatk.col(0) = snpk;

                    for(std::int64_t m = 0; m < m_num_bye; m++){
                            xmatk.col(1) = snpkME.col(m);
                            if(speed == 0){
                                // exact test
                                MatrixXd design_cross_product = xmatk.transpose() * m_Vi0 * xmatk;
                                MatrixXd design_cross_inverse = design_cross_product.inverse();
                                VectorXd beta_Vec = design_cross_inverse * (xmatk.transpose() * Vi0y);
                                double beta_var = design_cross_inverse(1, 1);
                                double beta = beta_Vec(1);
                                double p_val = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);
                                beta_Mat(k, m+1) = beta;
                                se_Mat(k, m+1) = sqrt(beta_var);
                                z_mat(k, m+1) = beta_Mat(k, m+1) / se_Mat(k, m+1);
                                p_Vec(m+1) = p_val;
                            }else{
                                // approximate test
                                MatrixXd beta_Vec_var = (gamma_Mat_vec[m] * snpk_norm2).inverse();
                                double beta = (beta_Vec_var * (xmatk.transpose() * Vi0y))(1);
                                double beta_var = beta_Vec_var(1, 1);
                                double chi_val = beta * beta / beta_var;
                                double p_val = gsl_cdf_chisq_Q(chi_val, 1);
                                beta_Mat(k, m+1) = beta;
                                se_Mat(k, m+1) = sqrt(beta_var);
                                z_mat(k, m+1) = beta_Mat(k, m+1) / se_Mat(k, m+1);
                                p_Vec(m+1) = p_val;
                            }
                    }

                }

                p_Mat.row(k) = p_Vec;
                Eigen::VectorXd result = acat_combine_pvalues(p_Vec.tail(m_num_bye));
                p_combined_Mat(k, 0) = result(1);

            }else{
                beta_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                se_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                p_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                p_combined_Mat(k, 0) = NAN;
            }
        }

        VectorXd score_Vec = z_mat.rightCols(m_num_bye).array().square().matrix().rowwise().sum();
        
        #pragma omp parallel for schedule(dynamic)
        for(std::int64_t k = 0; k < partition.num_snp_read; k++){
            double score_val = score_Vec(k);
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut) && score_val >= 0.0 && std::isfinite(score_val)){

                double p1 = saddle(score_val, corrE_eigenvalues);
                p_combined_Mat(k, 1) = p1;

                VectorXd p_Vec(2);
                p_Vec(0) = p_combined_Mat(k, 0);
                p_Vec(1) = p1;
                p_combined_Mat(k, 2) = acat_combine_pvalues(p_Vec)(1);
            }else{
                p_combined_Mat(k, 1) = p_combined_Mat(k, 2) = NAN;
            }
        }
        
        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        fastGxE::output_GxE(fout, snp_info_vec, freq_arr, missing_rate_arr,
            beta_Mat, se_Mat, p_Mat, p_combined_Mat,
            partition.start_snp, partition.num_snp_read);
    }
    fout.close();
}

void fastGxE::output_GxE(ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
            const MatrixXd& beta_mat, const MatrixXd& se_mat, const MatrixXd& p_mat, const MatrixXd& p_combined_Mat,
            std::int64_t start_snp, std::int64_t num_snp_read){
    
    std::stringstream ss;
    ss.precision(8);

    for (std::int64_t i = 0; i < num_snp_read; i++) {
        ss << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        
        for (std::int64_t m = 0; m < beta_mat.cols(); m++) {
            ss << " " << beta_mat(i, m) << " " << se_mat(i, m) << " " << p_mat(i, m);
        }

        for (std::int64_t m = 0; m < p_combined_Mat.cols(); m++) {
            ss << " " << p_combined_Mat(i, m);
        }

        ss << std::endl;
    }

    fout << ss.str();
}


int fastGxE::run(int argc, char **argv) {
    CLI::App app{"fastGxE - Scalable and fast multivariate GxE analysis"};
    RunOptions options;
    configure_run_cli(app, options);
    CLI11_PARSE(app, argc, argv);
    validate_run_options(options);
    log_run_options(options);
    
    // Set MKL and OpenMP threads
    mkl_set_num_threads(options.threads);
    omp_set_num_threads(options.threads);

    // Resolve all requested column names against the data-file header before any
    // expensive GRM or genotype work starts.
    ResolvedInputColumns resolved = resolve_input_columns(options);
    m_trait_name_vec = options.trait_vec;

    // prepare data
    this->pre_data(options.out_file, options.data_file, options.agrm_file,
            resolved.covariate_index_vec, resolved.class_index_vec,
            resolved.bye_index_vec, resolved.trait_index_vec, options.missing_in_data_vec);
    
    if (options.test_main_binary) {
        require_binary_trait_matrix(m_y);
    } else if (options.test_main_binary_continuous) {
        require_binary_continuous_trait_matrix(m_y);
    } else if (options.test_main_multitrait_continuous) {
        require_multitrait_continuous_trait_matrix(m_y);
    }

    VectorXd init_varcom;
    this->process_grm(options.test_main);
    if(options.test_main){
        this->varcom_main(init_varcom, options.maxiter, options.cc_gra, options.cc_gra, options.cc_logL);
    }else if(options.test_main_binary){
        this->varcom_main_binary(init_varcom, options.maxiter, options.cc_par, options.cc_gra, options.cc_logL);
    }else if(options.test_main_binary_continuous){
        this->varcom_main_binary_continuous(init_varcom, options.maxiter, options.cc_par, options.cc_gra, options.cc_logL);
    }else if(options.test_main_multitrait_continuous){
        this->varcom_main_multitrait_continuous(init_varcom, options.maxiter, options.cc_par, options.cc_gra, options.cc_logL);
    }else if(options.test_gxe){
        this->pre_data_GxE(options.standardize_env, true);  // Phenotype correction is always required here.
        this->varcom_GxE(init_varcom, options.no_noisebye, options.maxiter, options.cc_par, options.cc_gra, options.cc_logL);
    }
    
    
    if(options.bed_file.empty()){
        spdlog::info("If you want to perform GWAS, please give --bfile.");
        return 0;
    }
    
    vector<std::int64_t> tmp_vec = this->get_start_end_pos(options.bed_file, options.snp_range_vec, options.split_task_vec);
    std::int64_t start_pos = tmp_vec[0], end_pos = tmp_vec[1];
    if(end_pos - start_pos <= options.npart_snp) options.npart_snp = 1;

    if(options.test_main){
        this->test_main(options.bed_file, start_pos, end_pos, options.npart_snp,
                options.speed, options.num_random_snp, options.p_approx_cut, options.p_cut, options.maf_cut, options.missing_rate_cut,
                options.maxiter, options.cc_par, options.cc_gra);
    }else if(options.test_main_binary){
        this->test_main_binary(options.bed_file, start_pos, end_pos, options.npart_snp,
                options.speed, options.num_random_snp, options.p_approx_cut, options.p_cut,
                options.maf_cut, options.missing_rate_cut, options.maxiter, options.cc_par, options.cc_gra);
    }else if(options.test_main_binary_continuous){
        this->test_main_binary_continuous(options.bed_file, start_pos, end_pos, options.npart_snp,
                options.speed, options.num_random_snp, options.p_approx_cut, options.p_cut,
                options.maf_cut, options.missing_rate_cut, options.maxiter, options.cc_par, options.cc_gra);
    }else if(options.test_main_multitrait_continuous){
        this->test_main_multitrait_continuous(options.bed_file, start_pos, end_pos, options.npart_snp,
                options.speed, options.num_random_snp, options.p_approx_cut, options.p_cut,
                options.maf_cut, options.missing_rate_cut, options.maxiter, options.cc_par, options.cc_gra);
    }else if(options.test_gxe){
        this->test_GxE(options.bed_file, start_pos, end_pos, options.npart_snp,
                                 options.speed, options.num_random_snp, options.p_approx_cut, options.p_cut,
                                 options.maf_cut, options.missing_rate_cut);
    }
    return 0;
}


vector<std::int64_t> fastGxE::get_start_end_pos(const string& bed_file,
            const vector<string>& snp_range_vec, const vector<std::int64_t>& split_task_vec){
    GENO GenoA(bed_file);
    std::int64_t num_snp = GenoA.get_num_sid();
    std::int64_t start_pos = 0, end_pos = num_snp;
    if((!snp_range_vec.empty()) && (!split_task_vec.empty())){
        fatal_error("Options --snp-range and --split-task cannot be used together.");
    }else if(!snp_range_vec.empty()){
        vector<std::int64_t> tmp = GenoA.find_bim_index(snp_range_vec);
        start_pos = tmp[0];
        end_pos = tmp[1] + 1;
        this->reset_output_prefix(this->m_out_file, start_pos, end_pos);
    }else if(!split_task_vec.empty()){
        if(split_task_vec[0] < split_task_vec[1] || split_task_vec[1] <= 0){
            fatal_error("The second number in --split-task must be greater than 0 and not exceed the first number.");
        }
        std::int64_t num_snp_part = num_snp / split_task_vec[0];
        start_pos = num_snp_part * (split_task_vec[1] - 1);
        end_pos = num_snp_part * split_task_vec[1];
        if(split_task_vec[0] == split_task_vec[1]) end_pos = num_snp;
        this->reset_output_prefix(this->m_out_file, split_task_vec[0], split_task_vec[1]);
    }
    if(start_pos >= end_pos){
        fatal_error("The start SNP position must not exceed the end SNP position.");
    }
    vector<std::int64_t> vec = {start_pos, end_pos};
    return vec;
}
