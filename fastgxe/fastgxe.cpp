/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:40
 * LastEditors: Chao Ning
 * LastEditTime: 2026-04-11 20:42:12
 */

// ============================================================================
// fastGxE - Scalable genome-wide GxE / SNP main-effect association via LMM
// ----------------------------------------------------------------------------
// This file implements the `fastGxE` class together with a set of stateless
// helper functions in an anonymous namespace. fastGxE is a statistical-genetics
// tool that runs genome-wide association scans under a linear mixed model (LMM).
//
// Two mutually exclusive analysis modes are supported:
//   --test-main : a SNP main-effect association scan (one variance component
//                 for the polygenic genetic effect plus residual noise).
//   --test-gxe  : a scalable and effective genome-wide multi-environment GxE method.
//
// High-level pipeline (see fastGxE::run):
//   1. Parse CLI options (configure_run_cli / validate_run_options) and resolve
//      requested column names against the data-file header (resolve_input_columns).
//   2. pre_data: load phenotypes, restrict to samples present in the GRM, then
//      REORDER individuals by their relatedness group. After reordering the GRM
//      is block-diagonal and stored sparsely:
//        - block 0 = unrelated singletons, stored as a diagonal (N x 1 vector);
//        - blocks 1.. = dense small within-family kinship blocks.
//   3. process_grm: assemble the sparse GRM. For main-effect mode it also
//      eigen-decomposes the GRM (block-wise) so the null model becomes diagonal
//      algebra in the GRM eigenbasis.
//   4. Fit null-model variance components by AI-REML:
//        - varcom_main  : 2 components via the eigendecomposition.
//        - varcom_GxE   : up to 4 components (genetic / GxE-relationship /
//                         noise-by-environment / residual).
//      Cache V^{-1} (m_Vi0), the null log-likelihood (m_logL_null) and the
//      estimated components (m_varcom_null).
//   5. Per-SNP fast score tests reading the PLINK .bed in chunks
//      (test_main, test_GxE, test_GxE_multi). The GxE scans calibrate
//      approximate test statistics using a random panel of SNPs ("gamma"
//      calibration), and combine per-environment GxE p-values with ACAT.
//
// Recurring numerical helpers (anonymous namespace):
//   build_inverse_covariance_matrix - builds block-diagonal V^{-1} and log|V|.
//   build_projection_cache/apply_projection - evaluate the fixed-effect
//      projection P = V^{-1} - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1} without ever
//      forming the dense n x n matrix P.
//   find_positive_variance_update - AI/EM blended REML step that stays in the
//      positive-variance region.
//   build_gxe_relationship_data - environment-weighted GRM for the GxE model.
// ============================================================================

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

// ===== Anonymous-namespace helpers (CLI, linear algebra, GRM, REML, scan) =====
namespace {

// Numerical constants for internal algorithms.
// These are intentionally not CLI-configurable: changing them would require
// re-validation of statistical calibration.
constexpr double kCorrMatrixNudge        = 0.001; // ridge added to environment correlation matrix

// ===== CLI / option handling and header parsing =====

// Reads only the first (header) line of the data file and splits it into column
// names. Aborts via fatal_error if the file cannot be opened or the header is
// empty. Returns the vector of column-name tokens.
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

// Joins a list of strings with `delimiter`, or returns the literal "None" when
// the list is empty. Used purely for human-readable logging of options.
std::string join_or_none(const std::vector<std::string>& values,
                         const char* delimiter = ", ") {
    return values.empty() ? "None" : join_string(values, delimiter);
}

// All command-line options parsed from argv, with their default values. A single
// instance is filled by configure_run_cli and then validated/logged.
struct RunOptions {
    int threads = 10;
    bool test_main = false;
    bool test_gxe = false;
    bool standardize_env = true;
    bool no_noisebye = false;
    bool log_transform = false;
    bool inv_normal = false;
    bool keep_random = false;

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

    std::int64_t num_random_snp = 1000;
    int npart_snp = 20;
    double exact_cut = 1.0e-5;  // GxE: exact-refit when approximate p < this (0=off, 1=all)
    int maxiter = 100;
    double p_approx_cut = 1.0e-3;
    double maf_cut = 0.01;
    double missing_rate_cut = 0.05;
    double cc_par = 1.0e-7;
    double cc_gra = 1.0e-6;
    double cc_logL = 5.0e-5;
};

// ===== Internal data structures (option resolution and numerical caches) =====

// Column indices (into the data-file header) resolved from the user-supplied
// trait / covariate / class / interacting-environment names.
struct ResolvedInputColumns {
    std::vector<std::int64_t> trait_index_vec;
    std::vector<std::int64_t> covariate_index_vec;
    std::vector<std::int64_t> class_index_vec;
    std::vector<std::int64_t> bye_index_vec;
};

// The environment-weighted ("GxE") relationship matrix, held both as per-group
// dense blocks and as the assembled sparse block-diagonal matrix.
struct GxeRelationshipData {
    std::vector<MatrixXd> group_blocks;
    SparseMatrix<double> sparse_matrix;
};

// Result of inverting the block-diagonal covariance V: the sparse inverse and
// the accumulated log-determinant log|V|.
struct InverseCovarianceResult {
    SparseMatrix<double> inverse_matrix;
    double logdet = 0.0;
};

// Precomputed quantities for the fixed-effect projection with a matrix design X:
//   inverse_times_design = V^{-1} X
//   design_cross_inverse = (X' V^{-1} X)^{-1}
//   design_logdet        = log|X' V^{-1} X| (optional).
struct ProjectionCache {
    MatrixXd inverse_times_design;
    MatrixXd design_cross_inverse;
    double design_logdet = 0.0;
};

// Same projection cache but specialized to a single (vector) design column, so
// the cross term X'V^{-1}X is just a scalar.
struct ScalarProjectionCache {
    VectorXd inverse_times_design;
    double design_cross_inverse = 0.0;
};

// Result of one REML variance-component update: the step `delta` and the
// resulting `updated_varcom = varcom + delta`.
struct VarianceUpdateResult {
    VectorXd delta;
    VectorXd updated_varcom;
};

// Per-iteration REML quantities for the 2-component main-effect model:
// the -2 log-likelihood, the score vector (fd_mat) and the AI matrix.
struct MainModelIterationState {
    double logL = 0.0;
    VectorXd fd_mat = VectorXd::Zero(2);
    MatrixXd ai_mat = MatrixXd::Zero(2, 2);
};

// Cached metadata for scanning a PLINK .bed file: the mapping from analysis
// samples to .fam rows, and the counts of used samples and SNPs.
struct BedScanContext {
    std::vector<std::int64_t> id_index_in_bed_vec;
    std::int64_t num_id_used = 0;
    std::int64_t num_snp = 0;
    std::int64_t num_id_in_bed = 0;
};

// A logical slice of the SNP range to read in one streaming chunk: the starting
// SNP index and how many SNPs to read.
struct SnpPartition {
    std::int64_t start_snp = 0;
    std::int64_t num_snp_read = 0;
};

// Work buffers for scalar (single-trait) SNP scans used by test_main.
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

using GrmTriplet = Eigen::Triplet<double>;

// ===== GRM construction (block-diagonal sparse assembly) =====

// Computes the global row/column start offset of each GRM block (a prefix sum of
// block row counts). offsets[k] is where block k begins; offsets.back() == N.
std::vector<std::int64_t> build_grm_block_offsets(const std::vector<MatrixXd>& grm_blocks) {
    std::vector<std::int64_t> offsets;
    offsets.reserve(grm_blocks.size() + 1);
    offsets.push_back(0);

    for (const auto& block : grm_blocks) {
        offsets.push_back(offsets.back() + block.rows());
    }

    return offsets;
}

// Counts how many sparse triplets the block-diagonal GRM produces, so callers can
// reserve() once. Block 0 (singletons) contributes only its diagonal; every other
// block contributes its full dense entries.
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

// Appends the diagonal triplets for the singleton block (block 0), which is stored
// as an N x 1 column vector. With use_identity=true the diagonal is forced to 1
// (used to build the orthonormal rotation for unrelated samples in the eigenbasis);
// otherwise the stored value is used.
void append_singleton_block_triplets(const MatrixXd& singleton_block, std::int64_t start_index,
        std::vector<GrmTriplet>& triplets, bool use_identity) {
    const std::int64_t block_size = singleton_block.rows();
    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        const double value = use_identity ? 1.0 : singleton_block(row_idx, 0);
        triplets.emplace_back(start_index + row_idx, start_index + row_idx, value);
    }
}

// Appends every entry of a dense within-group block as a sparse triplet, shifted
// to its global position by start_index.
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

// For the singleton (unrelated) block, the GRM is already diagonal, so its
// "eigenvalues" are simply the diagonal entries; copy them into the global
// eigenvalue vector at the block offset.
void load_singleton_block_eigenvalues(const MatrixXd& singleton_block, std::int64_t start_index,
        VectorXd& eigenvalues) {
    const std::int64_t block_size = singleton_block.rows();
    for (std::int64_t row_idx = 0; row_idx < block_size; ++row_idx) {
        eigenvalues(start_index + row_idx) = singleton_block(row_idx, 0);
    }
}

// Eigen-decomposes one dense within-group GRM block. Writes the block's
// eigenvalues into the global eigenvalue vector, and appends its eigenvectors as
// sparse triplets so they assemble into the global block-diagonal rotation Q.
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

// ===== Numerical / linear-algebra helpers =====

// Computes a column-pivoting QR of `design_matrix` and aborts (fatal_error with
// `error_message`) if its rank is below `expected_rank`, i.e. if the design has
// linearly dependent columns. Returns the QR factorization for reuse.
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

// Standardizes each column of the environment matrix in place to mean 0 and
// (population) standard deviation 1.
void standardize_environment_matrix(MatrixXd& environment_matrix) {
    const VectorXd mean_vector = environment_matrix.colwise().mean();
    MatrixXd centered_matrix = environment_matrix.rowwise() - mean_vector.transpose();
    const VectorXd std_vector =
        (centered_matrix.cwiseProduct(centered_matrix).colwise().sum().array() / centered_matrix.rows()).sqrt();
    environment_matrix = centered_matrix.array().rowwise() * (1.0 / std_vector.transpose().array());
}

// Natural-log transform of the phenotype matrix in place (applied per column).
// Requires strictly positive values; aborts otherwise so the user is forced to
// fix or shift the trait rather than silently producing NaN/-inf.
void apply_log_transform(MatrixXd& phenotype) {
    if ((phenotype.array() <= 0.0).any()) {
        fatal_error("--log-transform requires all phenotype values to be strictly positive; "
                    "found a value <= 0. Shift the trait or drop non-positive samples first.");
    }
    phenotype = phenotype.array().log().matrix();
}

// Rank-based inverse-normal transform of a single vector. Ties get the average
// (mid) rank; the Blom offset c = 3/8 is used, i.e.
//   transformed_i = Phi^{-1}( (rank_i - 3/8) / (n + 1/4) ),
// which maps the values onto standard-normal quantiles while preserving order.
VectorXd rank_based_inverse_normal(const VectorXd& values) {
    const std::int64_t n = values.size();
    if (n == 0) {
        return values;
    }

    // Order indices by value so ties can be detected as equal-valued runs.
    std::vector<std::int64_t> order(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) {
        order[static_cast<std::size_t>(i)] = i;
    }
    std::sort(order.begin(), order.end(),
              [&values](std::int64_t a, std::int64_t b) { return values(a) < values(b); });

    // Assign 1-based ranks, averaging within each tie group.
    VectorXd ranks(n);
    std::int64_t i = 0;
    while (i < n) {
        std::int64_t j = i;
        while (j + 1 < n && values(order[static_cast<std::size_t>(j + 1)]) ==
                            values(order[static_cast<std::size_t>(i)])) {
            ++j;
        }
        const double average_rank = (static_cast<double>(i) + static_cast<double>(j)) / 2.0 + 1.0;
        for (std::int64_t k = i; k <= j; ++k) {
            ranks(order[static_cast<std::size_t>(k)]) = average_rank;
        }
        i = j + 1;
    }

    const double c = 3.0 / 8.0;  // Blom offset
    VectorXd transformed(n);
    for (std::int64_t k = 0; k < n; ++k) {
        const double quantile = (ranks(k) - c) / (static_cast<double>(n) - 2.0 * c + 1.0);
        transformed(k) = gsl_cdf_ugaussian_Pinv(quantile);
    }
    return transformed;
}

// Applies the rank-based inverse-normal transform to each column of the
// phenotype matrix in place.
void apply_inverse_normal_transform(MatrixXd& phenotype) {
    for (std::int64_t col = 0; col < phenotype.cols(); ++col) {
        phenotype.col(col) = rank_based_inverse_normal(phenotype.col(col));
    }
}

// Appends `extra_columns` to the right of `design_matrix` in place, growing the
// column count while preserving the existing columns.
void append_columns(MatrixXd& design_matrix, const MatrixXd& extra_columns) {
    const std::int64_t original_col_count = design_matrix.cols();
    design_matrix.conservativeResize(design_matrix.rows(), original_col_count + extra_columns.cols());
    design_matrix.rightCols(extra_columns.cols()) = extra_columns;
}

// Formats a vector as a single space-separated line for compact logging.
std::string format_vector_for_log(const VectorXd& values) {
    const Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "", "", "", "", "");
    std::ostringstream oss;
    oss << values.transpose().format(fmt);
    return oss.str();
}

// Returns a valid starting point for REML: the caller's initial_values if they
// have the right length and are all strictly positive, otherwise a vector of ones.
VectorXd initialize_variance_components(const VectorXd& initial_values, std::int64_t expected_size) {
    if (initial_values.size() == expected_size && initial_values.minCoeff() > 0) {
        return initial_values;
    }
    return VectorXd::Ones(expected_size);
}

// Builds one block of the GxE relationship matrix: the elementwise product of the
// GRM block with the environment Gram matrix E E'/num_bye. For the singleton block
// only the diagonal (row-wise squared norm of E / num_bye) is needed.
MatrixXd build_gxe_relationship_block(const MatrixXd& grm_block, const MatrixXd& bye_block, std::int64_t num_bye,
        bool is_singleton_block) {
    if (is_singleton_block) {
        MatrixXd gxe_block = (bye_block.cwiseProduct(bye_block)).rowwise().sum() / num_bye;
        return gxe_block.cwiseProduct(grm_block);
    }

    MatrixXd gxe_block = bye_block * bye_block.transpose() / num_bye;
    return gxe_block.cwiseProduct(grm_block);
}

// Builds the full GxE relationship matrix (the "environment-weighted GRM") used
// as a random-effect covariance in the GxE model. Returns both the per-group
// dense blocks and the assembled sparse block-diagonal matrix. Empty blocks are
// passed through unchanged.
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

        // Slice the rows of the environment matrix belonging to this group, then
        // form its environment-weighted GRM block.
        const MatrixXd bye_block = bye_mat.middleRows(start_index, block_size);
        const MatrixXd gxe_block = build_gxe_relationship_block(
            grm_block, bye_block, num_bye, block_idx == 0);
        relationship_data.group_blocks[block_idx] = gxe_block;

        // Block 0 is diagonal (singletons); all later blocks are dense.
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

// Builds the diagonal of the noise-by-environment covariance: for each sample,
// the mean squared environment value (||E_i||^2 / num_env). This scales the
// heteroscedastic residual variance component in the GxE model.
VectorXd build_noise_by_environment_vector(const MatrixXd& bye_mat) {
    return bye_mat.rowwise().squaredNorm() / bye_mat.cols();
}

// Wraps a vector as an Eigen sparse diagonal matrix.
SparseMatrix<double> build_sparse_diagonal_matrix(const VectorXd& diagonal_values) {
    SparseMatrix<double> diagonal_matrix(diagonal_values.size(), diagonal_values.size());
    diagonal_matrix.reserve(Eigen::VectorXi::Constant(diagonal_values.size(), 1));
    for (std::int64_t i = 0; i < diagonal_values.size(); ++i) {
        diagonal_matrix.insert(i, i) = diagonal_values(i);
    }
    diagonal_matrix.makeCompressed();
    return diagonal_matrix;
}

// Precomputes the pieces of the fixed-effect projection
//   P = V^{-1} - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}
// for a sparse V^{-1}. Stores V^{-1}X and (X'V^{-1}X)^{-1}, factorized via LDLT
// (aborting if X'V^{-1}X is not factorizable), and optionally log|X'V^{-1}X|
// which contributes to the REML log-likelihood.
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

// Overload of build_projection_cache for a diagonal V^{-1} (the GRM eigenbasis
// case): the same projection pieces, but V^{-1}X reduces to a column-wise scaling
// and the cross product uses a Cholesky (CustomLLT) factorization.
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

// Scalar version of build_projection_cache for a single design column x:
// stores V^{-1}x and the scalar (x'V^{-1}x)^{-1}. Used when projecting out a
// single SNP main effect.
ScalarProjectionCache build_scalar_projection_cache(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& design_vector) {
    ScalarProjectionCache projection_cache;
    projection_cache.inverse_times_design = inverse_matrix * design_vector;
    projection_cache.design_cross_inverse =
        1.0 / design_vector.dot(projection_cache.inverse_times_design);
    return projection_cache;
}

// The apply_projection family all evaluate P*rhs using the cached projection
// pieces (V^{-1}X and (X'V^{-1}X)^{-1}) without ever forming the n x n matrix P.
// The overloads cover combinations of: sparse vs diagonal V^{-1}, matrix vs
// scalar fixed-effect cross-inverse, and vector vs matrix right-hand side.

// Evaluates P·rhs = V⁻¹rhs - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹rhs without forming the n×n projection matrix P.
VectorXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const VectorXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Same projection P·rhs with a matrix right-hand side (columns projected jointly).
MatrixXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const MatrixXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Scalar-design variant: a single fixed-effect column x with scalar (x'V^{-1}x)^{-1}.
VectorXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& inverse_times_design, double design_cross_inverse, const VectorXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Scalar-design variant with a matrix right-hand side.
MatrixXd apply_projection(const SparseMatrix<double>& inverse_matrix,
        const VectorXd& inverse_times_design, double design_cross_inverse, const MatrixXd& rhs) {
    return inverse_matrix * rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Diagonal-V^{-1} variant (GRM eigenbasis): V^{-1}rhs is an elementwise product.
VectorXd apply_projection(const VectorXd& inverse_diagonal,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const VectorXd& rhs) {
    return inverse_diagonal.cwiseProduct(rhs) -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Diagonal-V^{-1} variant with a matrix right-hand side.
MatrixXd apply_projection(const VectorXd& inverse_diagonal,
        const MatrixXd& inverse_times_design, const MatrixXd& design_cross_inverse, const MatrixXd& rhs) {
    MatrixXd projected_rhs = rhs.array().colwise() * inverse_diagonal.array();
    return projected_rhs -
            inverse_times_design * (design_cross_inverse * (inverse_times_design.transpose() * rhs));
}

// Inverts the singleton (diagonal) part of V for the GxE model. The diagonal of V
// for each unrelated sample is
//   var = grm*σ²_g + gxe*σ²_gxe (+ nxe*σ²_noise) + σ²_resid,
// where the variance components are taken from `varcom` (last entry = residual).
// Appends 1/var as diagonal triplets of V^{-1} and returns sum(log var), the
// block's contribution to log|V|.
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

// Inverts one dense within-group block of V for the GxE model. Forms the block
//   V_g = σ²_g*GRM + σ²_gxe*GxE (+ σ²_noise*diag(nxe)) + σ²_resid*I,
// factorizes it with LDLT, appends the dense inverse as triplets of V^{-1}, and
// returns the block's log-determinant contribution.
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

// Assembles the full inverse covariance V^{-1} for the GxE model under the given
// variance components. Because V is block-diagonal, it is inverted block by block
// (singleton diagonal block + dense within-group blocks) and the log-determinants
// summed. Returns the sparse V^{-1} and log|V|.
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

// ===== Variance-component estimation (REML) helpers =====

// Computes a REML update step that keeps every variance component positive.
// It blends the average-information matrix (ai_mat) with the EM information
// (em_mat) as (1-γ)*AI + γ*EM, increasing γ from 0 in steps of 0.01 until the
// resulting varcom + delta has all positive entries. delta solves
// (blended information) * delta = score (fd_mat). Returns the step and the
// updated variance components.
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

// Writes the estimated variance components and their standard errors, one per
// line, to <out_file>.var.
void write_variance_components_with_se(const std::string& out_file, const VectorXd& varcom,
        const VectorXd& se_varcom, const std::vector<std::string>& component_names) {
    ofstream fout(out_file + ".var");
    if(!fout.is_open()){ fatal_error("Fail to open the output file: {}.var", out_file); }
    fout << "component\tvariance\tse\n";
    for (std::int64_t i = 0; i < varcom.size(); ++i) {
        fout << component_names[i] << "\t" << varcom(i) << "\t" << se_varcom(i) << "\n";
    }
}

// Write per-component heritability / variance proportion and delta-method SE,
// matching the rows of the .var file (G, GxE, [NxE], residual) plus h2_total.
//   prop_i  = varcom_i * d_i / V_p,   V_p = sum_k varcom_k d_k
//   d = per-kernel diagonal means [d_g, d_gxe, (d_nxe), 1]
//   SE via the delta method on the full AI^{-1}.
// For G and GxE the proportion IS the heritability; for NxE/residual it is the
// phenotypic-variance proportion (not counted as genetic h2).
// Writes per-component variance proportion (h2) and delta-method SE to <out>.h2.
//   prop_i = varcom_i * d_i / Vp,   Vp = sum_k varcom_k * d_k
//   d      = per-kernel diagonal means (e.g. [tr(K)/n, ..., 1] with 1 for residual)
//   SE via the delta method on the full AI^{-1}.
// component_names labels each row (length == varcom.size()); genetic_indices lists
// the components summed into the final "total" row (e.g. {0} = G for the main model,
// {0,1} = G+GxE for the GxE model).
void write_heritability_with_se(const std::string& out_file, const VectorXd& varcom,
        const MatrixXd& ai_mat_inv, const VectorXd& d,
        const std::vector<std::string>& component_names,
        const std::vector<int>& genetic_indices) {
    const int nc = static_cast<int>(varcom.size());
    const double Vp = varcom.dot(d);

    auto ratio_se = [&](const VectorXd& c) -> double {
        const double num = c.dot(varcom);
        VectorXd grad(nc);
        for (int k = 0; k < nc; ++k) grad(k) = (c(k) * Vp - num * d(k)) / (Vp * Vp);
        const double var = grad.dot(ai_mat_inv * grad);
        return var > 0.0 ? std::sqrt(var) : 0.0;
    };
    auto comp_se = [&](int k) -> double {
        VectorXd c = VectorXd::Zero(nc);
        c(k) = d(k);
        return ratio_se(c);
    };

    const VectorXd prop = (varcom.cwiseProduct(d)) / Vp;

    ofstream fout(out_file + ".h2");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.h2", out_file);
    }
    fout << "component\th2\tse\n";
    for (int k = 0; k < nc; ++k) {
        fout << component_names[k] << "\t" << prop(k) << "\t" << comp_se(k) << "\n";
    }
    // "total" = combined proportion of the genetic / interaction components.
    VectorXd c_tot = VectorXd::Zero(nc);
    for (int idx : genetic_indices) c_tot(idx) = d(idx);
    fout << "total\t" << (c_tot.dot(varcom) / Vp) << "\t" << ratio_se(c_tot) << "\n";
}

// Writes the estimated variance components (no standard errors) to <out_file>.var.
void write_variance_components(const std::string& out_file, const VectorXd& varcom) {
    ofstream fout(out_file + ".var");
    if(!fout.is_open()){
        fatal_error("Fail to open the output file: {}.var", out_file);
    }

    fout << varcom << std::endl;
}

// Performs one AI-REML iteration for the 2-component main-effect model, working
// entirely in the GRM eigenbasis where everything reduces to diagonal algebra.
// Inputs are the GRM eigenvalues and the rotated design / response (xmat_trans,
// y_trans). Returns the -2 log-likelihood, the score vector and the AI matrix for
// the two components (genetic σ²_g, residual σ²_e).
MainModelIterationState evaluate_main_model_iteration(const VectorXd& grm_eigenvals,
        const MatrixXd& xmat_trans, const MatrixXd& y_trans, const VectorXd& varcom) {
    MainModelIterationState iteration_state;

    // In the eigenbasis V is diagonal with entries λ_i σ²_g + σ²_e, so V^{-1} is
    // just the reciprocal; -log|V^{-1}| = log|V| contributes to -2 logL.
    const VectorXd inverse_diagonal =
        (1.0 / (grm_eigenvals.array() * varcom(0) + varcom(1))).matrix();
    iteration_state.logL = -(inverse_diagonal.array().log()).sum();

    // Add the fixed-effect determinant term log|X'V^{-1}X| (REML correction).
    ProjectionCache fixed_effect_projection =
        build_projection_cache(inverse_diagonal, xmat_trans, true);
    iteration_state.logL += fixed_effect_projection.design_logdet;

    // Py = projected response; y'Py completes the -2 logL.
    const VectorXd py = apply_projection(
        inverse_diagonal,
        fixed_effect_projection.inverse_times_design,
        fixed_effect_projection.design_cross_inverse,
        y_trans).col(0);
    iteration_state.logL += (y_trans.transpose() * py).sum();

    // dpy = D P y (D = diag(eigenvalues)) and the further projections PDPy, PPy
    // are the building blocks of the REML score and AI matrix.
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

    // Score for σ²_g: -0.5 [ tr(PD) - y'PDPy ], with tr(PD) computed as
    // tr(V^{-1}D) minus the fixed-effect correction term.
    iteration_state.fd_mat(0) = (inverse_diagonal.cwiseProduct(grm_eigenvals)).sum();
    MatrixXd inverse_times_genetic_design =
        fixed_effect_projection.inverse_times_design.array().colwise() * grm_eigenvals.array();
    iteration_state.fd_mat(0) -= (
        fixed_effect_projection.design_cross_inverse.cwiseProduct(
            inverse_times_genetic_design.transpose() * fixed_effect_projection.inverse_times_design)).sum();
    iteration_state.fd_mat(0) -= (py.transpose() * dpy).sum();
    iteration_state.fd_mat(0) *= -0.5;

    // Score for σ²_e: -0.5 [ tr(P) - y'PPy ].
    iteration_state.fd_mat(1) = inverse_diagonal.sum() - (
        fixed_effect_projection.design_cross_inverse.cwiseProduct(
            fixed_effect_projection.inverse_times_design.transpose() *
            fixed_effect_projection.inverse_times_design)).sum();
    iteration_state.fd_mat(1) -= (py.transpose() * py).sum();
    iteration_state.fd_mat(1) *= -0.5;

    // Average-information matrix: 0.5 * (K_i P y)' P (K_j P y) for each pair.
    iteration_state.ai_mat(0, 0) = 0.5 * dpy.dot(pdpy);
    iteration_state.ai_mat(1, 1) = 0.5 * py.dot(ppy);
    iteration_state.ai_mat(0, 1) = 0.5 * dpy.dot(ppy);
    iteration_state.ai_mat(1, 0) = iteration_state.ai_mat(0, 1);

    return iteration_state;
}

// ===== Per-SNP scanning helpers (.bed streaming, random-SNP calibration) =====

// Prepares metadata for scanning the .bed file: maps the analysis sample IDs to
// their .fam row indices, records the SNP/sample counts, and validates that the
// .bed file size matches those counts.
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

// Samples `num_random_snp` distinct SNP indices in [0, num_snp) for the gamma
// calibration panel.
std::vector<std::int64_t> sample_random_snp_indices(std::int64_t num_snp, std::int64_t num_random_snp) {
    return sample_unique_integers(0, num_snp, num_random_snp);
}

// Reads the random-SNP genotype panel from the .bed file and mean-centers each
// SNP column (subtracting 2*allele-frequency). Also returns per-SNP frequency,
// missing rate and non-missing genotype counts.
void load_centered_random_snp_panel(GENO& geno, const string& bed_file,
        const std::vector<std::int64_t>& snp_indices, const std::vector<std::int64_t>& id_index_in_bed_vec,
        MatrixXd& snp_mat, VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr) {
    geno.read_bed_by_snp_indices(
        bed_file + ".bed", snp_mat, freq_arr, missing_rate_arr, nobs_geno_arr, snp_indices, id_index_in_bed_vec);
    snp_mat.rowwise() -= 2 * freq_arr.transpose();
}

// Splits the SNP range [start_pos, end_pos) into npart_snp roughly equal chunks
// and returns the start index / count for the given partition_index. The last
// partition absorbs any remainder so the whole range is covered.
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

// Aborts (with a formatted error containing num_used) if the fraction of random
// SNPs that passed QC is below min_fraction, since gamma calibration would then be
// unreliable.
void require_random_snp_fraction(std::int64_t num_used, std::int64_t num_random_snp,
        double min_fraction, const char* error_template) {
    if (num_used * 1.0 / num_random_snp < min_fraction) {
        fatal_error("{}", fmt::format(fmt::runtime(error_template), num_used));
    }
}

// Opens an output file stream, aborting via fatal_error if it cannot be created.
ofstream open_output_stream(const string& file_path) {
    ofstream fout(file_path);
    if(!fout.is_open()){
        fatal_error("Fail to open output file: {}", file_path);
    }
    return fout;
}

// ===== CLI definition, validation and logging =====

// Declares all command-line flags/options on the CLI11 app and binds them to the
// fields of `options`, including help text and defaults. Establishes that
// --test-main and --test-gxe are mutually exclusive.
void configure_run_cli(CLI::App& app, RunOptions& options) {
    app.description(R"(
    Quick Start:

      Test GxE:
        fastgxe --test-gxe --grm grm_file --bfile bed_file --data data_file --trait BMI --env-int age smoking:alcohol --out test_gxe

      Test SNP Main Effects:
        fastgxe --test-main --grm grm_file --bfile bed_file --data data_file --trait BMI --out test_main
    )");

    app.add_option("-p,--threads", options.threads,
        "Number of threads to use (default: 10).")
        ->default_val(10);

    auto* test_main_flag = app.add_flag("--test-main", options.test_main,
        "Enable testing of SNP main effects for continuous traits (mutually exclusive with --test-gxe).");
    auto* test_gxe_flag = app.add_flag("--test-gxe", options.test_gxe,
        "Enable testing of SNP-environment interactions (mutually exclusive with --test-main).");
    test_main_flag->excludes(test_gxe_flag);
    test_gxe_flag->excludes(test_main_flag);

    app.add_flag("--no-noisebye", [&options](int count) {
        if (count > 0) options.no_noisebye = true;
    }, "Disable noise-by-environment interaction terms.\n"
       "  - Noise-by-environment interactions are ENABLED by default.\n"
       "  - Use --no-noisebye to turn it OFF.");

    app.add_flag("--keep-random", options.keep_random,
        "Also write the per-SNP random-SNP gamma-calibration files (*.random.res). Off by default; the calibration summary is always written to <out>.gamma.");

    app.add_option("--data", options.data_file, "Path to input data file (required).")->required();
    app.add_option("--trait", options.trait_vec,
        "Trait(s) to analyze.\n"
        "  - Use 1 trait for --test-main and --test-gxe.")
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

    app.add_flag("--log-transform", options.log_transform,
        "Natural-log transform the raw phenotype before analysis (requires positive values).\n"
        "  - DISABLED by default.");
    app.add_flag("--inv-normal", options.inv_normal,
        "Rank-based inverse-normal transform of the phenotype (Blom, mid-ranks for ties).\n"
        "  - DISABLED by default.\n"
        "  - For --test-main: applied to the (optionally log-transformed) phenotype.\n"
        "  - For --test-gxe: applied to the residual AFTER covariate/environment correction.");

    app.add_option("--missing-data", options.missing_in_data_vec,
        "List of missing value indicators for phenotype/covariates.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-data . -999 \"?\").")
        ->expected(-1);

    app.add_option("--maxiter", options.maxiter, "Maximum number of optimization iterations (default: 100).")
        ->default_val(100);
    app.add_option("--cc-par", options.cc_par, "Convergence threshold for parameter updates (default: 1e-7).")
        ->default_val(1e-7);
    app.add_option("--cc-gra", options.cc_gra, "Convergence threshold for gradient norm (default: 1e-6).")
        ->default_val(1e-6);
    app.add_option("--cc-logL", options.cc_logL, "Convergence threshold for log-likelihood change (default: 5e-5).")
        ->default_val(5e-5);
    app.add_option("--snp-range", options.snp_range_vec, "Specify start and end SNPs (expects 2 values).")->expected(2);
    app.add_option("--npart-SNP", options.npart_snp, "Number of SNP partitions (default: 1000).")->default_val(1000);
    app.add_option("--exact-cut", options.exact_cut,
        "Adaptive exact-refit threshold for --test-gxe (default: 1e-5).\n"
        "  - A promising SNP whose fast approximate p-value is below this value is\n"
        "    recomputed exactly (x'V^-1 x), so top hits get exact p-values while the\n"
        "    genome-wide scan stays fast.\n"
        "  - Range [0, 1]. Set to 0 to disable (pure fast approximation);\n"
        "    set to 1 to force exact refit for ALL promising SNPs (slowest, most accurate).")
        ->default_val(1e-5);
    app.add_option("--p-approx-cut", options.p_approx_cut, "P-value approximation threshold (default: 1e-3).")->default_val(1e-3);
    app.add_option("--maf", options.maf_cut, "Minor allele frequency (MAF) cutoff (default: 0.01).")->default_val(0.01);
    app.add_option("--missing-rate", options.missing_rate_cut, "Missing genotype rate cutoff (default: 0.05).")->default_val(0.05);
    app.add_option("--num-random-snp", options.num_random_snp, "Number of randomly selected SNPs (default: 1000).")->default_val(1000);
    app.add_option("--split-task", options.split_task_vec,
        "Partition the task for parallel execution (expects 2 values: total parts, current part).")
        ->expected(2);
}

// Validates parsed options: exactly one analysis mode, --env-int required for
// GxE, exactly one trait, and a sane thread count. Aborts on violation.
void validate_run_options(RunOptions& options) {
    const int num_analysis_modes =
        static_cast<int>(options.test_main) +
        static_cast<int>(options.test_gxe);
    if (num_analysis_modes != 1) {
        fatal_error("Exactly one of --test-main or --test-gxe must be specified.");
    }
    if (options.test_gxe && options.bye_vec.empty()) {
        fatal_error("--env-int is required when --test-gxe is specified.");
    }
    if (options.trait_vec.size() != 1) {
        fatal_error("This analysis mode requires exactly one trait.");
    }
    if (options.exact_cut < 0.0 || options.exact_cut > 1.0) {
        fatal_error("--exact-cut ({}) must be between 0 and 1.", options.exact_cut);
    }
    if (options.threads <= 0) {
        options.threads = 10;
    }
}

// Logs all resolved run options through spdlog for reproducibility/debugging.
void log_run_options(const RunOptions& options) {
    spdlog::info("=== Parsed Arguments ===");
    spdlog::info("Threads: {}", options.threads);
    if (options.test_main) {
        spdlog::info("Analysis Mode: Continuous SNP main effects");
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
    spdlog::info("Phenotype log-transform: {}", options.log_transform ? "ENABLED" : "DISABLED");
    spdlog::info("Phenotype inverse-normal transform: {}", options.inv_normal ? "ENABLED" : "DISABLED");
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
    spdlog::info("GxE exact-refit cutoff: {}", options.exact_cut <= 0.0 ? "disabled (pure approximate)"
            : (options.exact_cut >= 1.0 ? "all promising SNPs (exact)"
            : "p < " + std::to_string(options.exact_cut)));
    spdlog::info("Minor Allele Frequency Cutoff: {}", options.maf_cut);
    spdlog::info("Missing Rate Cutoff: {}", options.missing_rate_cut);
    spdlog::info("Number of Random SNPs: {}", options.num_random_snp);
    if(!options.split_task_vec.empty()) {
        spdlog::info("Task Partitioning: {} parts, running part {}", options.split_task_vec[0], options.split_task_vec[1]);
    }
    spdlog::info("========================");
}

// Resolves the user-supplied variable names to column indices in the data-file
// header. Expands range specifications (a:b) for covariates/classes/environments,
// verifies all required names exist, and enforces that interacting environments
// are not also listed as covariates or classes. Returns the resolved indices.
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

// ===== fastGxE: data loading and sample/GRM alignment =====

// Extracts the (column 0) sample IDs retained in the phenotype table after
// filtering, asserting they are unique and non-empty. Returns the ID list.
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

// Reads <agrm_file>.grm.group, keeping only records whose sample ID is in the
// retained set, recomputes each group's size, and sorts records by (group size,
// group id, GRM index). The resulting order places unrelated singletons first and
// groups related samples contiguously, which is what makes the reordered GRM
// block-diagonal. Returns the sorted group records.
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

    spdlog::info("Recomputing relatedness-group sizes");
    for(auto& record : group_records){
        record.group_size = group_size_map[record.group_id];
    }

    spdlog::info("Sorting samples by relatedness group (size, then group id)");
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

// Translates the sorted group records into the corresponding sample-ID order by
// looking up each record's (1-based) GRM index in the GRM ID list. This is the
// target order to which the phenotype table is then reordered.
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

// Logs which fixed-effect design columns were dropped as linearly dependent
// (no-op when none were removed).
void fastGxE::log_removed_design_columns(const vector<std::int64_t>& removed_column_indices) const {
    if (removed_column_indices.empty()) {
        return;
    }

    spdlog::info("Removing linearly dependent fixed-effect design columns:");
    std::string removed_columns;
    for (const auto& col : removed_column_indices) {
        removed_columns += std::to_string(col) + " ";
    }
    spdlog::info("{}", removed_columns);
}

// Builds the fixed-effect design matrix m_xmat (and response m_y) from the
// covariate/class columns, removes any linearly dependent design columns, and
// loads the interacting-environment columns into m_bye_mat.
void fastGxE::build_fixed_and_environment_matrices(PHEN& phenoA,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait) {
    spdlog::info("Building fixed-effect design matrix");
    vector<std::int64_t> index_fixed_effect_vec;
    vector<string> fixed_effect_name_vec;
    m_xmat = phenoA.build_fixed_effect_design_matrix(
        covariate_arr, class_arr, index_fixed_effect_vec, fixed_effect_name_vec, trait, m_y);

    vector<std::int64_t> removed_column_indices;
    phenoA.remove_dependent_design_columns(
        m_xmat, index_fixed_effect_vec, fixed_effect_name_vec, removed_column_indices);
    log_removed_design_columns(removed_column_indices);

    spdlog::info("Loading interacting-environment covariates");
    if(!bye_arr.empty()){
        phenoA.get_columns_as_matrix(bye_arr, m_bye_mat);
    }
}

// Builds the block-diagonal GRM layout from the sorted group records. Singletons
// (group size 1) all go into block 0 (stored as an N x 1 column); each multi-member
// relatedness group becomes its own dense block. Also fills sample_location_map,
// mapping each sample's GRM index to its (block_id, row_index) so the kinship
// values can later be placed into the correct block. Returns the GroupLayout.
fastGxE::GroupLayout fastGxE::build_group_layout(const vector<GroupRecord>& group_records) const {
    spdlog::info("Building block-diagonal GRM layout from relatedness groups");

    std::vector<std::int64_t> singleton_indices;
    std::vector<vector<std::int64_t>> dense_group_indices;
    singleton_indices.reserve(group_records.size());
    dense_group_indices.reserve(group_records.size());

    // Partition samples: singletons accumulate into one list; each contiguous run
    // of a shared group_id (records are already sorted) becomes a dense group.
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
    // Block 0 holds all singletons as a single N x 1 column.
    layout.grm_blocks.emplace_back(MatrixXd::Zero(singleton_indices.size(), 1));

    // Singletons are tagged with a negative block id (encoding their row only);
    // loader uses row_index to place each diagonal kinship value.
    for(std::int64_t i = 0; i < static_cast<std::int64_t>(singleton_indices.size()); ++i){
        layout.sample_location_map.emplace(singleton_indices[i], SampleBlockLocation{-(i + 1), i});
    }

    // Each dense group becomes its own square block (block_id >= 1).
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

// Streams the sparse GRM file (<agrm_file>.grm.sp.bin) of (index0, index1, value)
// records and scatters each value into the correct dense/singleton block. Pairs
// are skipped unless both samples are retained and fall in the same block (cross-
// block kinship is treated as zero); dense blocks are filled symmetrically and the
// singleton block only stores its diagonal.
void fastGxE::load_grouped_grm_blocks(const string& agrm_file,
        const std::unordered_map<std::int64_t, SampleBlockLocation>& sample_location_map,
        vector<MatrixXd>& grm_blocks) const {
    spdlog::info("Reading GRM blocks from {}.grm.sp.bin", agrm_file);

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

// Default constructor.
fastGxE::fastGxE(){

}

// Destructor.
fastGxE::~fastGxE() {
}

// Records the output prefix and the counts of traits / covariates / classes /
// interacting environments on the object, logging them for the user.
void fastGxE::initialize_analysis_metadata(const string& out_file,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait) {
    spdlog::info("Preparing data for analysis...");
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

// Reads the data file and keeps only rows whose sample ID is present in the GRM
// and whose used columns have no missing values. Populates m_id_in_data_vec /
// m_num_id with the surviving sample order and count.
void fastGxE::load_filtered_phenotype(PHEN& phenoA, const string& data_file,
        const vector<std::int64_t>& covariate_arr, const vector<std::int64_t>& class_arr,
        const vector<std::int64_t>& bye_arr, const vector<std::int64_t>& trait,
        const vector<string>& grm_sample_ids, const vector<string>& missing_in_data_vec) {
    spdlog::info("Reading data file: {}", data_file);
    const vector<std::int64_t> col_used_vec = merge_unique_vectors({covariate_arr, class_arr, bye_arr, trait});
    phenoA.read_and_filter(data_file, {0}, {grm_sample_ids}, col_used_vec, missing_in_data_vec);

    m_id_in_data_vec = extract_retained_sample_ids(phenoA);
    m_num_id = m_id_in_data_vec.size();
    spdlog::info("The number of used iids in data file: {}", m_num_id);
}

// Reorders the phenotype table so its sample order matches the relatedness-group
// ordering used to make the GRM block-diagonal. Loads/sorts the group records,
// reorders the PHEN rows accordingly, refreshes m_id_in_data_vec, and returns the
// sorted group records for later block construction.
vector<fastGxE::GroupRecord> fastGxE::reorder_phenotype_by_group(PHEN& phenoA, const string& agrm_file,
        const vector<string>& grm_sample_ids) {
    const vector<GroupRecord> group_records =
        load_and_sort_group_records(agrm_file, m_id_in_data_vec, m_num_id);

    spdlog::info("Reordering phenotype rows to match GRM group order");
    phenoA.reorder_records_by_sample_id_order(build_reordered_sample_ids(group_records, grm_sample_ids));
    m_id_in_data_vec = extract_retained_sample_ids(phenoA);
    return group_records;
}

// Top-level data-preparation step. Loads the GRM sample IDs, filters the
// phenotype data to those samples, reorders individuals by relatedness group so
// the GRM is block-diagonal, builds the fixed-effect / environment design
// matrices, and assembles the per-group GRM blocks (m_grm_mat_group_vec). All
// later variance-component and SNP-scan routines consume these cached structures.
// Returns 0 on success.
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
    spdlog::info("Reading sample IDs from GRM: {}.grm", agrm_file);
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

    spdlog::info("Data preparation completed: {} samples retained", m_num_id);
    return 0;
}








// ===== GRM assembly / eigen-decomposition =====

// Finalizes the GRM representation from the per-group blocks.
// - If use_eigen is false (GxE mode): assembles the block-diagonal sparse GRM
//   m_grm_mat directly.
// - If use_eigen is true (main-effect mode): block-wise eigen-decomposes the GRM,
//   producing the eigenvalues m_grm_eigenvals and the sparse rotation
//   m_grm_eigenvecs (Q), then pre-rotates the response and design into the
//   eigenbasis (m_y_trans, m_xmat_trans) so the null model becomes diagonal.
/**
 * @Description: Sparse GRM or perform Eigen decomposition of GRM
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
        spdlog::info("Eigen-decomposing the GRM (block-wise)");
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



// ===== GxE-specific preprocessing =====

// Prepares the environment design for the GxE model. Checks the environment
// matrix is full column rank, optionally standardizes it, and appends it to the
// fixed-effect design (checking the combined design stays full rank). When
// phen_correct is true it regresses the fixed effects + environments out of the
// response and collapses the fixed-effect design to a single intercept column, so
// downstream GxE fitting works on the residualized phenotype.
/**
 * @Description: Scale the interaction environment covariates and add to xmat
 */
void fastGxE::pre_data_GxE(bool standardize_env, bool phen_correct, bool inv_normal){
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

        // For GxE the inverse-normal transform is applied to the corrected
        // residual (not the raw trait), so spurious interactions driven by
        // phenotype non-normality are removed after covariate/environment adjustment.
        if(inv_normal){
            spdlog::info("Applying inverse-normal transform to the corrected phenotype residual");
            apply_inverse_normal_transform(m_y);
        }
    }
}


// ===== Variance-component estimation (REML) =====

// Fits the null-model variance components for the GxE (improved structLMM) model
// by AI-REML. The covariance is
//   V = σ²_g*GRM + σ²_gxe*GxE + [σ²_noise*diag(nxe) +] σ²_resid*I,
// giving 4 components (or 3 when no_noisebye is set). Each iteration rebuilds the
// block-diagonal V^{-1} and log|V|, evaluates the REML -2 log-likelihood, the
// score vector fd_mat and the average-information matrix ai_mat (entries
// 0.5*(K_i P y)' P (K_j P y)), then takes a positive-variance step blending AI and
// EM information. On convergence it caches m_Vi0 = V^{-1}, m_V0_logdet, the null
// -2 logL (m_logL_null) and the components (m_varcom_null), writes them with SEs,
// and returns the estimates. The parameters maxiter0/cc_par0/cc_gra0/cc_logL0 are
// the iteration cap and convergence thresholds (init_varcom is the starting point).
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
    MatrixXd ai_mat_inv = MatrixXd::Identity(num_cov, num_cov);  // full AI^-1 for h2 SE

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
        // Score for the genetic component σ²_g (K = GRM):
        //   fd = tr(P K) - y' P K P y, with tr(PK) = tr(V^{-1}K) minus the
        //   fixed-effect correction tr((X'V^{-1}X)^{-1} X'V^{-1} K V^{-1} X).
        MatrixXd KPy = m_grm_mat * projected_response;
        double tr_ViK = (vmat.cwiseProduct(m_grm_mat)).sum();
        double tr_design_cross_genetic = (
            fixed_effect_projection.design_cross_inverse.cwiseProduct(design_cross_genetic_design)).sum();
        fd_mat(0) = tr_ViK - tr_design_cross_genetic - (projected_response.transpose() * KPy).sum();
        KPy_vec[0] = KPy;

        // Score for the GxE-relationship component σ²_gxe (K = GxE matrix).
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

        // Score for the noise-by-environment component σ²_noise (K = diag(nxe)),
        // only present when noise-by-environment is enabled.
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



        // Score for the residual component σ²_resid (K = I): tr(P) - y'PPy.
        fd_mat(num_cov - 1) = vmat.diagonal().sum() - (
            fixed_effect_projection.design_cross_inverse.cwiseProduct(
                fixed_effect_projection.inverse_times_design.transpose() *
                fixed_effect_projection.inverse_times_design)).sum() -
            (projected_response.transpose() * projected_response).sum();
        fd_mat *= -0.5;
        KPy_vec[num_cov - 1] = projected_response;

        // Project each K_i P y once; these P K_i P y vectors feed the AI matrix.
        vector <MatrixXd> PKPy_vec(num_cov);
        for(std::int64_t m = 0; m < num_cov; m++){
            const MatrixXd& KPy_tmp = KPy_vec[m];
            PKPy_vec[m] = apply_projection(
                vmat,
                fixed_effect_projection.inverse_times_design,
                fixed_effect_projection.design_cross_inverse,
                KPy_tmp);
        }

        // AI: average-information matrix, ai(m,n) = 0.5 (K_m P y)' P (K_n P y).
        MatrixXd ai_mat(num_cov, num_cov);
        for(std::int64_t m = 0; m < num_cov; m++){
            ai_mat(m, m) = (KPy_vec[m].transpose() * PKPy_vec[m]).sum();
            for(std::int64_t n = 0; n < m; n++){
                ai_mat(m, n) = ai_mat(n, m) = (KPy_vec[m].transpose() * PKPy_vec[n]).sum();
            }
        }
        ai_mat *= 0.5;
        ai_mat_inv = ai_mat.inverse();
        ai_mat_inv_diag = ai_mat_inv.diagonal();

        // EM information diagonal: n/(2σ⁴) per variance parameter; a valid lower bound when AI is ill-conditioned.
        MatrixXd em_mat = MatrixXd::Zero(num_cov, num_cov);
        em_mat.diagonal() = m_num_id / (2 * varcom.array() * varcom.array());
        
        // Update: AI/EM-blended Newton step kept inside the positive orthant.
        const VarianceUpdateResult update_result = find_positive_variance_update(ai_mat, em_mat, fd_mat, varcom);
        spdlog::info("Updated variances: {}", format_vector_for_log(update_result.updated_varcom));

        varcom = update_result.updated_varcom;
        
        // Convergence: relative parameter change and gradient norm.
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
    
    // Rebuild V^{-1} at the converged components and cache it (m_Vi0) plus log|V|
    // for reuse by the per-SNP GxE score tests.
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
    const std::vector<std::string> comp_names = no_noisebye
        ? std::vector<std::string>{"G", "GxE", "residual"}
        : std::vector<std::string>{"G", "GxE", "NxE", "residual"};
    write_variance_components_with_se(m_out_file, varcom, se_varcom, comp_names);

    // Per-kernel diagonal means for the heritability decomposition:
    //   d_g = tr(K)/n, d_gxe = tr(K_gxe)/n, d_nxe = mean(||W_i||^2/q), d_e = 1.
    VectorXd d_means(num_cov);
    d_means(0) = m_grm_mat.diagonal().sum() / static_cast<double>(m_num_id);
    d_means(1) = gxe_relationship.sparse_matrix.diagonal().sum() / static_cast<double>(m_num_id);
    if(!no_noisebye){
        d_means(2) = nxe_vec.sum() / static_cast<double>(m_num_id);
        d_means(3) = 1.0;
    } else {
        d_means(2) = 1.0;
    }
    write_heritability_with_se(m_out_file, varcom, ai_mat_inv, d_means, comp_names, {0, 1});

    return m_varcom_null;
}




// Fits the 2-component null model for the SNP main-effect scan by AI-REML in the
// GRM eigenbasis, where V is diagonal (λ_i σ²_g + σ²_e). Each iteration calls
// evaluate_main_model_iteration for the -2 logL, score and AI matrix, then takes
// an AI/EM-blended positive-variance step. On convergence it caches the estimates
// in m_varcom_null, writes them out, and returns them. init is reset to ones (the
// model always starts from a unit initialization); maxiter0 and the cc_* values
// are the iteration cap and convergence thresholds.
/**
 * @Description: Estimate variances for model with only one random effect
 * @param {VectorXd&} init
 * @param {int} maxiter0
 * @param {double} cc_par0
 * @param {double} cc_gra0
 * @param {double} cc_logL0
 */
VectorXd fastGxE::varcom_main(VectorXd& init, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0){

    spdlog::info("Estimating variance components via GRM eigendecomposition...");

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

        // EM information diagonal n/(2σ⁴), used to keep the blended step positive.
        MatrixXd em_mat = MatrixXd::Zero(2, 2);
        em_mat(0, 0) = m_num_id / (2 * varcom(0) * varcom(0));
        em_mat(1, 1) = m_num_id / (2 * varcom(1) * varcom(1));

        // AI/EM-blended Newton step (stays in the positive-variance region).
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

    // Recompute the AI matrix at the converged estimates so its inverse gives the
    // sampling covariance of the variance components (for the SEs and the .h2 SE).
    const MainModelIterationState final_state = evaluate_main_model_iteration(
        m_grm_eigenvals, m_xmat_trans, m_y_trans, varcom);
    const MatrixXd ai_mat_inv = final_state.ai_mat.inverse();
    const VectorXd se_varcom = ai_mat_inv.diagonal().cwiseSqrt();
    const std::vector<std::string> comp_names = std::vector<std::string>{"G", "residual"};
    write_variance_components_with_se(m_out_file, varcom, se_varcom, comp_names);

    // Heritability: G = sigma2_g * d_g / Vp (SNP h2), residual = sigma2_e / Vp;
    // total = G. d_g = tr(GRM)/n = mean of the GRM eigenvalues in the eigenbasis.
    VectorXd d_means(2);
    d_means(0) = m_grm_eigenvals.sum() / static_cast<double>(m_num_id);
    d_means(1) = 1.0;
    write_heritability_with_se(m_out_file, varcom, ai_mat_inv, d_means,
            comp_names, {0});

    return m_varcom_null;
}

// ===== Per-SNP association scans =====

// Genome-wide SNP main-effect scan under the 2-component null model. Working in
// the GRM eigenbasis it (1) caches the null fixed-effect projection, (2) uses a
// random panel of SNPs to estimate a single calibration factor gamma =
// x'V^{-1}x / x'x (the score-test denominator approximation), then (3) streams the
// genotype matrix in npart_snp chunks and, per SNP, forms a fast score statistic
// (beta = se^2 * x' (P y), se^2 = 1/(gamma * x'x)) and its chi-square (1 df)
// p-value. Results pass QC (MAF / missing-rate) and are written to <out>.res,
// with the random panel written to <out>.random.res.
/**
 * @Description: Test SNPs with one random-effect model
 */
void fastGxE::test_main(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                bool keep_random,
                std::int64_t num_random_snp,
                double maf_cut, double missing_rate_cut){
    spdlog::info("Association: scanning SNPs for main effects");
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
    // Rotate the random SNPs into the GRM eigenbasis so projections are diagonal.
    snp_mat_part_random = m_grm_eigenvecs.transpose() * snp_mat_part_random;

    // Per random SNP: exact score-test p-value and the calibration ratio
    // gamma_k = x' P x / x'x, later averaged into the single gamma factor.
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
    ofstream fout;
    if(keep_random){
        fout = open_output_stream(m_out_file + ".random.assoc");
        fout << "order\taf\tN\tgamma\tp" << std::endl;
    }
    // Average gamma over QC-passing random SNPs to get one calibration factor.
    for(std::int64_t k = 0; k < num_random_snp; k++){
        if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
            if(keep_random){
                fout << random_snp_index_vec[k] << "\t" << freq_arr(k) << "\t" <<
                        static_cast<long long>(nobs_geno_arr(k)) << "\t" << gamma_correct_Vec(k) << "\t" << p_Vec(k) << std::endl;
            }
            num_random_snp_used++;
            gamma_correct += gamma_correct_Vec[k];
        }
    }
    if(keep_random){ fout.close(); }

    require_random_snp_fraction(
        num_random_snp_used, num_random_snp, 0.1,
        "The number of used random SNPs to calculate gamma is: {}, which is less than 10% of random SNPs");
    gamma_correct /= num_random_snp_used;
    spdlog::info("Gamma factor: {}", gamma_correct);
    spdlog::info("The number of used random SNPs is: {}", num_random_snp_used);

    // Always write a compact calibration summary to <out>.gamma.
    {
        ofstream fout_gamma = open_output_stream(m_out_file + ".gamma");
        fout_gamma << "term\tvalue" << std::endl;
        fout_gamma << "gamma_main\t" << gamma_correct << std::endl;
        fout_gamma << "n_random_used\t" << num_random_snp_used << std::endl;
        fout_gamma.close();
    }

    // Move back to the sample space for fast genome-wide SNP scanning once the
    // random-SNP calibration factor has been estimated.
    spdlog::info("Scanning SNPs for main effects; writing results to {}.assoc", m_out_file);
    fout = open_output_stream(m_out_file + ".assoc");
    fout << "order\tchrom\tSNP\tbase\tallele1\tallele2\taf\tN\tbeta\tse\tp" << std::endl;

    // P y mapped back to sample space; dotting each centered SNP with it gives
    // the score numerator x'(Py) in O(n) per SNP.
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

        // Test: per SNP compute x'x, the calibrated variance se^2 = 1/(gamma*x'x),
        // the effect beta = se^2 * x'(Py), and the chi-square(1) p-value.
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
        fastGxE::output(fout, snp_info_vec, freq_arr, nobs_geno_arr,
            bufs.eff.data(), bufs.se.data(), bufs.p, partition.start_snp, partition.num_snp_read);
    }
    fout.close();
}

// ===== Result output =====

// Writes one result row per SNP for the main-effect scan (annotation, frequency,
// missing rate, effect, SE, p-value) to the output stream.
void fastGxE::output(ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& nobs_geno_arr,
                const double* eff_arr, const double* se_arr, const VectorXd& p_arr,
                std::int64_t start_snp, std::int64_t num_snp_read){
    for(std::int64_t i = 0; i < num_snp_read; i++){
        fout << start_snp + i << "\t" << snp_info_vec[i] << "\t" << freq_arr(i) << "\t" << static_cast<long long>(nobs_geno_arr(i))
                << "\t" << eff_arr[i] << "\t" << se_arr[i]
                << "\t" << p_arr(i) << std::endl;
    }
}

// Appends a "<first>_<second>" suffix to the output prefix so that range- or
// task-split runs write to distinct files.
void fastGxE::reset_output_prefix(const string& out_file, std::int64_t first, std::int64_t second){
    m_out_file = out_file + "." + std::to_string(first) + "_"  + std::to_string(second);
}

// SNP QC predicate: passes if the allele frequency is within [maf_cut, 1-maf_cut]
// and the missing-genotype rate is below missing_rate_cut.
bool fastGxE::isValidSnp(double af, double missing_rate, double maf_cut, double missing_rate_cut) {
    return af > maf_cut && af < 1 - maf_cut && missing_rate < missing_rate_cut;
}



// Multi-environment GxE scan that produces a single combined GxE statistic per
// SNP (structLMM-style), as opposed to the per-environment breakdown of test_GxE.
// Using the cached null V^{-1} (m_Vi0): (1) a random SNP panel calibrates the
// main-effect gamma and the matrices gamma_GxE (with and without the SNP main
// effect) whose eigenvalues weight a mixture-of-chi-squares; (2) per genome-wide
// SNP it computes the main-effect test, then a score statistic
//   score = || E' diag(x) P y ||^2
// whose null distribution is a weighted sum of chi-squares evaluated with the
// saddlepoint approximation (saddle). When the approximate main- or GxE-p crosses
// p_approx_cut it refits exactly (exact == true) or with an approximate calibration.
// Results go to <out>.res; the random panel to <out>.random.res.
/**
 * @Description: Association of GxE using an improved structLMM model
 */
void fastGxE::test_GxE_multi(const string& bed_file, std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                bool exact, std::int64_t num_random_snp,
                double p_approx_cut, double maf_cut, double missing_rate_cut){
    spdlog::info("Association: combined multi-environment GxE scan");
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
    // Mean-center each random SNP column.
    snp_mat_part_random.rowwise() -= 2*freq_arr.transpose();

    // Per-random-SNP calibration accumulators: main-effect gamma, and the GxE
    // gamma matrices (with / without the SNP main effect) plus their statistics.
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

            // main: exact LMM main-effect estimate and gamma = (x'V^{-1}x)/(x'x).
            double beta_main_var = 1 / (snpk.transpose() * this->m_Vi0 * snpk);
            se_main_Vec(i) = sqrt(beta_main_var);
            double beta_main = beta_main_var * snpk.transpose() * Vi0y;
            beta_main_Vec(i) = beta_main;
            double p_main = gsl_cdf_chisq_Q(beta_main * beta_main / beta_main_var, 1);
            p_main_Vec(i) = p_main;
            gamma_main_Vec(i) = (1 / beta_main_var) / snp_norm2;

            // GxE without main effect: score = ||(E∘x)' V^{-1} y||^2; its null is
            // a weighted chi-square mixture with weights = eig(EoG' V^{-1} EoG),
            // and saddle() gives the p-value. EoG = E elementwise-scaled by the SNP.
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

            // GxE with main effect: same score test but with the SNP main effect
            // projected out of y (P from the scalar projection on snpk).
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
    
    // Eigenvalues of the averaged gamma matrices become the chi-square mixture
    // weights used (scaled by x'x per SNP) in the genome-wide saddlepoint tests.
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

                // GxE without main effect: fast combined score using the
                // precomputed chi-square weights (scaled by x'x) for the p-value.
                vdMul(num_id_used, Vi0y_pt, snp_mat_part.data() + k*num_id_used, EP0y_part[thread_id].data());
                cblas_dgemv(CblasColMajor, CblasTrans, num_id_used, m_num_bye, 1.0, bye_mat_pt, num_id_used,
                            EP0y_part[thread_id].data(), 1, 0.0, WEP0y_part[thread_id].data(), 1);
                double score = cblas_ddot(m_num_bye, WEP0y_part[thread_id].data(), 1, WEP0y_part[thread_id].data(), 1);
                VectorXd chi_coef_Vec_kth = chi2_coef_GxEnoMain_Vec * snpk_norm2;
                double p_gxe = saddle(score, chi_coef_Vec_kth);

                // Only the promising SNPs are refined with the (more expensive)
                // main-effect-adjusted GxE test: exact (exact == true) or approximate.
                if(p_main < p_approx_cut || p_gxe < p_approx_cut){
                    // GxE with main effect
                    Eigen::Map<Eigen::VectorXd> snpk(snp_mat_part.data() + k * num_id_used, num_id_used);
                    ScalarProjectionCache snp_main_projection = build_scalar_projection_cache(this->m_Vi0, snpk);
                    VectorXd projected_null_response = apply_projection(
                        this->m_Vi0,
                        snp_main_projection.inverse_times_design,
                        snp_main_projection.design_cross_inverse,
                        m_y).col(0);
                    
                    if(exact){
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

// Writes one result row per SNP for the multi-environment GxE scan: annotation,
// frequency, missing rate, the main-effect estimate and the combined GxE score /
// p-value.
void fastGxE::output_GxE_multi(std::ofstream& fout, const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& missing_rate_arr,
            const VectorXd& beta_main_Vec, const VectorXd& se_main_Vec, const VectorXd& p_main_Vec,
            const VectorXd& score_Vec, const VectorXd& p_Vec, std::int64_t start_snp, std::int64_t num_snp_read){
    for(std::int64_t i = 0; i < num_snp_read; i++){
        fout << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) << " "
            << beta_main_Vec(i) << " " << se_main_Vec(i) << " " << p_main_Vec(i) << " "
            << score_Vec(i) << " " << p_Vec(i) << std::endl;
    }
}

// Main per-environment GxE genome-wide scan (the default --test-gxe path).
// Using the cached null V^{-1} (m_Vi0) it runs four stages, all calibrated on a
// random SNP panel so the genome-wide pass is cheap:
//   1) main effect: average gamma_main = (x'V^{-1}x)/(x'x);
//   2) marginal GxE per environment without the SNP main effect: gamma_GxEnoMain;
//   3) joint SNP + single-environment fits: per-environment 2x2 gamma matrices;
//   4) genome-wide scan: for each SNP compute the main-effect and per-environment
//      interaction z-scores from the calibrated denominators, refit exactly /
//      approximately only when a p-value passes p_approx_cut, combine the
//      per-environment GxE p-values with ACAT (p_single), compute a correlation-
//      based combined score test (p_multi) via the saddlepoint method on the
//      environment-correlation eigenvalues, and ACAT-combine the two into p_gxe.
// Results are written to <out>.res with per-stage random panels alongside.
void fastGxE::test_GxE(const string& bed_file,
                std::int64_t start_pos, std::int64_t end_pos, int npart_snp,
                bool keep_random,
                double exact_cut, std::int64_t num_random_snp,
                double p_approx_cut, double maf_cut, double missing_rate_cut){

    spdlog::info("Association: SNP-by-environment (GxE) scan for individual environments");
    spdlog::info("GxE exact-refit cutoff for promising SNPs: {}",
            exact_cut <= 0.0 ? "disabled" : (exact_cut >= 1.0 ? "all" : "p < " + std::to_string(exact_cut)));
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
    spdlog::info("Calibrating gamma: SNP main effect");
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

    ofstream fout;
    if(keep_random){
        fout = open_output_stream(m_out_file + ".main.random.assoc");
        fout << "order\taf\tN\tgamma\tbeta\tse\tp" << std::endl;
        for(std::int64_t i = 0; i < num_random_snp; i++){
            fout << random_snp_index_vec[i] << "\t" << freq_arr(i) << "\t" << static_cast<long long>(nobs_geno_arr(i)) << "\t"
                  << gamma_main_Vec[i] << "\t" << beta_main_Vec[i] << "\t" << se_main_Vec[i] << "\t" << p_main_Vec[i]
                   << std::endl;
        }
        fout.close();
    }
    gamma_main_Vec.resize(0); beta_main_Vec.resize(0); se_main_Vec.resize(0); p_main_Vec.resize(0); beta_main_var_Vec.resize(0);

    // 2) Random-SNP calibration for marginal GxE terms without the SNP main effect.
    spdlog::info("Calibrating gamma: GxE (no main effect)");
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

    if(keep_random){
        fout = open_output_stream(m_out_file + ".GxEnoMain.random.assoc");
        fout << "order\taf\tN";
        for(std::int64_t i = 0; i < m_num_bye; i++){
            fout << "\tgamma_" << m_bye_name_vec[i] << "\tbeta_" << m_bye_name_vec[i] << "\tse_" << m_bye_name_vec[i] << "\tp_" << m_bye_name_vec[i];
        }
        fout << std::endl;

        for(std::int64_t i = 0; i < num_random_snp; i++){
            fout << random_snp_index_vec[i] << "\t" << freq_arr(i) << "\t" << static_cast<long long>(nobs_geno_arr(i));
            for(std::int64_t j = 0; j < m_num_bye; j++){
                fout <<  "\t" << gamma_GxEnoMain_Mat(i, j) << "\t" << beta_GxEnoMain_Mat(i, j) << "\t"
                << se_GxEnoMain_Mat(i, j) << "\t" << p_GxEnoMain_Mat(i, j);
            }
            fout << std::endl;
        }
        fout.close();
    }
    gamma_GxEnoMain_Mat.resize(0, 0); beta_GxEnoMain_Mat.resize(0, 0); se_GxEnoMain_Mat.resize(0, 0); p_GxEnoMain_Mat.resize(0, 0);

    // 3) Random-SNP calibration for joint SNP + single-environment fits.
    spdlog::info("Calibrating gamma: GxE (with main effect)");
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

    
    if(keep_random){
        fout = open_output_stream(m_out_file + ".GxE.random.assoc");
        fout << "order\taf\tN";
        for(std::int64_t i = 0; i < m_num_bye; i++){
            fout << "\tbeta_" << m_bye_name_vec[i] << "\tse_" << m_bye_name_vec[i] << "\tp_" << m_bye_name_vec[i];
        }
        fout << std::endl;
        for(std::int64_t i = 0; i < num_random_snp; i++){
            fout << random_snp_index_vec[i] << "\t" << freq_arr(i) << "\t" << static_cast<long long>(nobs_geno_arr(i));
            for(std::int64_t j = 0; j < m_num_bye; j++){
                fout << "\t" << beta_Mat(i, j) << "\t" << se_Mat(i, j) << "\t" << p_Mat(i, j);
            }
            fout << std::endl;
        }
        fout.close();
    }

    // Always write a compact calibration summary to <out>.gamma.
    {
        ofstream fout_gamma = open_output_stream(m_out_file + ".gamma");
        fout_gamma << "term\tvalue" << std::endl;
        fout_gamma << "gamma_main\t" << gamma_main << std::endl;
        for(std::int64_t j = 0; j < m_num_bye; j++){
            fout_gamma << "gamma_noMain_" << m_bye_name_vec[j] << "\t" << gamma_GxEnoMain_Vec(j) << std::endl;
        }
        for(std::int64_t j = 0; j < m_num_bye; j++){
            fout_gamma << "gamma_withMain_" << m_bye_name_vec[j] << "\t" << gamma_Mat_vec[j](1, 1) << std::endl;
        }
        fout_gamma << "n_random_used\t" << num_random_snp_used << std::endl;
        fout_gamma.close();
    }

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
    spdlog::info("Genome-wide GxE scan; writing results to {}.assoc", m_out_file);
    // Main result: SNP main effect + the three combined p-values.
    fout = open_output_stream(m_out_file + ".assoc");
    fout << "order\tchrom\tSNP\tbase\tallele1\tallele2\taf\tN\tbeta_main\tse_main\tp_main\tp_single\tp_multi\tp_gxe" << std::endl;

    // Single-environment interaction tests (one beta/se/p triple per environment)
    // go to a separate <out>.perE.assoc.
    ofstream fout_single = open_output_stream(m_out_file + ".perE.assoc");
    fout_single << "order\tchrom\tSNP\tbase\tallele1\tallele2\taf\tN";
    for(std::int64_t j = 0; j < m_num_bye; j++){
        fout_single << "\tbeta_" << m_bye_name_vec[j] << "\tse_" << m_bye_name_vec[j] << "\tp_" << m_bye_name_vec[j];
    }
    fout_single << std::endl;

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
        // Convert calibrated variances to z-scores (effect / SE) for column 0
        // (main effect) and columns 1.. (per-environment interactions).
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

                // Refit each environment's interaction (exact or approximate) only
                // when the cheap approximation already looks significant.
                if(p_Vec.minCoeff() < p_approx_cut){
                    Eigen::Map<Eigen::VectorXd> snpk(
                        snp_mat_part.data() + k * bed_context.num_id_used, bed_context.num_id_used);
                    double snpk_norm2 = geno_var_Vec(k);
                    MatrixXd snpkME = m_bye_mat.array().colwise() * snpk.array();
                    MatrixXd xmatk = MatrixXd::Zero(bed_context.num_id_used, 2);
                    xmatk.col(0) = snpk;

                    for(std::int64_t m = 0; m < m_num_bye; m++){
                            xmatk.col(1) = snpkME.col(m);

                            // Fast approximate refit (gamma-calibrated information matrix).
                            MatrixXd beta_Vec_var = (gamma_Mat_vec[m] * snpk_norm2).inverse();
                            double beta = (beta_Vec_var * (xmatk.transpose() * Vi0y))(1);
                            double beta_var = beta_Vec_var(1, 1);
                            double p_val = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);

                            // Adaptive exact refit: recompute x'V^-1 x exactly for this
                            // SNP-environment when the approximate hit is significant enough
                            // (p < --exact-cut; exact_cut == 1 forces it for all promising SNPs).
                            if(p_val < exact_cut){
                                MatrixXd design_cross_product = xmatk.transpose() * m_Vi0 * xmatk;
                                MatrixXd design_cross_inverse = design_cross_product.inverse();
                                VectorXd beta_Vec = design_cross_inverse * (xmatk.transpose() * Vi0y);
                                beta_var = design_cross_inverse(1, 1);
                                beta = beta_Vec(1);
                                p_val = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);
                            }

                            beta_Mat(k, m+1) = beta;
                            se_Mat(k, m+1) = sqrt(beta_var);
                            z_mat(k, m+1) = beta_Mat(k, m+1) / se_Mat(k, m+1);
                            p_Vec(m+1) = p_val;
                    }

                }

                // p_single: ACAT combination of the per-environment GxE p-values.
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

        // Combined GxE score = sum of squared interaction z-scores across
        // environments; its null is a chi-square mixture weighted by the
        // environment-correlation eigenvalues, evaluated with saddle().
        VectorXd score_Vec = z_mat.rightCols(m_num_bye).array().square().matrix().rowwise().sum();

        #pragma omp parallel for schedule(dynamic)
        for(std::int64_t k = 0; k < partition.num_snp_read; k++){
            double score_val = score_Vec(k);
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut) && score_val >= 0.0 && std::isfinite(score_val)){

                // p_multi: saddlepoint p-value of the combined score.
                double p1 = saddle(score_val, corrE_eigenvalues);
                p_combined_Mat(k, 1) = p1;

                // p_gxe: ACAT combination of p_single and p_multi.
                VectorXd p_Vec(2);
                p_Vec(0) = p_combined_Mat(k, 0);
                p_Vec(1) = p1;
                p_combined_Mat(k, 2) = acat_combine_pvalues(p_Vec)(1);
            }else{
                p_combined_Mat(k, 1) = p_combined_Mat(k, 2) = NAN;
            }
        }
        
        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        fastGxE::output_GxE(fout, fout_single, snp_info_vec, freq_arr, nobs_geno_arr,
            beta_Mat, se_Mat, p_Mat, p_combined_Mat,
            partition.start_snp, partition.num_snp_read);
    }
    fout.close();
    fout_single.close();
}

// Writes the genome-wide GxE results, split across two files. The main stream
// `fout` (<out>.assoc) gets the SNP annotation, the SNP main effect (beta/se/p)
// and the three combined p-values (p_single, p_multi, p_gxe). The `fout_single`
// stream (<out>.perE.assoc) gets the SNP annotation plus every single-
// environment interaction (beta/se/p, columns 1..num_env of the matrices) side
// by side. Rows are buffered per stream and flushed once for throughput.
void fastGxE::output_GxE(ofstream& fout, ofstream& fout_single,
            const vector<string>& snp_info_vec, const VectorXd& freq_arr, const VectorXd& nobs_geno_arr,
            const MatrixXd& beta_mat, const MatrixXd& se_mat, const MatrixXd& p_mat, const MatrixXd& p_combined_Mat,
            std::int64_t start_snp, std::int64_t num_snp_read){

    const std::int64_t num_env = beta_mat.cols() - 1;  // column 0 is the SNP main effect
    std::stringstream ss_main;
    std::stringstream ss_single;
    ss_main.precision(8);
    ss_single.precision(8);

    for (std::int64_t i = 0; i < num_snp_read; i++) {
        const long long n_i = static_cast<long long>(nobs_geno_arr(i));

        // Main file: SNP-identifying prefix + main effect + combined p-values.
        ss_main << start_snp + i << "\t" << snp_info_vec[i] << "\t" << freq_arr(i) << "\t" << n_i
                << "\t" << beta_mat(i, 0) << "\t" << se_mat(i, 0) << "\t" << p_mat(i, 0);
        for (std::int64_t m = 0; m < p_combined_Mat.cols(); m++) {
            ss_main << "\t" << p_combined_Mat(i, m);
        }
        ss_main << "\n";

        // Single-environment file: same prefix + every environment's interaction.
        ss_single << start_snp + i << "\t" << snp_info_vec[i] << "\t" << freq_arr(i) << "\t" << n_i;
        for (std::int64_t j = 0; j < num_env; j++) {
            ss_single << "\t" << beta_mat(i, j + 1) << "\t" << se_mat(i, j + 1) << "\t" << p_mat(i, j + 1);
        }
        ss_single << "\n";
    }

    fout << ss_main.str();
    fout_single << ss_single.str();
}


// ===== Top-level driver =====

// Program entry point for the fastGxE class. Parses and validates CLI options,
// sets the MKL/OpenMP thread counts, resolves the requested columns, runs the
// shared data-preparation pipeline (pre_data + process_grm), fits the null-model
// variance components for the selected mode (varcom_main or varcom_GxE), and, if a
// .bed file is supplied, runs the corresponding genome-wide scan (test_main or
// test_GxE) over the resolved SNP range. Returns 0 on success.
int fastGxE::run(int argc, char **argv) {
    CLI::App app{"fastGxE -  Scalable and fast genome-wide multi-environment GxE method"};
    RunOptions options;
    configure_run_cli(app, options);
    CLI11_PARSE(app, argc, argv);
    // Echo the user-supplied options (no defaults, no descriptions) so the exact
    // invocation can be reproduced from the log.
    spdlog::info("Options:\n{}", app.config_to_str(false, false));
    validate_run_options(options);
    log_run_options(options);
    
    // Set MKL and OpenMP threads
    mkl_set_num_threads(options.threads);
    // Let MKL trim its thread count for small problems and, crucially, run
    // single-threaded when called from inside an OpenMP parallel region (the
    // per-SNP scan loops in test_GxE). This is the MKL default; setting it
    // explicitly guards against an environment that forces MKL_DYNAMIC=false,
    // which could otherwise oversubscribe to threads^2 threads.
    mkl_set_dynamic(1);
    omp_set_num_threads(options.threads);

    // Resolve all requested column names against the data-file header before any
    // expensive GRM or genotype work starts.
    ResolvedInputColumns resolved = resolve_input_columns(options);
    m_trait_name_vec = options.trait_vec;
    m_bye_name_vec = options.bye_vec;

    // prepare data
    this->pre_data(options.out_file, options.data_file, options.agrm_file,
            resolved.covariate_index_vec, resolved.class_index_vec,
            resolved.bye_index_vec, resolved.trait_index_vec, options.missing_in_data_vec);

    // Optional log-transform is applied to the RAW phenotype (before any
    // covariate/environment correction), for both analysis modes.
    if(options.log_transform){
        spdlog::info("Applying natural-log transform to the raw phenotype");
        apply_log_transform(m_y);
    }

    VectorXd init_varcom;
    this->process_grm(options.test_main);
    if(options.test_main){
        // Main-effect mode keeps covariates in the model, so the inverse-normal
        // transform (if requested) is applied directly to the phenotype here.
        if(options.inv_normal){
            spdlog::info("Applying inverse-normal transform to the phenotype");
            apply_inverse_normal_transform(m_y);
        }
        this->varcom_main(init_varcom, options.maxiter, options.cc_gra, options.cc_gra, options.cc_logL);
    }else if(options.test_gxe){
        // For GxE the inverse-normal transform is applied inside pre_data_GxE,
        // AFTER the covariate/environment correction (see its inv_normal arg).
        this->pre_data_GxE(options.standardize_env, true, options.inv_normal);  // Phenotype correction is always required here.
        this->varcom_GxE(init_varcom, options.no_noisebye, options.maxiter, options.cc_par, options.cc_gra, options.cc_logL);
    }
    
    
    if(options.bed_file.empty()){
        spdlog::info("No --bfile provided; variance components estimated only (no genome-wide scan).");
        spdlog::info("Analysis completed.");
        return 0;
    }

    vector<std::int64_t> tmp_vec = this->get_start_end_pos(options.bed_file, options.snp_range_vec, options.split_task_vec);
    std::int64_t start_pos = tmp_vec[0], end_pos = tmp_vec[1];
    if(end_pos - start_pos <= options.npart_snp) options.npart_snp = 1;
    spdlog::info("Scanning SNPs [{}, {}) ({} SNPs) in {} partition(s)",
            start_pos, end_pos, end_pos - start_pos, options.npart_snp);

    if(options.test_main){
        this->test_main(options.bed_file, start_pos, end_pos, options.npart_snp,
                options.keep_random,
                options.num_random_snp, options.maf_cut, options.missing_rate_cut);
    }else if(options.test_gxe){
        this->test_GxE(options.bed_file, start_pos, end_pos, options.npart_snp,
                                 options.keep_random,
                                 options.exact_cut, options.num_random_snp, options.p_approx_cut,
                                 options.maf_cut, options.missing_rate_cut);
    }
    spdlog::info("Analysis completed.");
    return 0;
}


// Determines the [start, end) SNP index range to scan. Defaults to the whole .bed
// file, or honors a --snp-range (resolved from SNP names via the .bim) or a
// --split-task partition (total parts, current part); the two are mutually
// exclusive. Updates the output prefix accordingly and aborts on an invalid range.
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
