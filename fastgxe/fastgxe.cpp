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
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <string>
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
#include "fastgxe.hpp"


using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::endl;
using Eigen::HouseholderQR;
using Eigen::ColPivHouseholderQR;

namespace {

std::vector<std::string> read_header_columns(const std::string& data_file) {
    std::ifstream fin(data_file);
    if (!fin.is_open()) {
        spdlog::error("Failed to open data file: {}", data_file);
        exit(1);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        spdlog::error("Failed to read the head line from file: {}", data_file);
        exit(1);
    }

    process_line(line);
    if (line.empty()) {
        spdlog::error("Head line is empty or starts with #");
        exit(1);
    }

    return split_string(line);
}

std::string join_or_none(const std::vector<std::string>& values,
                         const char* delimiter = ", ") {
    return values.empty() ? "None" : join_string(values, delimiter);
}

struct RunOptions {
    int threads = 10;
    bool test_main = false;
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

using GrmTriplet = Eigen::Triplet<double>;

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
        spdlog::error("Eigen decomposition failed for GRM block starting at index {}.", start_index);
        exit(1);
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
        spdlog::error(error_message);
        exit(1);
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

MatrixXd build_gxe_relationship_block(const MatrixXd& grm_block, const MatrixXd& bye_block, std::int64_t num_bye,
        bool is_singleton_block) {
    if (is_singleton_block) {
        MatrixXd gxe_block = (bye_block.cwiseProduct(bye_block)).rowwise().sum() / num_bye;
        return (gxe_block.cwiseProduct(grm_block)).eval();
    }

    MatrixXd gxe_block = bye_block * bye_block.transpose() / num_bye;
    return (gxe_block.cwiseProduct(grm_block)).eval();
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
        spdlog::error("Failed to factorize covariance block starting at index {}.", start_index);
        exit(1);
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
        spdlog::error("Fail to open the output file: {}.var", out_file);
        exit(1);
    }

    for(std::int64_t i = 0; i < varcom.size(); ++i){
        fout << varcom(i) << " " << se_varcom(i) << std::endl;
    }
}

void write_variance_components(const std::string& out_file, const VectorXd& varcom) {
    ofstream fout(out_file + ".var");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file: {}.var", out_file);
        exit(1);
    }

    fout << varcom << std::endl;
}

MainModelIterationState evaluate_main_model_iteration(const VectorXd& grm_eigenvals,
        const MatrixXd& xmat_trans, const MatrixXd& y_trans, const VectorXd& varcom) {
    MainModelIterationState iteration_state;

    const VectorXd vmat = (1.0 / (grm_eigenvals.array() * varcom(0) + varcom(1))).matrix();
    iteration_state.logL = -(vmat.array().log()).sum();

    Eigen::ArrayXXd vxmat_arr = xmat_trans.array();
    vxmat_arr.colwise() *= vmat.array();
    const MatrixXd vxmat = vxmat_arr.matrix();

    MatrixXd xvxmat = xmat_trans.transpose() * vxmat;
    CustomLLT llt_solver;
    llt_solver.compute(xvxmat);
    iteration_state.logL += llt_solver.logDeterminant();
    xvxmat = llt_solver.inverse();

    VectorXd xvy = vxmat.transpose() * y_trans;
    const VectorXd py = vmat.cwiseProduct(y_trans - xmat_trans * (xvxmat * xvy));
    iteration_state.logL += (y_trans.transpose() * py).sum();

    const VectorXd dpy = grm_eigenvals.cwiseProduct(py);
    xvy = vxmat.transpose() * dpy;
    const VectorXd pdpy = vmat.cwiseProduct(dpy - xmat_trans * (xvxmat * xvy));
    xvy = vxmat.transpose() * py;
    const VectorXd ppy = vmat.cwiseProduct(py - xmat_trans * (xvxmat * xvy));

    iteration_state.fd_mat(0) = (vmat.cwiseProduct(grm_eigenvals)).sum();
    vxmat_arr.colwise() *= grm_eigenvals.array();
    iteration_state.fd_mat(0) -= (xvxmat.cwiseProduct(vxmat_arr.matrix().transpose() * vxmat)).sum();
    iteration_state.fd_mat(0) -= (py.transpose() * dpy).sum();
    iteration_state.fd_mat(0) *= -0.5;

    iteration_state.fd_mat(1) = vmat.sum() - (xvxmat.cwiseProduct(vxmat.transpose() * vxmat)).sum();
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
        spdlog::error(error_template, num_used);
        exit(1);
    }
}

ofstream open_output_stream(const string& file_path) {
    ofstream fout(file_path);
    if(!fout.is_open()){
        spdlog::error("Fail to open output file: {}", file_path);
        exit(1);
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
    )");

    app.add_option("-p,--threads", options.threads,
        "Number of threads to use (default: 10).")
        ->default_val(10);

    auto* test_main_flag = app.add_flag("--test-main", options.test_main,
        "Enable testing of SNP main effects (mutually exclusive with --test-gxe).");
    auto* test_gxe_flag = app.add_flag("--test-gxe", options.test_gxe,
        "Enable testing of SNP-environment interactions (mutually exclusive with --test-main).");
    test_main_flag->excludes(test_gxe_flag);
    test_gxe_flag->excludes(test_main_flag);

    app.add_flag("--no-noisebye", [&options](int count) {
        if (count > 0) options.no_noisebye = true;
    }, "Disable noise-by-environment interaction terms.\n"
       "  - Noise-by-environment interactions are ENABLED by default.\n"
       "  - Use --no-noisebye to turn it OFF.");

    app.add_option("--data", options.data_file, "Path to input data file (required).")->required();
    app.add_option("--trait", options.trait_vec, "Trait to analyze (required, expects 1 value).")
        ->expected(1)->required();
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
    if (options.test_main == options.test_gxe) {
        spdlog::error("Exactly one of --test-main or --test-gxe must be specified.");
        exit(1);
    }
    if (options.test_gxe && options.bye_vec.empty()) {
        spdlog::error("--env-int is required when --test-gxe is specified.");
        exit(1);
    }
    if (options.threads <= 0) {
        options.threads = 10;
    }
}

void log_run_options(const RunOptions& options) {
    spdlog::info("=== Parsed Arguments ===");
    spdlog::info("Threads: {}", options.threads);
    spdlog::info("Input Data File: {}", options.data_file);
    spdlog::info("Output File: {}", options.out_file);
    spdlog::info("GRM File: {}", options.agrm_file);
    spdlog::info("BED File: {}", options.bed_file.empty() ? "Not Provided" : options.bed_file);
    spdlog::info("Trait: {}", options.trait_vec[0]);
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
        spdlog::error("Trait names not found in the header: {}", join_string(missing_names));
        exit(1);
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
        spdlog::error("Interacting environments not found in the header: {}", join_string(missing_names));
        exit(1);
    }
    missing_names.clear();

    if (!resolved.bye_index_vec.empty()) {
        const std::vector<std::int64_t> covar_overlap = find_index(options.covariate_vec, options.bye_vec, missing_names);
        const std::vector<std::int64_t> class_overlap = find_index(options.class_vec, options.bye_vec, missing_names);
        if (!covar_overlap.empty() || !class_overlap.empty()) {
            spdlog::error("Interacting environments should not be included in --covar or --class. They are treated as covariates by default.");
            exit(1);
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
            spdlog::error("Duplicated sample IDs exist in the data file!");
            exit(1);
        }
    }

    if (sample_ids.empty()) {
        spdlog::error("No samples remain after aligning the phenotype data with the GRM IDs.");
        exit(1);
    }

    return sample_ids;
}

vector<fastGxE::GroupRecord> fastGxE::load_and_sort_group_records(const string& agrm_file,
        const vector<string>& retained_sample_ids, std::int64_t expected_num_samples) const {
    spdlog::info("Read the group files: {}.grm.group", agrm_file);
    std::unordered_set<string> retained_id_set(retained_sample_ids.begin(), retained_sample_ids.end());

    ifstream fin(agrm_file + ".grm.group");
    if(!fin.is_open()){
        spdlog::error("Fail to open the {}.grm.group", agrm_file);
        exit(1);
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
            spdlog::error("Malformed record in {}.grm.group: {}", agrm_file, line);
            exit(1);
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
        spdlog::error("The number of retained records in {}.grm.group ({}) does not match the retained phenotype sample count ({}).",
                      agrm_file, group_records.size(), expected_num_samples);
        exit(1);
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
            spdlog::error("GRM index {} in group file is out of range.", record.grm_index);
            exit(1);
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

    std::ifstream fin(agrm_file + ".agrm.sp.bin", std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}.agrm.sp.bin", agrm_file);
        exit(1);
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
        m_y = (m_y - m_xmat * design_qr.solve(m_y)).eval();
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
        
        MatrixXd ViX = vmat * m_xmat;
        MatrixXd XViX = m_xmat.transpose() * ViX;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XViX);
        Eigen::VectorXd diag = ldlt.vectorD();
        double XViX_logdet = diag.array().log().sum();
        MatrixXd XViXi = XViX.inverse();
        
        MatrixXd XViKViX = ViX.transpose() * m_grm_mat * ViX;
        MatrixXd Py = vmat * m_y - ViX * (XViXi * (ViX.transpose() * m_y));
        const double logL = inverse_result.logdet + XViX_logdet + (m_y.transpose() * Py).sum();
        spdlog::info("-2logL: {}", logL);
        MatrixXd KPy = m_grm_mat * Py;
        double tr_ViK = (vmat.cwiseProduct(m_grm_mat)).sum();
        double tr_XViXi_XViKViX = (XViXi.cwiseProduct(XViKViX)).sum();
        fd_mat(0) = tr_ViK - tr_XViXi_XViKViX - (Py.transpose() * KPy).sum();
        KPy_vec[0] = KPy;

        double tr_ViKe = (vmat.cwiseProduct(gxe_relationship.sparse_matrix)).sum();
        MatrixXd XViKeViX = ViX.transpose() * gxe_relationship.sparse_matrix * ViX;
        double tr_XViXi_XViKeViX = (XViXi.cwiseProduct(XViKeViX)).sum();
        MatrixXd KePy = gxe_relationship.sparse_matrix * Py;
        fd_mat(1) = tr_ViKe - tr_XViXi_XViKeViX - (Py.transpose() * KePy).sum();
        KPy_vec[1] = KePy;

        if(!no_noisebye){
            double tr_ViNxe = (vmat.cwiseProduct(nxe_diag)).sum();
            MatrixXd XViNxeViX = ViX.transpose() * nxe_diag * ViX;
            double tr_XViXi_XViNxeViX = (XViXi.cwiseProduct(XViNxeViX)).sum();
            MatrixXd NxePy = nxe_diag * Py;
            fd_mat(2) = tr_ViNxe - tr_XViXi_XViNxeViX - (Py.transpose() * NxePy).sum();
            KPy_vec[2] = NxePy;
        }



        fd_mat(num_cov - 1) = vmat.diagonal().sum() - (XViXi.cwiseProduct(ViX.transpose() * ViX)).sum() - (Py.transpose() * Py).sum();
        fd_mat *= -0.5;
        KPy_vec[num_cov - 1] = Py;

        vector <MatrixXd> PKPy_vec(num_cov);
        for(std::int64_t m = 0; m < num_cov; m++){
            const MatrixXd& KPy_tmp = KPy_vec[m];
            PKPy_vec[m] = vmat * KPy_tmp - ViX * (XViXi * (ViX.transpose() * KPy_tmp));
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

        // EM
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
    MatrixXd ViX = m_Vi0 * m_xmat;
    MatrixXd XViX = m_xmat.transpose() * ViX;

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XViX);
    Eigen::VectorXd diag = ldlt.vectorD();
    double logL = m_V0_logdet + diag.array().log().sum();
    MatrixXd XViXi = XViX.inverse();
    
    MatrixXd Py = this->m_Vi0 * m_y - ViX * (XViXi * (ViX.transpose() * m_y));
    logL += (m_y.transpose() * Py).sum();
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

    const VectorXd vmat0 = 1.0 / (m_grm_eigenvals.array() * m_varcom_null(0) + m_varcom_null(1));
    Eigen::ArrayXXd vxmat_arr0 = m_xmat_trans.array();
    vxmat_arr0.colwise() *= vmat0.array();
    const MatrixXd vxmat0 = vxmat_arr0.matrix();
    const MatrixXd xvxmat_inv0 = (m_xmat_trans.transpose() * vxmat0).inverse();
    const MatrixXd xvxmat_inv_xvmat0 = xvxmat_inv0 * vxmat0.transpose();
    
    spdlog::info("Randomly select SNPs and calculate the gamma");
    const vector<std::int64_t> random_snp_index_vec = sample_random_snp_indices(bed_context.num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
    load_centered_random_snp_panel(
        GenoA, bed_file, random_snp_index_vec, bed_context.id_index_in_bed_vec,
        snp_mat_part_random, freq_arr, missing_rate_arr, nobs_geno_arr);
    snp_mat_part_random = (m_grm_eigenvecs.transpose() * snp_mat_part_random).eval();

    VectorXd gamma_correct_Vec(num_random_snp);
    VectorXd p_Vec(num_random_snp);
    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t k = 0; k < num_random_snp; k++){
        VectorXd isnp_arr = snp_mat_part_random.col(k);
        // !score test
        VectorXd psnp = vmat0.cwiseProduct(isnp_arr - m_xmat_trans * (xvxmat_inv_xvmat0 * isnp_arr));
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
    
    const VectorXd pyarr0 = m_grm_eigenvecs * (vmat0.cwiseProduct(m_y_trans) - vxmat0 * (xvxmat_inv0 * (vxmat0.transpose() * m_y_trans)));
    const double* pyarr0_pt = pyarr0.data();

    const std::int64_t initial_partition_size = get_snp_partition(start_pos, end_pos, npart_snp, 0).num_snp_read;
    std::vector<double> snp_mat_part(initial_partition_size * bed_context.num_id_used);
    std::vector<double> geno_var_arr(initial_partition_size);
    std::vector<double> se_arr(initial_partition_size);
    std::vector<double> prod_arr(initial_partition_size);
    std::vector<double> eff_arr(initial_partition_size);
    VectorXd p_arr = VectorXd::Constant(initial_partition_size, std::numeric_limits<double>::quiet_NaN());
    
    for(std::int64_t i = 0; i < npart_snp; i++){
        const SnpPartition partition = get_snp_partition(start_pos, end_pos, npart_snp, i);
        if(i == npart_snp - 1){
            snp_mat_part.resize(partition.num_snp_read * bed_context.num_id_used);
            geno_var_arr.resize(partition.num_snp_read);
            se_arr.resize(partition.num_snp_read);
            prod_arr.resize(partition.num_snp_read);
            eff_arr.resize(partition.num_snp_read);
            p_arr = VectorXd::Constant(partition.num_snp_read, std::numeric_limits<double>::quiet_NaN());
        } else if (partition.num_snp_read != static_cast<std::int64_t>(p_arr.size())) {
            snp_mat_part.resize(partition.num_snp_read * bed_context.num_id_used);
            geno_var_arr.resize(partition.num_snp_read);
            se_arr.resize(partition.num_snp_read);
            prod_arr.resize(partition.num_snp_read);
            eff_arr.resize(partition.num_snp_read);
            p_arr = VectorXd::Constant(partition.num_snp_read, std::numeric_limits<double>::quiet_NaN());
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}", 
             i + 1, npart_snp, partition.start_snp, partition.num_snp_read);
        
        // Each partition reuses the same work buffers; only the logical SNP span
        // changes as we stream the genotype matrix block by block.
        GenoA.read_bed_centered_to_buffer_omp(bed_file + ".bed", snp_mat_part.data(), freq_arr, missing_rate_arr, nobs_geno_arr,
                partition.start_snp, partition.num_snp_read, bed_context.id_index_in_bed_vec);
        
        // Test
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            geno_var_arr[j] = cblas_ddot(
                bed_context.num_id_used, snp_mat_part.data() + j * bed_context.num_id_used, 1,
                snp_mat_part.data() + j * bed_context.num_id_used, 1);
        }

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            se_arr[j] = 1.0 / (gamma_correct * geno_var_arr[j]);
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, partition.num_snp_read, bed_context.num_id_used, 1.0, snp_mat_part.data(), 
                bed_context.num_id_used, pyarr0_pt, 1, 0.0, prod_arr.data(), 1);
        
        vdMul(partition.num_snp_read, se_arr.data(), prod_arr.data(), eff_arr.data());

        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < partition.num_snp_read; ++j) {
            se_arr[j] = sqrt(se_arr[j]);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for (std::int64_t k = 0; k < partition.num_snp_read; k++) {
            if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                double eff = eff_arr[k];
                double se = se_arr[k];
                double p_wald = gsl_cdf_chisq_Q(eff*eff/(se*se), 1);
                p_arr(k) = p_wald;
            }else{
                eff_arr[k] = NAN;
                se_arr[k] = NAN;
                p_arr(k) = NAN;
            }
        }
        const vector<string> snp_info_vec = GenoA.snp_anno(partition.start_snp, partition.num_snp_read);
        fastGxE::output(fout, snp_info_vec, freq_arr, missing_rate_arr,  
            eff_arr.data(), se_arr.data(), p_arr, partition.start_snp, partition.num_snp_read, p_cut);
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
            VectorXd Vi0X = this->m_Vi0 * snpk; // the phenotypes are corrected for fixed effects
            double XVi0Xi = 1 / (snpk.transpose() * Vi0X);
            VectorXd P0y = Vi0y - Vi0X * (XVi0Xi * (Vi0X.transpose() * m_y));

            VectorXd EoGtPy = EoG.transpose() * P0y;
            score = EoGtPy.transpose() * EoGtPy;
            score_GxEwithMain_Vec(i) = score;
            EoGtPEoG = EoG.transpose() * (this->m_Vi0 * EoG - Vi0X * (XVi0Xi * (Vi0X.transpose() * EoG)));
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
        spdlog::error("Fail to open the output file: {}.random.res", m_out_file);
        exit(1);
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
        spdlog::error("The number of used random SNPs to calculate gamma is: {}, which is less than 50% of random SNPs", 
                num_random_snp_used);
        exit(1);
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
        spdlog::error("Fail to open the output file: {}.res", m_out_file);
        exit(1);
    }
    fout << "order chrom SNP cm base allele1 allele2 af missing beta se p_main score p_gxe" << std::endl;

    // Allocate memory
    std::int64_t _omp_max_threads = omp_get_max_threads();

    double **EP0y_part = new double*[_omp_max_threads];
    for(int i = 0; i < _omp_max_threads; i++){
        EP0y_part[i] = new double[num_id_used];
    }

    double **WEP0y_part = new double*[_omp_max_threads];
    for(int i = 0; i < _omp_max_threads; i++){
        WEP0y_part[i] = new double[m_num_bye];
    }

    std::int64_t num_snp_read = (std::int64_t)((end_pos - start_pos)/npart_snp);
    double *snp_mat_part = new double[num_snp_read * num_id_used];
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
            delete [] snp_mat_part;
            snp_mat_part = new double[num_snp_read * num_id_used];
            beta_main_Vec.resize(num_snp_read);
            se_main_Vec.resize(num_snp_read);
            p_main_Vec.resize(num_snp_read);
            score_Vec.resize(num_snp_read);
            p_gxe_Vec.resize(num_snp_read);
        }
        
        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}", 
             ipart + 1, npart_snp, start_snp, num_snp_read);
        
        GenoA.read_bed_centered_to_buffer_omp(bed_file + ".bed", snp_mat_part, freq_arr, missing_rate_arr, nobs_geno_arr,
                start_snp, num_snp_read, id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for(std::int64_t k = 0; k < num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                std::int64_t thread_id = omp_get_thread_num();

                // snpk.T * snpk
                double snpk_norm2 = cblas_ddot(num_id_used, snp_mat_part + k*num_id_used, 1, 
                                    snp_mat_part + k*num_id_used, 1);

                // main
                double beta_main_var = 1 / (gamma_main * snpk_norm2);
                se_main_Vec(k) = sqrt(beta_main_var);
                double beta_main = cblas_ddot(num_id_used, snp_mat_part + k*num_id_used, 1, 
                                    Vi0y_pt, 1) * beta_main_var;
                beta_main_Vec(k) = beta_main;
                double p_main = gsl_cdf_chisq_Q(beta_main * beta_main / beta_main_var, 1);
                p_main_Vec(k) = p_main;

                // GxE without main effect
                vdMul(num_id_used, Vi0y_pt, snp_mat_part + k*num_id_used, EP0y_part[thread_id]);
                cblas_dgemv(CblasColMajor, CblasTrans, num_id_used, m_num_bye, 1.0, bye_mat_pt, num_id_used,
                            EP0y_part[thread_id], 1, 0.0, WEP0y_part[thread_id], 1);
                double score = cblas_ddot(m_num_bye, WEP0y_part[thread_id], 1, WEP0y_part[thread_id], 1);
                VectorXd chi_coef_Vec_kth = chi2_coef_GxEnoMain_Vec * snpk_norm2;
                double p_gxe = saddle(score, chi_coef_Vec_kth);

                if(p_main < p_approx_cut || p_gxe < p_approx_cut){
                    // GxE with main effect
                    Eigen::Map<Eigen::VectorXd> snpk(snp_mat_part + k * num_id_used, num_id_used);
                    VectorXd Vi0X = this->m_Vi0 * snpk;
                    double XVi0Xi = 1 / (snpk.transpose() * Vi0X);
                    VectorXd P0y = Vi0y - Vi0X * (Vi0X.transpose() * m_y) * XVi0Xi;
                    
                    if(speed == 0){
                        // exact test
                        MatrixXd EoG = this->m_bye_mat.array().colwise() * snpk.array();
                        VectorXd EoGtPy = EoG.transpose() * P0y;
                        score = EoGtPy.transpose() * EoGtPy;
                        MatrixXd EoGtPEoG = EoG.transpose() * (this->m_Vi0 * EoG - Vi0X * (XVi0Xi * (Vi0X.transpose() * EoG)));
                        
                        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
                        eigensolver.compute(EoGtPEoG);
                        chi_coef_Vec_kth = eigensolver.eigenvalues().real();
                        p_gxe = saddle(score, chi_coef_Vec_kth);
                    }else{
                        // approximate test
                        VectorXd EtGoPy = this->m_bye_mat.transpose() * (snpk.cwiseProduct(P0y));
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

    // free memory
    delete [] snp_mat_part;
    for(int i = 0; i < _omp_max_threads; i++){
        delete [] EP0y_part[i];
    }
    delete [] EP0y_part;

    for(int i = 0; i < _omp_max_threads; i++){
        delete [] WEP0y_part[i];
    }
    delete [] WEP0y_part;
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
                MatrixXd XVi0X = xmati.transpose() * m_Vi0 * xmati;
                MatrixXd XVi0Xi = XVi0X.inverse();
                VectorXd beta_Vec = XVi0Xi * (xmati.transpose() * Vi0y);
                double beta_var = XVi0Xi(1, 1);
                double beta = beta_Vec(1);
                double p = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);
                beta_Mat(i, j) = beta;
                se_Mat(i, j) = sqrt(beta_var);
                p_Mat(i, j) = p;
                gamma_Mat_2Dvec[j][i] = XVi0X / snp_norm2;
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
    corrE.diagonal() += Eigen::VectorXd::Constant(corrE.rows(), 0.001);

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(corrE);
    if (eigensolver.info() != Eigen::Success) {
        spdlog::error("Fail to eigen decomposition of GRM");
        exit(1);
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
                                MatrixXd XVi0X = xmatk.transpose() * m_Vi0 * xmatk;
                                MatrixXd XVi0Xi = XVi0X.inverse();
                                VectorXd beta_Vec = XVi0Xi * (xmatk.transpose() * Vi0y);
                                double beta_var = XVi0Xi(1, 1);
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

    // prepare data
    this->pre_data(options.out_file, options.data_file, options.agrm_file,
            resolved.covariate_index_vec, resolved.class_index_vec,
            resolved.bye_index_vec, resolved.trait_index_vec, options.missing_in_data_vec);
    
    VectorXd init_varcom;
    this->process_grm(options.test_main);
    if(options.test_main){
        this->varcom_main(init_varcom, options.maxiter, options.cc_gra, options.cc_gra, options.cc_logL);
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
        spdlog::error("Options --snp-range and --split-task cannot be used together.");
        exit(1);
    }else if(!snp_range_vec.empty()){
        vector<std::int64_t> tmp = GenoA.find_bim_index(snp_range_vec);
        start_pos = tmp[0];
        end_pos = tmp[1] + 1;
        this->reset_output_prefix(this->m_out_file, start_pos, end_pos);
    }else if(!split_task_vec.empty()){
        if(split_task_vec[0] < split_task_vec[1] || split_task_vec[1] <= 0){
            spdlog::error("The second number in --split-task must be greater than 0 and not exceed the first number.");
            exit(1);
        }
        std::int64_t num_snp_part = num_snp / split_task_vec[0];
        start_pos = num_snp_part * (split_task_vec[1] - 1);
        end_pos = num_snp_part * split_task_vec[1];
        if(split_task_vec[0] == split_task_vec[1]) end_pos = num_snp;
        this->reset_output_prefix(this->m_out_file, split_task_vec[0], split_task_vec[1]);
    }
    if(start_pos >= end_pos){
        spdlog::error("The start SNP position must not exceed the end SNP position.");
        exit(1);
    }
    vector<std::int64_t> vec = {start_pos, end_pos};
    return vec;
}
