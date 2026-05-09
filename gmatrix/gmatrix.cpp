#include <cstdint>

#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen.

#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>  // Include CLI11

#include "mkl.h"
#include "omp.h"

#include "gmatrix.hpp"

#include "string_utils.hpp"
#include "EigenMatrix_utils.hpp"
#include "partitionLowerTriangle.hpp"
#include "geno.hpp"
#include "../utils/fatal_error.hpp"


using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using Eigen::RowMajor;


/**
 * @brief Compute the row range handled by the current matrix partition.
 *
 * The result is stored in `m_start_id` and `m_num_id_read`, which define the
 * contiguous block of individuals covered by the current lower-triangular GRM
 * partition. A single-partition run is handled directly to avoid the extra
 * helper call on the common full-matrix path.
 *
 * @param npart Total number of output partitions.
 * @param ipart One-based index of the current partition.
 */
void Gmatrix::iid_part(std::int64_t npart, std::int64_t ipart) {
    if (m_num_id <= 0) {
        spdlog::error("Cannot partition individuals because m_num_id is {}", m_num_id);
        throw std::invalid_argument("Gmatrix::iid_part: m_num_id must be positive.");
    }

    if (npart == 1) {
        if (ipart != 1) {
            spdlog::error("Single-partition mode expects ipart=1, got {}", ipart);
            throw std::invalid_argument("Gmatrix::iid_part: invalid ipart for npart=1.");
        }

        m_start_id = 0;
        m_num_id_read = m_num_id;
    } else {
        partitionLowerTriangle(m_num_id, npart, ipart, m_start_id, m_num_id_read);
    }

    spdlog::info("Partition {}/{} -> Start row: {}, Row count: {}", ipart, npart, m_start_id, m_num_id_read);
    spdlog::info("GRM elements in partition: {}", (m_start_id + m_start_id + m_num_id_read + 1) * m_num_id_read / 2);
}



/**
 * @brief Write the IDs corresponding to the current partition to `<out>.id`.
 */
void Gmatrix::output_iid() {
    const std::string out_id_file = m_out_file + ".id";

    // Open the sidecar ID file for the current GRM block.
    std::ofstream fout(out_id_file, std::ios::out);
    if (!fout) {
        spdlog::error("Failed to open output file for individual IDs: {}", out_id_file);
        throw std::runtime_error("File open failed: " + out_id_file);
    }

    // Guard against accidental out-of-range access before writing.
    if (m_start_id < 0 || m_num_id_read < 0) {
        spdlog::error("Index out of range: start={} + num_read={} exceeds vector size={}",
                      m_start_id, m_num_id_read, m_id_in_geno_vec.size());
        throw std::out_of_range("Invalid ID range");
    }

    const std::size_t start_idx = static_cast<std::size_t>(m_start_id);
    const std::size_t count = static_cast<std::size_t>(m_num_id_read);
    const std::size_t total_ids = m_id_in_geno_vec.size();

    if (start_idx > total_ids || count > total_ids - start_idx) {
        spdlog::error("Index out of range: start={} + num_read={} exceeds vector size={}",
                      m_start_id, m_num_id_read, total_ids);
        throw std::out_of_range("Invalid ID range");
    }

    // Preserve the original one-ID-per-line write pattern.
    const std::size_t end_idx = start_idx + count;
    for (std::size_t i = start_idx; i < end_idx; ++i) {
        fout << m_id_in_geno_vec[i] << '\n';
    }

    spdlog::info("Successfully wrote {} individual IDs to {}", m_num_id_read, out_id_file);
}


/**
 * @brief Update the output prefix used by downstream writers.
 *
 * @param file Output prefix without the binary/text suffixes added later.
 */
void Gmatrix::set_out_file(const std::string& file) {
    m_out_file = file;
}


/**
 * @brief Write the lower triangle of the GRM in dense binary format.
 *
 * Output stores only the lower-triangular values in row-major traversal.
 *
 * @param gmat GRM block to write.
 * @param val Value added to the diagonal before serialization.
 */
void Gmatrix::out_gmat(const MatrixXd& gmat, const double val) {
    const std::string filename = m_out_file + ".bin";

    // Binary output matches the downstream GRM readers.
    std::ofstream fout(filename, std::ios::binary);
    if (!fout) {
        spdlog::error("Failed to open output file: {}", filename);
        throw std::runtime_error("File open failed: " + filename);
    }
    spdlog::info("Writing matrix to file: {}", filename);

    // Only the lower triangle is materialized on disk.
    for (std::int64_t i = 0; i < m_num_id_read; ++i) {
        const std::int64_t global_row = m_start_id + i;

        for (std::int64_t j = 0; j < global_row; ++j) {
            const double value = gmat(i, j);
            fout.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }

        const double diag_value = gmat(i, global_row) + val;
        fout.write(reinterpret_cast<const char*>(&diag_value), sizeof(double));
    }

    if (!fout) {
        spdlog::error("Failed while writing matrix to file: {}", filename);
        throw std::runtime_error("File write failed: " + filename);
    }

    fout.close();
    spdlog::info("Successfully wrote matrix to file: {}", filename);
}


/**
 * @brief Write GRM metadata to `<out>.N`.
 *
 * The file stores two whitespace-separated values:
 * `m_num_snp_used` and `m_scale`.
 */
void Gmatrix::out_num_snp_used() {
    std::string filename = m_out_file + ".N";

    // Metadata is written as plain text for easy inspection.
    std::ofstream fout(filename);
    if (!fout) {
        spdlog::error("Failed to open output file: {}", filename);
        throw std::runtime_error("File open failed: " + filename);
    }

    fout << m_num_snp_used << " " << m_scale << "\n";

    fout.close();
    spdlog::info("Successfully wrote SNP count and scale to file: {}", filename);
}


/**
 * @brief Compute the additive genomic relationship matrix.
 *
 * SNPs are streamed in chunks from disk, centered, optionally filtered by MAF,
 * and accumulated into the target GRM block.
 *
 * @param npart_snp Number of SNP chunks processed sequentially.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param missing_rate_cut Maximum allowed missing rate per SNP/feature.
 * @param code_type `1` for variance scaling, `2` for SNP-count scaling.
 * @return Lower-triangular GRM block for the current IID partition.
 */
MatrixXd Gmatrix::agmatrix(std::int64_t npart_snp, double maf, double missing_rate_cut, int code_type) {
    spdlog::info("Computing Additive Relationship Matrix");
    if (npart_snp <= 0) {
        fatal_error("--npart-snp must be a positive integer");
    }

    // The output stores the current row block against all preceding columns.
    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read);
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;
    const bool is_bed = (m_input_geno_fmt == 0);
    const bool full_grm_block = (m_start_id == 0 && m_num_id_read == m_num_id);
    const std::int64_t base_num_snp_read = m_num_snp / npart_snp;
    const std::int64_t trailing_num_snp = m_num_snp % npart_snp;

    // Keep the genotype reader open across SNP chunks.
    GENO GENO_on_disk(m_geno_file);

    for (std::int64_t i = 0; i < npart_snp; ++i) {
        std::int64_t num_snp_read = base_num_snp_read;
        const std::int64_t start_snp = i * base_num_snp_read;
        if (i == npart_snp - 1) {
            num_snp_read += trailing_num_snp;
        }
        if (num_snp_read <= 0) {
            continue;
        }

        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);

        // Read one SNP chunk plus per-SNP summary statistics.
        MatrixXd snp_mat_part;
        VectorXd freq_arr, missing_rate_arr, nobs_geno_arr;
        std::vector<std::int64_t> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, nobs_geno_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        // BED input is allele dosage, so the additive mean is `2p`.
        // Generic feature matrices are centered by their observed column mean.
        if (is_bed) {
            snp_mat_part.rowwise() -= 2 * freq_arr.transpose();
        } else {
            snp_mat_part.rowwise() -= freq_arr.transpose();
        }

        Eigen::ArrayXd scale_base;
        if (is_bed) {
            scale_base = 2.0 * freq_arr.array() * (1.0 - freq_arr.array());
        } else {
            scale_base = snp_mat_part.array().square().colwise().sum() / static_cast<double>(m_num_id);
        }

        // Low-MAF BED SNPs and high-missing SNPs are removed from scaling and accumulation.
        std::vector<Eigen::Index> local_filtered_indices;
        local_filtered_indices.reserve(static_cast<std::size_t>(num_snp_read));
        for (Eigen::Index it_snp = 0; it_snp < freq_arr.size(); ++it_snp) {
            const bool low_maf = is_bed && (freq_arr(it_snp) < maf);
            const bool high_missing = missing_rate_arr(it_snp) > missing_rate_cut;
            const bool zero_variance = (scale_base(it_snp) == 0.0);
            if (low_maf || high_missing || zero_variance) {
                m_num_snp_used--;
                local_filtered_indices.push_back(it_snp);
            }
        }

        VectorXd col_scale = VectorXd::Ones(num_snp_read);
        // `code_type == 1` keeps an explicit global scale factor.
        // `code_type == 2` rescales columns by the inverse square root variance.
        if (code_type == 1) {
            if (!local_filtered_indices.empty()) {
                scale_base(local_filtered_indices).setZero(); // Filtered SNPs have zero contribution to the global scale.
                col_scale(local_filtered_indices).setZero();
            }
            m_scale += scale_base.sum();
        } else {
            if (!local_filtered_indices.empty()) {
                scale_base(local_filtered_indices).setOnes(); // Filtered SNPs must not be scaled to avoid NaN/Inf weights.
            }
            col_scale = scale_base.sqrt().inverse().matrix();
            if (!local_filtered_indices.empty()) {
                col_scale(local_filtered_indices).setZero();
            }
        }

        scale_cols_inplace(snp_mat_part, col_scale);

        // Either form the full block or only the requested lower-triangular slice.
        if (full_grm_block) {
            mat += snp_mat_part * snp_mat_part.transpose();
        } else {
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }

    // Under SNP-count scaling, the denominator is simply the number of retained SNPs.
    if (code_type == 2) {
        m_scale = m_num_snp_used;
    }
    if (m_scale <= 0.0) {
        fatal_error("No SNPs/features remain after filtering; cannot scale the GRM");
    }

    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);

    return mat / m_scale;
}

/**
 * @brief Compute the dominance GRM using the biological genotype coding.
 *
 * Homozygous alternate genotypes are first recoded from `2` to `0`, then each
 * SNP is centered by `2pq`, where `p` is the reference-allele frequency and
 * `q = 1 - p`.
 *
 * @param npart_snp Number of SNP chunks processed sequentially.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param missing_rate_cut Maximum allowed missing rate per SNP.
 * @param code_type `1` for variance scaling, `2` for SNP-count scaling.
 * @return Lower-triangular GRM block for the current IID partition.
 */
MatrixXd Gmatrix::dgmatrix_as(std::int64_t npart_snp, double maf, double missing_rate_cut, int code_type){
    spdlog::info("Biological genotypic dominance relationship matrix");
    if(m_input_geno_fmt != 0){
        fatal_error("Input genotype file must be plink bed file given by --bfile");
    }
    if (npart_snp <= 0) {
        fatal_error("--npart-snp must be a positive integer");
    }

    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read);
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;
    const bool full_grm_block = (m_start_id == 0 && m_num_id_read == m_num_id);
    const std::int64_t base_num_snp_read = m_num_snp / npart_snp;
    const std::int64_t trailing_num_snp = m_num_snp % npart_snp;
    constexpr double homo_alt_tol = 1e-5;
    GENO GENO_on_disk(m_geno_file);

    for(std::int64_t i = 0; i < npart_snp; i++){
        std::int64_t num_snp_read = base_num_snp_read;
        const std::int64_t start_snp = i * base_num_snp_read;
        if(i == npart_snp - 1) {
            num_snp_read += trailing_num_snp;
        }
        if (num_snp_read <= 0) {
            continue;
        }
        
        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);
        
        MatrixXd snp_mat_part;
        VectorXd freq_arr;
        VectorXd missing_rate_arr, nobs_geno_arr;
        vector<std::int64_t> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, nobs_geno_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        // Recode homozygous alternate genotypes to the biological dominance coding.
        for(Eigen::Index m = 0; m < snp_mat_part.rows(); ++m){
            for(Eigen::Index n = 0; n < snp_mat_part.cols(); ++n){
                if(std::fabs(snp_mat_part(m, n) - 2.0) < homo_alt_tol){
                    snp_mat_part(m, n) = 0.0;
                }
            }
        }

        std::vector<Eigen::Index> local_filtered_indices;
        local_filtered_indices.reserve(static_cast<std::size_t>(num_snp_read));
        for (Eigen::Index it_snp = 0; it_snp < freq_arr.size(); ++it_snp) {
            const bool low_maf = freq_arr(it_snp) < maf;
            const bool high_missing = missing_rate_arr(it_snp) > missing_rate_cut;
            if (low_maf || high_missing) {
                freq_arr(it_snp) = 0.0;
                snp_mat_part.col(it_snp).setZero();
                m_num_snp_used--;
                local_filtered_indices.push_back(it_snp);
            }
        }

        Eigen::ArrayXd arr_2pq = 2 * freq_arr.array() * (1 - freq_arr.array());
        snp_mat_part.rowwise() -= arr_2pq.matrix().transpose();
        Eigen::ArrayXd dominance_scale_base = arr_2pq * (1 - arr_2pq);
        VectorXd col_scale = VectorXd::Ones(num_snp_read);

        if(code_type == 1){
            if (!local_filtered_indices.empty()) {
                dominance_scale_base(local_filtered_indices).setZero();
                col_scale(local_filtered_indices).setZero();
            }
            m_scale += dominance_scale_base.sum();
        }else{
            if (!local_filtered_indices.empty()) {
                dominance_scale_base(local_filtered_indices).setOnes();
            }
            col_scale = dominance_scale_base.sqrt().inverse().matrix();
            if (!local_filtered_indices.empty()) {
                col_scale(local_filtered_indices).setZero();
            }
        }
        scale_cols_inplace(snp_mat_part, col_scale);

        if(full_grm_block){
            mat += snp_mat_part * snp_mat_part.transpose();
        }else{
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }
    if(code_type == 2){
        m_scale = m_num_snp_used;
    }
    if (m_scale <= 0.0) {
        fatal_error("No SNPs remain after filtering; cannot scale the GRM");
    }
    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);
    return mat / m_scale;
}

/**
 * @brief Compute the dominance GRM using the breeding-value coding.
 *
 * Each genotype is converted to the standard dominance contrast:
 * `0 -> -2p^2`, `1 -> 2pq`, `2 -> -2q^2`.
 *
 * @param npart_snp Number of SNP chunks processed sequentially.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param missing_rate_cut Maximum allowed missing rate per SNP.
 * @param code_type `1` for variance scaling, `2` for SNP-count scaling.
 * @return Lower-triangular GRM block for the current IID partition.
 */
MatrixXd Gmatrix::dgmatrix_gs(std::int64_t npart_snp, double maf, double missing_rate_cut, int code_type){
    spdlog::info("Statistical breeding dominance relationship matrix");
    if(m_input_geno_fmt != 0){
        fatal_error("Input genotype file must be plink bed file given by --bfile");
    }
    if (npart_snp <= 0) {
        fatal_error("--npart-snp must be a positive integer");
    }

    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read);
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;
    const bool full_grm_block = (m_start_id == 0 && m_num_id_read == m_num_id);
    const std::int64_t base_num_snp_read = m_num_snp / npart_snp;
    const std::int64_t trailing_num_snp = m_num_snp % npart_snp;
    GENO GENO_on_disk(m_geno_file);

    for(std::int64_t i = 0; i < npart_snp; i++){
        std::int64_t num_snp_read = base_num_snp_read;
        const std::int64_t start_snp = i * base_num_snp_read;
        if(i == npart_snp - 1) {
            num_snp_read += trailing_num_snp;
        }
        if (num_snp_read <= 0) {
            continue;
        }
        
        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);

        MatrixXd snp_mat_part;
        VectorXd freq_arr;
        VectorXd missing_rate_arr, nobs_geno_arr;
        vector<std::int64_t> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, nobs_geno_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        // Convert additive dosages to the breeding-value dominance coding.
        for(Eigen::Index m = 0; m < snp_mat_part.rows(); ++m){
            for(Eigen::Index n = 0; n < snp_mat_part.cols(); ++n){
                if (snp_mat_part(m, n) < 0.5) {
                    snp_mat_part(m, n) = -2 * freq_arr(n) * freq_arr(n);
                }
                else if (snp_mat_part(m, n) > 1.5) {
                    snp_mat_part(m, n) = -2 * (1 - freq_arr(n)) * (1 - freq_arr(n));
                }
                else {
                    snp_mat_part(m, n) = 2 * (1 - freq_arr(n)) * freq_arr(n);
                }
            }
        }

        std::vector<Eigen::Index> local_filtered_indices;
        local_filtered_indices.reserve(static_cast<std::size_t>(num_snp_read));
        for (Eigen::Index it_snp = 0; it_snp < freq_arr.size(); ++it_snp) {
            const bool low_maf = freq_arr(it_snp) < maf;
            const bool high_missing = missing_rate_arr(it_snp) > missing_rate_cut;
            if (low_maf || high_missing) {
                freq_arr(it_snp) = 0.0;
                snp_mat_part.col(it_snp).setZero();
                m_num_snp_used--;
                local_filtered_indices.push_back(it_snp);
            }
        }

        Eigen::ArrayXd arr_2pq = 2 * freq_arr.array() * (1 - freq_arr.array());
        Eigen::ArrayXd dominance_scale_base = arr_2pq;
        VectorXd col_scale = VectorXd::Ones(num_snp_read);
        if(code_type == 1){
            if (!local_filtered_indices.empty()) {
                dominance_scale_base(local_filtered_indices).setZero();
                col_scale(local_filtered_indices).setZero();
            }
            m_scale += dominance_scale_base.square().sum();
        }else{
            if (!local_filtered_indices.empty()) {
                dominance_scale_base(local_filtered_indices).setOnes();
            }
            col_scale = dominance_scale_base.inverse().matrix();
            if (!local_filtered_indices.empty()) {
                col_scale(local_filtered_indices).setZero();
            }
        }
        scale_cols_inplace(snp_mat_part, col_scale);

        if(full_grm_block){
            mat += snp_mat_part * snp_mat_part.transpose();
        }else{
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }
    if(code_type == 2){
        m_scale = m_num_snp_used;
    }
    if (m_scale <= 0.0) {
        fatal_error("No SNPs remain after filtering; cannot scale the GRM");
    }

    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);

    return mat / m_scale;
}



/**
 * @brief Parse CLI arguments and dispatch GRM computation.
 */
int Gmatrix::run(int argc, char* argv[]) {
    CLI::App app{"Gmatrix - Genomic Relationship Matrix Tool"};
    const std::vector<std::string> supported_grm_types = {"agrm", "dgrm_as", "dgrm_gs"};

    app.description(R"(
Compute dense genomic relationship matrices (GRMs) and write lower-triangular
binary output to `<out>.grm.bin`, with companion `.grm.id` and `.grm.N` files.

Defaults:
    --grm-type agrm
    --code-type 2
    --maf 0.01
    --missing-rate 0.05

Quick start:
    Additive GRM from PLINK BED
        gmat --make-grm --bfile test --out test

    Additive GRM from feature matrix
        gmat --make-grm --ffile test --out test

    Biological dominance GRM
        gmat --make-grm --grm-type dgrm_as --bfile test --out test

    Statistical breeding dominance GRM
        gmat --make-grm --grm-type dgrm_gs --bfile test --out test

    Variance-scaled alternative
        gmat --make-grm --bfile test --code-type 1 --out test
    )");

    bool make_grm = false;
    app.add_flag("--make-grm", make_grm, "Calculate GRM");

    int threads = 10;
    app.add_option("--threads", threads, "Number of threads (default: 10)");

    // GRM type.
    std::string grm_type = "agrm";
    app.add_option("--grm-type", grm_type, "GRM type (agrm, dgrm_as, dgrm_gs)")
        ->default_val("agrm")
        ->check(CLI::IsMember(supported_grm_types));

    // Exactly one input source must be provided.
    std::string bfile, ffile;
    auto opt_bfile = app.add_option("--bfile", bfile, "Prefix for PLINK binary PED files");
    auto opt_ffile = app.add_option("--ffile", ffile, "Prefix for feature-by-sample matrix files");
    opt_bfile->excludes(opt_ffile);
    opt_ffile->excludes(opt_bfile);

    // Output prefix.
    app.add_option("--out", m_out_file, "Prefix for output file")->required();

    // Normalization strategy.
    int code_type = 2;
    app.add_option("--code-type", code_type, "Genotype coding type (1 or 2)")
        ->default_val(2)
        ->check(CLI::IsMember({1, 2}));

    // Output partition: total number of pieces and current one-based index.
    std::vector<std::int64_t> npart_values(2, 1);
    app.add_option("--npart", npart_values, "Partition GRM elements: <npart> <ipart>")
    ->expected(2)
    ->check(CLI::PositiveNumber)
    ->default_val(std::vector<std::int64_t>{1, 1});

    // Number of SNP chunks streamed from disk.
    std::int64_t npart_snp = 20;
    app.add_option("--npart-snp", npart_snp, "Partition SNPs into n parts")
        ->default_val(20)
        ->check(CLI::PositiveNumber);

    // MAF cutoff for BED-based GRMs.
    double maf = 0.01;
    app.add_option("--maf", maf, "MAF cutoff value")->default_val(0.01);

    // Maximum allowed missing rate per SNP/feature.
    double missing_rate_cut = 0.05;
    app.add_option("--missing-rate", missing_rate_cut, "Maximum allowed missing rate")->default_val(0.05);

    // Diagonal jitter added during serialization.
    double val = 0.001;
    app.add_option("--val", val, "Add value to diagonal to ensure positive matrix")->default_val(0.001);

    // String tokens treated as missing values in text-based inputs.
    std::vector<std::string> missing_geno = {
        "NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan",
        "<NA>", "<na>", "N/A", "n/a"
    };
    app.add_option("--missing-geno", missing_geno, "Missing genotype values (space-separated)")
        ->expected(-1)
        ->default_val(missing_geno);

    CLI11_PARSE(app, argc, argv);

    if (!make_grm) {
        spdlog::error("No task selected. Please pass --make-grm.");
        return 1;
    }

    spdlog::info("Using {} threads", threads);
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);
    
    // Persist parsed missing-value markers for later genotype reads.
    m_missing_in_geno_vec = missing_geno;

    // Infer the on-disk genotype format from the chosen input flag.
    if (!bfile.empty()) {
        m_geno_file = bfile;
        m_input_geno_fmt = 0;
    } else if (!ffile.empty()) {
        m_geno_file = ffile;
        m_input_geno_fmt = 1;
    } else {
        spdlog::error("One of --bfile or --ffile must be provided.");
        return 1;
    }

    // Validate the requested output partition.
    const std::int64_t npart = npart_values[0];
    const std::int64_t ipart = npart_values[1];
    if(ipart <= 0 || ipart > npart){
        spdlog::error("Check the argument for --npart");
        return 1;
    }

    // Log the final runtime configuration.
    spdlog::info("GRM Type: {}", grm_type);
    spdlog::info("Genotype File: {}", m_geno_file);
    spdlog::info("Output File: {}", m_out_file);
    spdlog::info("Partitions: npart={}, ipart={}", npart, ipart);
    spdlog::info("SNP Partitions: {}", npart_snp);
    spdlog::info("MAF Cutoff: {}", maf);
    spdlog::info("Missing Rate Cutoff: {}", missing_rate_cut);
    spdlog::info("Diagonal Addition Value: {}", val);
    spdlog::info("Missing Genotype Values: {}", join_string(m_missing_in_geno_vec, " "));

    // Read dataset dimensions and individual IDs once before streaming SNP chunks.
    GENO GENO_on_disk(m_geno_file);
    m_num_id = GENO_on_disk.get_num_iid();
    m_num_snp = GENO_on_disk.get_num_sid();
    m_id_in_geno_vec = GENO_on_disk.iid_vec();

    // Log the discovered genotype metadata.
    spdlog::info("Genotype file: {}", m_geno_file);
    spdlog::info("Output file: {}", m_out_file);
    spdlog::info("Input genotype format: {}", m_input_geno_fmt);
    spdlog::info("The number of individuals: {}", m_num_id);
    spdlog::info("The number of SNPs: {}", m_num_snp);

    this->iid_part(npart, ipart);
    MatrixXd gmat;

    if (grm_type == "agrm") {
        spdlog::info("Computing Additive Genomic Relationship Matrix");
        gmat = this->agmatrix(npart_snp, maf, missing_rate_cut, code_type);
    } else if (grm_type == "dgrm_as") {
        spdlog::info("Computing Dominance GRM (association)");
        gmat = this->dgmatrix_as(npart_snp, maf, missing_rate_cut, code_type);
    } else if (grm_type == "dgrm_gs") {
        spdlog::info("Computing Dominance GRM (Genomic Selection)");
        gmat = this->dgmatrix_gs(npart_snp, maf, missing_rate_cut, code_type);
    } else {
        spdlog::error("Invalid argument for --grm-type: {}", grm_type);
        return 1;
    }

    // Partitioned runs append `<npart>_<ipart>` before the GRM suffix.
    string out_file_name = m_out_file;
    if (npart > 1) {
        out_file_name += "." + std::to_string(npart) + "_" + std::to_string(ipart);
    }

    out_file_name += ".grm";
    this->set_out_file(out_file_name);
    spdlog::info("Output file set to: {}", out_file_name);

    this->output_iid();
    this->out_num_snp_used();
    this->out_gmat(gmat, val);

    return 0;
}
