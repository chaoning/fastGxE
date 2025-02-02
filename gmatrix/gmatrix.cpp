#define EIGEN_USE_MKL_ALL  // must be before include Eigen

#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>  // Include CLI11

#include "mkl.h"
#include "omp.h"

#include "gmatrix.hpp"

#include "../utils/string_utils.hpp"
#include "../utils/EigenMatrix_utils.hpp"
#include "../utils/partitionLowerTriangle.hpp"
#include "../utils/geno.hpp"


using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using Eigen::RowMajor;
using Eigen::Triplet;


/**
 * @brief Construct a new Gmatrix::Gmatrix object
 * 
 */
Gmatrix::Gmatrix() {
    m_num_snp_used = 0;
    m_scale = 0; 
    m_start_id = 0; 
    m_num_id_read = 0;  // Initialize correctly
    m_num_id = 0;       // Ensure this is initialized
    m_num_snp = 0;      // Ensure this is initialized
}


Gmatrix::~Gmatrix() {
}


/**
 * @brief Determines the partition range for individuals (m_start_id, m_num_id_read)
 * 
 * @param npart Number of partitions
 * @param ipart Current partition index (1-based)
 */
void Gmatrix::iid_part(long long npart, long long ipart) {
    spdlog::info("Subset: {}/{}", ipart, npart);
    partitionLowerTriangle(m_num_id, npart, ipart, 
                            m_start_id, m_num_id_read);
}



/**
 * @brief Output the individual IDs to a file
 */
void Gmatrix::output_iid() {
    std::string out_id_file = m_out_file + ".id"; // Output file for individuals' ID
    
    // Open file safely
    std::ofstream fout(out_id_file, std::ios::out);
    if (!fout) {
        spdlog::error("Failed to open output file for individuals' ID: {}", out_id_file);
        exit(1);
    }

    // Write IDs to file
    for (long long i = m_start_id; i <  m_start_id + m_num_id_read; i++) {
        fout << m_id_in_geno_vec[i] << '\n';  // '\n' is more efficient than std::endl
    }

    spdlog::info("Successfully wrote {} individual IDs to {}", m_num_id_read, out_id_file);
}


/**
 * @brief Sets the name of the output file.
 * 
 * This function updates the output file name used by the Gmatrix class.
 * 
 * @param file The new name of the output file.
 */
void Gmatrix::set_out_file(const std::string& file) {
    m_out_file = file;
}


/**
 * @brief Outputs the matrix in different formats (dense or sparse).
 * 
 * Writes the given matrix `gmat` to a binary file in either dense or sparse format.
 * 
 * @param gmat The matrix to be output.
 * @param sparse If 0, outputs as dense; otherwise, outputs as sparse format.
 * @param cut_val The threshold for sparsity; values below it are not written (except diagonals).
 * @param val The value to be added to diagonal elements before writing.
 */
void Gmatrix::out_gmat(const MatrixXd& gmat, const int sparse, const double cut_val, const double val) {
    // Determine output file name
    std::string filename = m_out_file + (sparse ? ".sp.bin" : ".bin");

    // Open file in binary mode
    std::ofstream fout(filename, std::ios::binary);
    if (!fout) {
        spdlog::error("Failed to open output file: {}", filename);
        exit(1);
    }
    spdlog::info("Writing matrix to file: {}", filename);


    // Writing process
    for (long long i = 0; i < m_num_id_read; i++) {
        long long ii = i + m_start_id;
        for (long long j = 0; j <= m_start_id + i; j++) {
            double value = gmat(i, j);
            if (ii == j) value += val;  // Add diagonal adjustment

            if (sparse == 0) {  
                // Dense format: Write directly
                fout.write(reinterpret_cast<const char*>(&value), sizeof(double));
            } else if (value > cut_val || ii == j) {  
                // Sparse format: Write only non-zero values
                fout.write(reinterpret_cast<const char*>(&ii), sizeof(long long));
                fout.write(reinterpret_cast<const char*>(&j), sizeof(long long));
                fout.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
    }

    // Close file
    fout.close();
    spdlog::info("Successfully wrote matrix to file: {}", filename);
}


/**
 * @brief Outputs the number of SNPs used and scale factor in binary format.
 */
void Gmatrix::out_num_snp_used() {
    std::string filename = m_out_file + ".N.bin";

    // Open file in binary mode
    std::ofstream fout(filename, std::ios::binary);
    if (!fout) {
        spdlog::error("Failed to open output file: {}", filename);
        exit(1);
    }

    // Write values to binary file
    fout.write(reinterpret_cast<const char*>(&m_num_snp_used), sizeof(double));
    fout.write(reinterpret_cast<const char*>(&m_scale), sizeof(double));

    fout.close();
    spdlog::info("Successfully wrote SNP count and scale to file: {}", filename);
}


/**
 * @brief Normalizes the weights from a file or defaults to 1 if no file is provided.
 * 
 * @param weight_arr Reference to the weight vector to be normalized.
 * @param weight_file Path to the weight file.
 */
void Gmatrix::normalized_weights(VectorXd& weight_arr, const std::string& weight_file) {
    if (!weight_file.empty()) {
        spdlog::info("Reading weight file: {}", weight_file);
        std::ifstream fin(weight_file);
        if (!fin) {
            spdlog::error("Failed to open weight file: {}", weight_file);
            exit(1);
        }

        std::vector<double> weights((std::istream_iterator<double>(fin)), std::istream_iterator<double>());
        fin.close();

        weight_arr = Eigen::Map<VectorXd>(weights.data(), weights.size());
    }

    if (weight_arr.size() == 0) {
        weight_arr.setOnes(m_num_snp);
    }

    if (weight_arr.size() != m_num_snp) {
        spdlog::error("Mismatch: weights ({}) vs SNPs ({})", weight_arr.size(), m_num_snp);
        exit(1);
    }

    if (weight_arr.minCoeff() <= 0) {
        spdlog::error("All weights must be greater than 0.");
        exit(1);
    }

    weight_arr = (m_num_snp * weight_arr.array() / weight_arr.sum()).sqrt();
    spdlog::info("Weights successfully normalized.");
}


/**
 * @brief Computes the additive genetic relationship matrix (GRM).
 * 
 * @param npart_snp Number of SNP partitions.
 * @param weight_arr SNP weight array.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param code_type Scaling type: 1 for variance-based, 2 for SNP count-based.
 * @return MatrixXd The computed GRM.
 */
MatrixXd Gmatrix::agmatrix(long long npart_snp, VectorXd& weight_arr, double maf, int code_type) {
    spdlog::info("Computing Additive Relationship Matrix");

    // Initialize G matrix
    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read);
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;

    // Load genotype data
    GENO GENO_on_disk(m_geno_file);

    for (long long i = 0; i < npart_snp; i++) {
        long long num_snp_read = m_num_snp / npart_snp;
        long long start_snp = i * num_snp_read;
        if (i == npart_snp - 1) {
            num_snp_read += m_num_snp % npart_snp;
        }

        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);

        // Read SNP data
        MatrixXd snp_mat_part;
        VectorXd freq_arr, missing_rate_arr;
        std::vector<long long> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        // Centralize SNP data (-2p)
        snp_mat_part.rowwise() -= 2 * freq_arr.transpose();

        // Apply MAF filter
        std::vector<long long> snp_index;
        if (m_input_geno_fmt != 2) {
            for (long long it_snp = 0; it_snp < freq_arr.size(); ++it_snp) {
                if (freq_arr(it_snp) < maf) {
                    freq_arr(it_snp) = 0;
                    weight_arr(it_snp + start_snp) = 0;
                    snp_mat_part.col(it_snp).setZero();
                    m_num_snp_used--;
                    snp_index.push_back(it_snp + start_snp);
                }
            }
        }

        // Scaling
        if (code_type == 1) {
            if (m_input_geno_fmt != 2) {
                m_scale += (2 * freq_arr.array() * (1.0 - freq_arr.array())).sum();
            } else {
                m_scale += (snp_mat_part.array().square().colwise().sum() / m_num_id).sum();
            }
        } else {
            Eigen::ArrayXd tmp_arr = weight_arr.segment(start_snp, num_snp_read);
            if (m_input_geno_fmt != 2) {
                tmp_arr /= (2 * freq_arr.array() * (1.0 - freq_arr.array())).sqrt();
            } else {
                tmp_arr /= (snp_mat_part.array().square().colwise().sum() / m_num_id).sqrt();
            }
            weight_arr.segment(start_snp, num_snp_read) = tmp_arr;
        }

        // Apply weighting
        for (auto idx : snp_index) {
            weight_arr(idx) = 0;
        }
        mat_row_elementwise_dot_vec(snp_mat_part, weight_arr.segment(start_snp, num_snp_read));

        // Compute G matrix
        if (m_start_id == 0 && m_num_id_read == m_num_id) {
            mat += snp_mat_part * snp_mat_part.transpose();
        } else {
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }

    // Final scaling
    if (code_type == 2) {
        m_scale = m_num_snp_used;
    }

    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);

    return mat / m_scale;
}

/**
 * @brief Computes the biological genotypic dominance relationship matrix.
 * 
 * @param npart_snp Number of SNP partitions.
 * @param weight_arr SNP weight array.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param code_type Scaling type: 1 for variance-based, 2 for SNP count-based.
 * @return MatrixXd The computed GRM.
 */
MatrixXd Gmatrix::dgmatrix_as(long long npart_snp, VectorXd& weight_arr, double maf, int code_type){
    spdlog::info("Biological genotypic dominance relationship matrix");
    if(m_input_geno_fmt != 0){
        spdlog::error("Input genotype file must be plink bed file given by --bfile");
        exit(1);
    }

    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read); // G matrix
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;
    GENO GENO_on_disk(m_geno_file);

    for(long long i = 0; i < npart_snp; i++){
        long long num_snp_read = (long long)(m_num_snp/npart_snp);
        long long start_snp = i*num_snp_read;
        if(i == npart_snp - 1)
            num_snp_read = num_snp_read + m_num_snp % npart_snp;
        
        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);
        
        MatrixXd snp_mat_part; // SNP matrix stored in the memory
        VectorXd freq_arr; // SNP frequence
        VectorXd missing_rate_arr;
        vector<long long> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        for(long long m = 0; m < m_num_id; m++){ // 2 -> 0
            for(long long n = 0; n < num_snp_read; n++){
                if(std::fabs(snp_mat_part(m, n) - 2.0) < 0.00001){
                    snp_mat_part(m, n) = 0.0;
                }
            }
        }


        vector <long long> snp_index;
        for(auto it = freq_arr.begin(); it != freq_arr.end(); it++){
            if(*it < 0.05){
                long long it_snp = std::distance(freq_arr.begin(), it);
                freq_arr(it_snp) = 0;
                weight_arr(it_snp + start_snp) = 0;
                snp_mat_part.col(it_snp).setZero();
                m_num_snp_used--;
                snp_index.push_back(it_snp + start_snp);
            }
        }

        Eigen::ArrayXd arr_2pq = 2 * freq_arr.array() * (1 - freq_arr.array());
        snp_mat_part.rowwise() -= arr_2pq.matrix().transpose();

        if(code_type == 1){
            m_scale += (arr_2pq * (1 - arr_2pq)).sum(); // not include the low maf SNP
        }else{
            Eigen::ArrayXd tmp_arr = weight_arr.segment(start_snp, num_snp_read);
            tmp_arr = tmp_arr / (arr_2pq * (1 - arr_2pq)).sqrt();
            weight_arr.segment(start_snp, num_snp_read) = tmp_arr;
        }
        weight_arr(snp_index).setZero(); // set weight for low maf SNP to zero, not contribute the GRM
        mat_row_elementwise_dot_vec(snp_mat_part, weight_arr.segment(start_snp, num_snp_read));  // genotypes are weighted

        if(m_start_id == 0 && m_num_id_read == m_num_id){
            mat += snp_mat_part * snp_mat_part.transpose();
        }else{
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }
    if(code_type == 2){
        m_scale = m_num_snp_used;
    }
    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);
    return mat / m_scale;
}

/**
 * @brief Computes the statistical breeding dominance relationship matrix.
 * 
 * @param npart_snp Number of SNP partitions.
 * @param weight_arr SNP weight array.
 * @param maf Minor allele frequency (MAF) threshold.
 * @param code_type Scaling type: 1 for variance-based, 2 for SNP count-based.
 * @return MatrixXd The computed GRM.
 */
MatrixXd Gmatrix::dgmatrix_gs(long long npart_snp, VectorXd& weight_arr, double maf, int code_type){
    spdlog::info("Statistical breeding dominance relationship matrix");
    if(m_input_geno_fmt != 0){
        spdlog::error("Input genotype file must be plink bed file given by --bfile");
        exit(1);
    }

    MatrixXd mat = MatrixXd::Zero(m_num_id_read, m_start_id + m_num_id_read); // G matrix
    m_scale = 0.0;
    m_num_snp_used = m_num_snp;
    GENO GENO_on_disk(m_geno_file);

    for(long long i = 0; i < npart_snp; i++){
        long long num_snp_read = (long long)(m_num_snp/npart_snp);
        long long start_snp = i*num_snp_read;
        if(i == npart_snp - 1)
            num_snp_read = num_snp_read + m_num_snp % npart_snp;
        
        spdlog::info("Processing Part {}/{}: Start SNP {}, Number of SNPs {}", i + 1, npart_snp, start_snp, num_snp_read);

        
        MatrixXd snp_mat_part; // SNP matrix stored in the memory
        VectorXd freq_arr; // SNP frequence
        VectorXd missing_rate_arr;
        vector<long long> index_vec;
        GENO_on_disk.read_geno(m_input_geno_fmt, snp_mat_part, freq_arr, missing_rate_arr, start_snp, num_snp_read, index_vec, m_missing_in_geno_vec);

        for(long long m = 0; m < m_num_id; m++){
            for(long long n = 0; n < num_snp_read; n++){
                if (snp_mat_part(m, n) < 0.5) { // 0 -> -2p^2
                    snp_mat_part(m, n) = -2 * freq_arr(n) * freq_arr(n);
                }
                else if (snp_mat_part(m, n) > 1.5) { // 2 -> -2q^2
                    snp_mat_part(m, n) = -2 * (1 - freq_arr(n)) * (1 - freq_arr(n));
                }
                else { // 1 -> 2pq
                    snp_mat_part(m, n) = 2 * (1 - freq_arr(n)) * freq_arr(n);
                }
            }
        }


        vector <long long> snp_index;
        for(auto it = freq_arr.begin(); it != freq_arr.end(); it++){
            if(*it < 0.05){
                long long it_snp = std::distance(freq_arr.begin(), it);
                freq_arr(it_snp) = 0;
                weight_arr(it_snp + start_snp) = 0;
                snp_mat_part.col(it_snp).setZero();
                m_num_snp_used--;
                snp_index.push_back(it_snp + start_snp);
            }
        }

        Eigen::ArrayXd arr_2pq = 2 * freq_arr.array() * (1 - freq_arr.array());
        if(code_type == 1){
            m_scale += arr_2pq.square().sum();
        }else{
            Eigen::ArrayXd tmp_arr = weight_arr.segment(start_snp, num_snp_read);
            tmp_arr = tmp_arr / arr_2pq;
            weight_arr.segment(start_snp, num_snp_read) = tmp_arr;
        }
        weight_arr(snp_index).setZero();
        mat_row_elementwise_dot_vec(snp_mat_part, weight_arr.segment(start_snp, num_snp_read));  // genotypes are weighted

        if(m_start_id == 0 && m_num_id_read == m_num_id){
            mat += snp_mat_part * snp_mat_part.transpose();
        }else{
            mat += snp_mat_part.middleRows(m_start_id, m_num_id_read) * snp_mat_part.topRows(m_start_id + m_num_id_read).transpose();
        }
    }
    if(code_type == 2){
        m_scale = m_num_snp_used;
    }

    spdlog::info("Number of used SNPs: {}", m_num_snp_used);
    spdlog::info("Scale factor of GRM: {}", m_scale);

    return mat / m_scale;
}



/**
 * @brief Parse command-line arguments using CLI11
 */
int Gmatrix::run(int argc, char* argv[]) {
    CLI::App app{"Gmatrix - Genomic Relationship Matrix Tool"};

    app.description(R"(
    Quick start:
    Additive GRM
        gmat --make-grm --bfile test --out test
        gmat --make-grm --bfile test --code-type 2 --out test
        gmat --make-grm --bfile test --out test --sparse

    Biological genotypic dominance GRM
        gmat --make-grm --grm-type dgrm_as --bfile test --out test
        gmat --make-grm --grm-type dgrm_as --bfile test --code-type 2 --out test

    Statistical breeding dominance GRM
        gmat --make-grm --grm-type dgrm_gs --bfile test --out test
        gmat --make-grm --grm-type dgrm_gs --bfile test --code-type 2 --out test
    )");

    bool make_grm = false;
    app.add_flag("--make-grm", make_grm, "Calculate GRM");
    
    // **Threads**
    int threads = 10;
    app.add_option("-p,--threads", threads, "Number of threads")->default_val(10);

    // **GRM Type**
    std::string grm_type = "agrm";
    app.add_option("--grm-type", grm_type, "GRM type (agrm, dgrm_as, dgrm_gs)")->default_val("agrm");

    // **File Inputs (Only One Allowed)**
    std::string bfile, dfile, cfile;
    auto opt_bfile = app.add_option("--bfile", bfile, "Prefix for PLINK binary PED files");
    auto opt_dfile = app.add_option("--dfile", dfile, "Prefix for genotype dosage files");
    auto opt_cfile = app.add_option("--cfile", cfile, "Prefix for general covariate files");
    opt_bfile->excludes(opt_dfile)->excludes(opt_cfile);
    opt_dfile->excludes(opt_bfile)->excludes(opt_cfile);
    opt_cfile->excludes(opt_bfile)->excludes(opt_dfile);

    // **Output File**
    app.add_option("--out", m_out_file, "Prefix for output file")->required();

    // **Genotype Weight Files**
    std::string awfile, dwfile;
    app.add_option("--awfile", awfile, "Genotype weights file for additive");
    app.add_option("--dwfile", dwfile, "Genotype weights file for dominance");

    // **Genotype Code Type**
    int code_type = 1;
    app.add_option("--code-type", code_type, "Genotype coded type (1 or 2)")->default_val(1);

    // **npart (Two Required Values)**
    std::vector<long long> npart_values(2, 1);
    app.add_option("--npart", npart_values, "Partition GRM elements: <npart> <ipart>")
    ->expected(2)  // **Exactly two values required**
    ->check(CLI::PositiveNumber)
    ->default_val(std::vector<long long>{1, 1});

    // **SNP Partition**
    long long npart_snp = 20;
    app.add_option("--npart-snp", npart_snp, "Partition SNPs into n parts")->default_val(20);

    // **MAF Cutoff**
    double maf = 0.00001;
    app.add_option("--maf", maf, "MAF cutoff value")->default_val(0.00001);

    // **Sparse Matrix Cutoff**
    double cut_value = 0.05;
    app.add_option("--cut-value", cut_value, "Cutoff value for sparse matrix")->default_val(0.05);

    // **Diagonal Addition Value**
    double val = 0.001;
    app.add_option("--val", val, "Add value to diagonal to ensure positive matrix")->default_val(0.001);

    // **Sparse GRM**
    bool sparse = false;
    app.add_flag("--sparse", sparse, "Enable sparse GRM output");

    // **Missing Genotype Values**
    std::vector<std::string> missing_geno = {
        "NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan",
        "<NA>", "<na>", "N/A", "n/a"
    };
    app.add_option("--missing-geno", missing_geno, "Missing genotype values (space-separated)")
        ->expected(-1)  // Accept multiple values
        ->default_val(missing_geno);

    CLI11_PARSE(app, argc, argv);

    // **Apply Settings**
    m_missing_in_geno_vec = missing_geno;

    // **Set OpenMP & MKL Threads**
    if (threads <= 0) threads = 10;
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);

    // **Set Genotype File**
    if (!bfile.empty()) {
        m_geno_file = bfile;
        m_input_geno_fmt = 0;
    } else if (!dfile.empty()) {
        m_geno_file = dfile;
        m_input_geno_fmt = 1;
    } else if (!cfile.empty()) {
        m_geno_file = cfile;
        m_input_geno_fmt = 2;
    } else {
        spdlog::error("One of --bfile, --dfile, or --cfile must be provided.");
        exit(1);
    }

    // **Assign `npart` values**
    long long npart = npart_values[0];
    long long ipart = npart_values[1];
    if(ipart <= 0 || ipart > npart){
        spdlog::error("Check the argument for --npart");
        exit(1);
    }

    // **Logging**
    spdlog::info("Threads: {}", threads);
    spdlog::info("GRM Type: {}", grm_type);
    spdlog::info("Genotype File: {}", m_geno_file);
    spdlog::info("Output File: {}", m_out_file);
    spdlog::info("Partitions: npart={}, ipart={}", npart, ipart);
    spdlog::info("SNP Partitions: {}", npart_snp);
    spdlog::info("MAF Cutoff: {}", maf);
    spdlog::info("Sparse GRM: {}", sparse ? "Enabled" : "Disabled");
    spdlog::info("Diagonal Addition Value: {}", val);
    spdlog::info("Missing Genotype Values: {}", join_string(m_missing_in_geno_vec, " "));

    // Load genotype data
    GENO GENO_on_disk(m_geno_file);
    m_num_id = GENO_on_disk.get_num_iid();
    m_num_snp = GENO_on_disk.get_num_sid();
    m_id_in_geno_vec = GENO_on_disk.iid_vec();

    // Logging with spdlog
    spdlog::info("Genotype file: {}", m_geno_file);
    spdlog::info("Output file: {}", m_out_file);
    spdlog::info("Input genotype format: {}", m_input_geno_fmt);
    spdlog::info("The number of individuals: {}", m_num_id);
    spdlog::info("The number of SNPs: {}", m_num_snp);

    this->iid_part(npart, ipart);
    VectorXd add_weight_arr, dom_weight_arr;
    MatrixXd gmat;
    // Process Additive GRM
    string out_file_name;
    if (grm_type == "agrm") {
        spdlog::info("Computing Additive Genomic Relationship Matrix");
        this->normalized_weights(add_weight_arr, awfile);
        gmat = this->agmatrix(npart_snp, add_weight_arr, maf, code_type);
        out_file_name = m_out_file + ".agrm";
    } 
    // Process Dominance GRM (association)
    else if (grm_type == "dgrm_as") {
        spdlog::info("Computing Dominance GRM (association)");
        this->normalized_weights(dom_weight_arr, dwfile);
        gmat = this->dgmatrix_as(npart_snp, dom_weight_arr, maf, code_type);
        out_file_name = m_out_file + ".dgrm_as";
    } 
    // Process Dominance GRM (Genomic Selection)
    else if (grm_type == "dgrm_gs") {
        spdlog::info("Computing Dominance GRM (Genomic Selection)");
        this->normalized_weights(dom_weight_arr, dwfile);
        gmat = this->dgmatrix_gs(npart_snp, dom_weight_arr, maf, code_type);
        out_file_name = m_out_file + ".dgrm_gs";
    } 
    // Handle invalid input
    else {
        spdlog::error("Invalid argument for --grm-fmt: {}", grm_type);
        exit(1);
    }

    // Handle partitioned output file naming
    if (npart > 1) {
        out_file_name += "." + std::to_string(npart) + "_" + std::to_string(ipart);
    }

    // Set the final output file
    this->set_out_file(out_file_name);
    spdlog::info("Output file set to: {}", out_file_name);

    this->output_iid();
    this->out_num_snp_used();
    this->out_gmat(gmat, sparse, cut_value, val);

    return 0;
}


