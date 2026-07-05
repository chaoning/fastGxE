/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:13:31
 * LastEditTime: 2026-04-01 17:47:23
 * LastEditors: Chao Ning
 */


#include <cstdint>

#define EIGEN_USE_MKL_ALL  // !must be before include Eigen

#include "geno.hpp"
#include "string_utils.hpp"
#include "EigenMatrix_utils.hpp"
#include "iterator_utils.hpp"
#include "fatal_error.hpp"

namespace {

constexpr double kMissingGenotypeValue = -9.0;
constexpr unsigned char kBedGenotypeMask = 0x03;

// Returns a .bim annotation line with the genetic-distance (cM) field removed,
// i.e. "chrom SNP cm base allele1 allele2" -> "chrom SNP base allele1 allele2".
inline std::string bim_line_without_cm(const std::string& bim_line) {
    std::vector<std::string> f = split_string(bim_line);
    if (f.size() >= 6) f.erase(f.begin() + 2);
    return join_string(f, " ");
}

inline std::int64_t bed_bytes_per_snp(std::int64_t num_id) {
    return (num_id + 3) / 4;
}

inline double decode_bed_sample(const char* snp_bytes, std::int64_t sample_idx) {
    // sample_idx>>2: byte index (4 samples/byte); (sample_idx&3)<<1: bit offset within that byte.
    const unsigned char byte = static_cast<unsigned char>(snp_bytes[sample_idx >> 2]);
    const unsigned char genotype = (byte >> ((sample_idx & 3) << 1)) & kBedGenotypeMask;
    switch (genotype) {
        case 0: return 2.0;
        case 1: return kMissingGenotypeValue;
        case 2: return 1.0;
        case 3: return 0.0;
        default: throw std::logic_error("Invalid genotype value");
    }
}

inline bool is_missing_genotype(double geno_value) {
    return geno_value < 0.0;
}

std::int64_t count_lines_or_throw(const string& file_path, const char* file_label) {
    ifstream fin(file_path);
    if (!fin.is_open()) {
        spdlog::error("Failed to open the {} file: {}", file_label, file_path);
        throw std::runtime_error("File open error");
    }

    std::int64_t line_count = 0;
    string one_line;
    while (getline(fin, one_line)) {
        ++line_count;
    }
    return line_count;
}

vector<string> read_id_column_or_throw(const string& file_path,
                                       const char* file_label,
                                       size_t expected_num_columns,
                                       size_t id_column_index,
                                       size_t reserve_count) {
    ifstream fin(file_path);
    if (!fin.is_open()) {
        spdlog::error("Failed to open the {} file: {}", file_label, file_path);
        throw std::runtime_error("File open error");
    }

    vector<string> id_vec;
    id_vec.reserve(reserve_count);

    string one_line;
    vector<string> fields;
    while (getline(fin, one_line)) {
        process_line(one_line);
        fields = split_string(one_line);
        if (fields.size() != expected_num_columns) {
            spdlog::error("Line: {} should have {} columns!", one_line, expected_num_columns);
            throw std::runtime_error("Incorrect genotype sidecar file format");
        }
        id_vec.push_back(std::move(fields[id_column_index]));
    }

    return id_vec;
}

vector<std::int64_t> locate_ids_or_throw(const vector<string>& ids_in_need_vec,
                                      const vector<string>& ids_in_geno_vec) {
    std::unordered_map<string, std::int64_t> id_in_geno_map;
    id_in_geno_map.reserve(ids_in_geno_vec.size());

    std::int64_t val = 0;
    for (const auto& tmp : ids_in_geno_vec) {
        id_in_geno_map.emplace(tmp, val++);
    }

    set<string> missing_id_set;
    for (const auto& tmp : ids_in_need_vec) {
        if (id_in_geno_map.find(tmp) == id_in_geno_map.end()) {
            missing_id_set.insert(tmp);
        }
    }
    if (!missing_id_set.empty()) {
        vector<string> missing_ids(missing_id_set.begin(), missing_id_set.end());
        fatal_error("Some requested IDs were not found: {}", join_string(missing_ids, " "));
    }

    vector<std::int64_t> index_vec;
    index_vec.reserve(ids_in_need_vec.size());
    for (const auto& tmp : ids_in_need_vec) {
        index_vec.push_back(id_in_geno_map.at(tmp));
    }
    return index_vec;
}

void parse_feature_row_or_throw(const vector<string>& fields,
                                const vector<string>& missing_in_geno_vec,
                                vector<double>& parsed_values,
                                double& observed_sum,
                                std::int64_t& num_missing_geno) {
    parsed_values.resize(fields.size());
    observed_sum = 0.0;
    num_missing_geno = 0;

    for (size_t i = 0; i < fields.size(); ++i) {
        if (is_nan(fields[i], missing_in_geno_vec)) {
            parsed_values[i] = kMissingGenotypeValue;
            ++num_missing_geno;
        } else {
            parsed_values[i] = string_to_double(fields[i]);
            observed_sum += parsed_values[i];
        }
    }
}

}  // namespace

/**
 * @brief Initialize genotype metadata from the PLINK prefix.
 *
 * The constructor records the file prefix, counts individuals from `.fam`,
 * and counts markers from `.bim`. The actual genotype matrix is still read
 * lazily by the downstream `read_*` helpers.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
GENO::GENO(const string& geno_file) : _geno_file(geno_file) {
    _num_id = count_lines_or_throw(_geno_file + ".fam", "fam");
    _num_snp = count_lines_or_throw(_geno_file + ".bim", "bim");
}

/**
 * @brief Return the number of samples recorded in the backing genotype dataset.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
std::int64_t GENO::get_num_iid() {
    return _num_id;
}

/**
 * @brief Return the number of markers recorded in the backing genotype dataset.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
std::int64_t GENO::get_num_sid() {
    return _num_snp;
}


/**
 * @brief Read the sample IDs from the `.fam` file in file order.
 *
 * The returned IDs are the second column of the PLINK `.fam` file and define
 * the row order used by genotype reads.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `fastgxe/mom.cpp`
 * - `utils/geno.cpp`
 */
vector<string> GENO::iid_vec(){
    if (!_iid_cache_loaded) {
        _iid_cache = read_id_column_or_throw(_geno_file + ".fam", "fam", 6, 1,
                                             static_cast<size_t>(_num_id));
        _iid_cache_loaded = true;
    }
    return _iid_cache;
}

/**
 * @brief Read the marker IDs from the `.bim` file in file order.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
vector<string> GENO::sid_vec(){
    if (!_sid_cache_loaded) {
        _sid_cache = read_id_column_or_throw(_geno_file + ".bim", "bim", 6, 1,
                                             static_cast<size_t>(_num_snp));
        _sid_cache_loaded = true;
    }
    return _sid_cache;
}



/**
 * @brief Read a contiguous block of raw `.bim` annotation lines.
 *
 * Each returned string is the processed `.bim` line for one SNP in the
 * requested range.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 */
vector<string> GENO::snp_anno(std::int64_t start_snp, std::int64_t num_snp_read){
    if (start_snp < 0 || num_snp_read <= 0 || start_snp + num_snp_read > _num_snp) {
        fatal_error("The start SNP or end SNP exceeds the boundary, please check!");
    }

    const string bim_file = _geno_file + ".bim";
    ifstream fin(bim_file);
    if (!fin.is_open()) {
        fatal_error("Fail to open the bim file: {}", bim_file);
    }

    string one_line;
    for (std::int64_t snp_index = 0; snp_index < start_snp; ++snp_index) {
        if (!getline(fin, one_line)) {
            fatal_error("Failed to skip to SNP index {} in {}", start_snp, bim_file);
        }
    }

    vector<string> snp_anno_vec;
    snp_anno_vec.reserve(static_cast<size_t>(num_snp_read));
    for (std::int64_t offset = 0; offset < num_snp_read; ++offset) {
        if (!getline(fin, one_line)) {
            fatal_error("Expected {} SNP annotation rows from index {} in {}",
                          num_snp_read, start_snp, bim_file);
        }
        process_line(one_line);
        snp_anno_vec.push_back(bim_line_without_cm(one_line));
    }
    return snp_anno_vec;
}


/**
 * @brief Read raw `.bim` annotation lines for an arbitrary SNP index set.
 *
 * Used in:
 * - currently no in-tree call sites were found
 */
vector<string> GENO::snp_anno_by_snp_index(const vector<std::int64_t>& snp_index_vec){
    if (snp_index_vec.empty()) {
        return {};
    }

    const string bim_file = _geno_file + ".bim";
    ifstream fin(bim_file);
    if(!fin.is_open()){
        fatal_error("Fail to open the bim file:{}", bim_file);
    }
    vector<std::pair<std::int64_t, size_t>> requested_indices;
    requested_indices.reserve(snp_index_vec.size());
    for (size_t i = 0; i < snp_index_vec.size(); ++i) {
        const std::int64_t snp_index = snp_index_vec[i];
        if (snp_index < 0 || snp_index >= _num_snp) {
            fatal_error("SNP index {} is out of range [0, {})", snp_index, _num_snp);
        }
        requested_indices.emplace_back(snp_index, i);
    }
    std::sort(requested_indices.begin(), requested_indices.end());

    vector<string> snp_anno_vec(snp_index_vec.size());
    string one_line;
    std::int64_t current_snp_index = -1;
    for (size_t request_idx = 0; request_idx < requested_indices.size(); ++request_idx) {
        const std::int64_t target_snp_index = requested_indices[request_idx].first;
        while (current_snp_index < target_snp_index) {
            if (!getline(fin, one_line)) {
                fatal_error("Failed to read SNP annotation row {} from {}", target_snp_index, bim_file);
            }
            ++current_snp_index;
        }
        process_line(one_line);
        const string anno = bim_line_without_cm(one_line);
        snp_anno_vec[requested_indices[request_idx].second] = anno;
        while (request_idx + 1 < requested_indices.size() &&
               requested_indices[request_idx + 1].first == target_snp_index) {
            ++request_idx;
            snp_anno_vec[requested_indices[request_idx].second] = anno;
        }
    }
    return snp_anno_vec;
}


/**
 * @brief Compute per-SNP allele frequency and observed-genotype counts from a PLINK BED file.
 *
 * The BED file is decoded SNP by SNP. Each 2-bit PLINK genotype code is
 * translated into the project's 0/1/2 dosage scale, while missing calls are
 * excluded from both the dosage sum and the observed-count denominator.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
void GENO::bed_allele_freq(const string& in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr){
    const std::int64_t num_byte_for_one_snp = bed_bytes_per_snp(_num_id);

    ifstream fin(in_file, std::ios::binary);
    if (!fin.is_open()) {
        fatal_error("Fail to open the plink bed file: {}", in_file);
    }

    vector<char> bytes_vec(num_byte_for_one_snp);
    freq_arr.setZero(_num_snp);
    nobs_geno_arr.setZero(_num_snp);
    if (!fin.read(bytes_vec.data(), sizeof(char) * 3)) {
        fatal_error("Failed to read the BED header from {}", in_file);
    }

    std::int64_t isnp = 0;
    while (fin.read(bytes_vec.data(), sizeof(char) * num_byte_for_one_snp)) {
        double dosage_sum = 0.0;
        std::int64_t num_missing_geno = 0;
        for (std::int64_t iid = 0; iid < _num_id; ++iid) {
            const double geno_value = decode_bed_sample(bytes_vec.data(), iid);
            if (is_missing_genotype(geno_value)) {
                ++num_missing_geno;
            } else {
                dosage_sum += geno_value;
            }
        }

        const std::int64_t nobs_geno = _num_id - num_missing_geno;
        nobs_geno_arr(isnp) = nobs_geno;
        freq_arr(isnp) = (nobs_geno == 0) ? 0.0 : dosage_sum / (2 * nobs_geno);
        ++isnp;
    }
}




/**
 * @brief Compute per-feature means and observed-value counts from a feature-by-sample text matrix.
 *
 * Each row is one feature and each column is one sample value. The matrix has
 * no header row and no leading feature-name column.
 * Missing values are excluded when computing the row mean.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
void GENO::fs_feature_mean(const string& in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr, const vector<string>& missing_in_geno_vec){
    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        fatal_error("Fail to open the genotype file:{}", in_file);
    }
    string line;
    vector<string> tmp;
    vector<double> parsed_values;
    std::int64_t igeno = 0;
    freq_arr.setZero(_num_snp);
    nobs_geno_arr.setZero(_num_snp);
    while(getline(fin, line)){
        process_line(line);
        tmp = split_string(line);
        if(tmp.size() != _num_id){
            fatal_error("The number of genotyped individuals for locus " + std::to_string(igeno + 1) + " should be " + std::to_string(_num_id));
        }

        double observed_sum = 0.0;
        std::int64_t num_missing_geno = 0;
        parse_feature_row_or_throw(tmp, missing_in_geno_vec, parsed_values, observed_sum, num_missing_geno);
        const std::int64_t nobs_geno = _num_id - num_missing_geno;
        nobs_geno_arr(igeno) = nobs_geno;
        freq_arr(igeno) = (nobs_geno == 0) ? 0.0 : observed_sum / nobs_geno;
        ++igeno;
        if(igeno > _num_snp){
            fatal_error("The number of loci exceeds " + std::to_string(_num_snp));
        }
    }
    if(igeno < _num_snp){
        fatal_error("The number of loci less than " + std::to_string(_num_snp));
    }
}


/**
 * @brief Dispatch summary-statistic computation by genotype input format.
 *
 * BED input (`0`) returns allele frequency, while feature-by-sample input (`1`)
 * returns per-feature row means in `freq_arr`.
 *
 * Used in:
 * - currently no in-tree call sites were found
 */
void GENO::allele_freq(int input_geno_fmt, VectorXd& freq_arr, VectorXd& nobs_geno_arr, const vector<string>& missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::bed_allele_freq(in_file, freq_arr, nobs_geno_arr);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".fmat";
        GENO::fs_feature_mean(in_file, freq_arr, nobs_geno_arr, missing_in_geno_vec);
    }else{
        fatal_error("The input genotypic format is not supported");
    }
}





/**
 * @brief Locate sample indices in `.fam` order with pre-validation and a hash map.
 *
 * The output preserves the order of `id_in_need_vec`.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
vector<std::int64_t> GENO::find_fam_index(const vector<string>& id_in_need_vec){
    return locate_ids_or_throw(id_in_need_vec, GENO::iid_vec());
}

/**
 * @brief Locate SNP indices in `.bim` order with pre-validation and a hash map.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
vector<std::int64_t> GENO::find_bim_index(const vector<string>& snp_in_need_vec){
    return locate_ids_or_throw(snp_in_need_vec, GENO::sid_vec());
}


/**
 * @brief Check that the `.bed` file size matches the expected sample-by-marker layout.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 */
void GENO::validate_bed_size(){
    const string bed_file = _geno_file + ".bed";
    std::ifstream fin(bed_file, std::ios::binary | std::ios::ate);
    if (!fin.is_open()) {
        fatal_error("Fail to open the plink bed file: {}", bed_file);
    }

    const std::streampos actual_file_size = fin.tellg();
    fin.close();

    const std::int64_t expected_data_bytes = (_num_id + 3) / 4 * _num_snp;
    const std::int64_t expected_file_size = expected_data_bytes + 3;
    if (actual_file_size != expected_file_size) {
        fatal_error("The size of {} doesn't correspond to the number of iids and SNPs",
                      bed_file);
    }
}

/**
 * @brief Read a contiguous BED block into a centered raw buffer.
 *
 * This variant centers genotypes on the fly and writes them into a caller-owned
 * contiguous buffer instead of an Eigen matrix.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 */
void GENO::read_bed_centered_to_buffer_omp(const string& in_file, double* snp_mat_part_pt, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr, std::int64_t start_snp, std::int64_t num_snp_read, 
                                           const vector<std::int64_t>& index_vec){
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        fatal_error("The start SNP or end SNP exceeds the boundary, please check");
    }
    const std::int64_t num_byte_for_one_snp = bed_bytes_per_snp(_num_id); // One byte stores four samples for one SNP.
    const std::int64_t total_bytes_to_read = num_byte_for_one_snp * num_snp_read;

    FILE* fin = fopen(in_file.c_str(), "rb");
    if (!fin) {
        fatal_error("Fail to open the plink bed file: " + in_file);
    }

    // Skip the 3-byte BED header, then jump to the first requested SNP block.
    const std::int64_t start_pos = start_snp * num_byte_for_one_snp + 3;
    if (fseek(fin, start_pos, SEEK_SET) != 0){
        fclose(fin);
        fatal_error("Fail to seek the predefined position");
    }

    vector<char> bytes_vec(total_bytes_to_read);
    const size_t read_count = fread(bytes_vec.data(), sizeof(char), total_bytes_to_read, fin);
    if (read_count != static_cast<size_t>(total_bytes_to_read)){
        fclose(fin);
        fatal_error("fread failed to read the expected amount");
    }
    if(ferror(fin)){
        fclose(fin);
        fatal_error("Failed to read from file: " + in_file);
    }
    fclose(fin);

    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    nobs_geno_arr.setZero(num_snp_read);

    const bool use_all_ids = index_vec.empty();
    const std::int64_t num_used_id = use_all_ids ? _num_id : static_cast<std::int64_t>(index_vec.size());

    #pragma omp parallel for schedule(dynamic)
    for(std::int64_t isnp = 0; isnp < num_snp_read; isnp++){
        const char* snp_bytes = bytes_vec.data() + isnp * num_byte_for_one_snp;
        double* output_row = snp_mat_part_pt + isnp * num_used_id;
        double freq = 0.0; // Sum of observed allele dosages before scaling to MAF.
        std::int64_t num_missing = 0;
        double code_val;

        if(use_all_ids){
            // Fast path: decode directly into the caller-provided row buffer,
            // then center that row in place once the SNP mean is known.
            for(std::int64_t i = 0; i < _num_id; i++){
                code_val = decode_bed_sample(snp_bytes, i);
                if(is_missing_genotype(code_val)){
                    num_missing++;
                }else{
                    freq += code_val;
                }
                output_row[i] = code_val;
            }
        }else{
            // Subset path: estimate statistics from the requested individuals
            // only, so the returned MAF/missingness/nobs all match the output.
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes, index_vec[i]);
                if(is_missing_genotype(code_val)){
                    num_missing++;
                }else{
                    freq += code_val;
                }
            }
        }

        const std::int64_t nobs_geno = num_used_id - num_missing;
        missing_rate_arr(isnp) = num_missing * 1.0 / num_used_id;
        nobs_geno_arr(isnp) = nobs_geno;

        // The allele frequency only uses observed genotypes.
        if(num_missing != num_used_id){
            freq /= (2 * num_used_id - 2 * num_missing);
        }else{
            freq = 0.0;
        }
        maf_arr(isnp) = freq;

        const double center = 2 * freq;
        if(use_all_ids){
            for(std::int64_t i = 0; i < _num_id; i++){
                if(output_row[i] < 0.0){
                    output_row[i] = 0.0;
                }else{
                    output_row[i] -= center;
                }
            }
        }else{
            // Only materialize the requested samples, preserving input order.
            // Missing values are mapped to 0 after centering, matching the
            // current downstream logic.
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes, index_vec[i]);
                if(is_missing_genotype(code_val)){
                    output_row[i] = 0.0;
                }else{
                    output_row[i] = code_val - center;
                }
            }
        }
    }
}

/**
 * @brief Read a contiguous BED block in parallel into a dense Eigen matrix.
 *
 * This is the main BED reader used by current GRM-building paths.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
void GENO::read_bed_omp(const string& in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr, std::int64_t start_snp, std::int64_t num_snp_read, 
                const vector<std::int64_t>& index_vec){
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        fatal_error("The start SNP or end SNP exceeds the boundary, please check");
    }
    const std::int64_t num_byte_for_one_snp = bed_bytes_per_snp(_num_id);
    const std::int64_t total_bytes_to_read = num_byte_for_one_snp * num_snp_read;
    ifstream fin;
    fin.open(in_file, std::ios::binary);
    if (!fin.is_open()) {
        fatal_error("Fail to open the plink bed file: " + in_file);
    }

    vector<char> bytes_vec(total_bytes_to_read);
    fin.read(bytes_vec.data(), sizeof(char) * 3);
    fin.seekg(start_snp * num_byte_for_one_snp, std::ios::cur);
    fin.read(bytes_vec.data(), sizeof(char) * total_bytes_to_read);

    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    nobs_geno_arr.setZero(num_snp_read);

    const bool use_all_ids = index_vec.empty();
    const std::int64_t num_used_id = use_all_ids ? _num_id : static_cast<std::int64_t>(index_vec.size());
    snp_mat_part.resize(num_used_id, num_snp_read);
    double* snp_mat_part_pt = snp_mat_part.data();

    #pragma omp parallel for schedule(static)
    for(std::int64_t isnp = 0; isnp < num_snp_read; isnp++){
        const char* snp_bytes = bytes_vec.data() + isnp * num_byte_for_one_snp;
        double* output_col = snp_mat_part_pt + isnp * num_used_id;
        double freq = 0.0;
        std::int64_t num_missing_geno = 0;
        double code_val;

        if(use_all_ids){
            // Fast path: decode directly into the output column, then replace
            // missing values with the SNP mean once frequency is known.
            for(std::int64_t i = 0; i < _num_id; i++){
                code_val = decode_bed_sample(snp_bytes, i);
                if(is_missing_genotype(code_val)){
                    num_missing_geno++;
                }else{
                    freq += code_val;
                }
                output_col[i] = code_val;
            }
        }else{
            // Subset path: estimate statistics from the requested individuals
            // only, so the returned summaries match the output matrix.
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes, index_vec[i]);
                if(is_missing_genotype(code_val)){
                    num_missing_geno++;
                }else{
                    freq += code_val;
                }
            }
        }

        nobs_geno_arr(isnp) = num_used_id - num_missing_geno;
        missing_rate_arr(isnp) = num_missing_geno * 1.0 / num_used_id;

        if(std::fabs(missing_rate_arr(isnp)) >= 0.99999){
            freq = 0.0;
        }else{
            freq /= (2 * num_used_id - 2 * num_missing_geno);
        }
        maf_arr(isnp) = freq;

        const double imputed_value = 2 * freq;
        if(use_all_ids){
            for(std::int64_t i = 0; i < _num_id; i++){
                if(is_missing_genotype(output_col[i])){
                    output_col[i] = imputed_value;
                }
            }
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes, index_vec[i]);
                output_col[i] = is_missing_genotype(code_val) ? imputed_value : code_val;
            }
        }
    }

    fin.close();
}



/**
 * @brief Read a contiguous block from a feature-by-sample text matrix into a dense matrix.
 *
 * Each row stores one feature and each column stores one sample value, with no
 * header row and no leading feature-name column. Missing values are mean-
 * imputed feature by feature. When `index_vec` is provided, rows are reordered
 * and subset to the requested sample order, and summary statistics are computed
 * on that same sample subset.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
void GENO::read_fs(const string& in_file, MatrixXd& geno_mat_part, VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
            std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec){
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        fatal_error("The start SNP or end SNP exceeds the boundary, please check");
    }

    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        fatal_error("Fail to open the genotypic file: " + in_file);
    }

    string line;
    // Skip lines before start SNP
    for(std::int64_t i = 0; i < start_snp; i++){
        getline(fin, line);
    }

    const bool use_all_ids = index_vec.empty();
    const std::int64_t num_used_id = use_all_ids ? _num_id : static_cast<std::int64_t>(index_vec.size());
    freq_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    nobs_geno_arr.setZero(num_snp_read);
    geno_mat_part.resize(num_used_id, num_snp_read);

    std::int64_t igeno = 0;
    vector<string> tmp;
    vector<double> parsed_values;
    vector<std::int64_t> missing_geno_index;
    while(igeno < num_snp_read && getline(fin, line)){
        process_line(line);
        tmp = split_string(line);
        if(tmp.size() != _num_id){
            fatal_error("The number of genotyped individuals for locus " + std::to_string(igeno + 1) + " should be " + std::to_string(_num_id));
        }

        missing_geno_index.clear();
        double* output_col = geno_mat_part.data() + igeno * num_used_id;
        double observed_sum = 0.0;
        std::int64_t num_missing_geno = 0;
        parse_feature_row_or_throw(tmp, missing_in_geno_vec, parsed_values, observed_sum, num_missing_geno);

        double selected_sum = 0.0;
        std::int64_t selected_missing = 0;
        if(use_all_ids){
            selected_sum = observed_sum;
            selected_missing = num_missing_geno;
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                const double geno_value = parsed_values[index_vec[i]];
                if(is_missing_genotype(geno_value)){
                    selected_missing++;
                }else{
                    selected_sum += geno_value;
                }
            }
        }

        const std::int64_t nobs_geno = num_used_id - selected_missing;
        const double feature_mean = (nobs_geno == 0) ? 0.0 : selected_sum / nobs_geno;

        if(use_all_ids){
            // Fast path: decode directly into the output matrix and remember
            // which entries need to be mean-imputed after the row mean is known.
            for(std::int64_t i = 0; i < _num_id; i++){
                if(is_missing_genotype(parsed_values[i])){
                    missing_geno_index.push_back(i);
                }else{
                    output_col[i] = parsed_values[i];
                }
            }
        }

        if(use_all_ids){
            for(std::int64_t idx : missing_geno_index){
                output_col[idx] = feature_mean;
            }
            freq_arr(igeno) = feature_mean;
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                const double geno_value = parsed_values[index_vec[i]];
                output_col[i] = is_missing_genotype(geno_value) ? feature_mean : geno_value;
            }
            freq_arr(igeno) = feature_mean;
        }

        nobs_geno_arr(igeno) = nobs_geno;
        missing_rate_arr(igeno) = static_cast<double>(selected_missing) / num_used_id;
        igeno++;
    }

    if(igeno < num_snp_read){
        geno_mat_part.conservativeResize(num_used_id, igeno);
        freq_arr.conservativeResize(igeno);
        missing_rate_arr.conservativeResize(igeno);
        nobs_geno_arr.conservativeResize(igeno);
    }
    fin.close();
}


/**
 * @brief Dispatch contiguous genotype-block reading by input format.
 *
 * Supported formats are BED (`0`) and feature-by-sample text input (`1`).
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 */
void GENO::read_geno(int input_geno_fmt, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
            std::int64_t start_snp, std::int64_t num_snp_read, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::read_bed_omp(in_file, snp_mat_part, maf_arr, missing_rate_arr, nobs_geno_arr, start_snp, num_snp_read, index_vec);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".fmat"; // Feature-by-sample matrix, e.g. SNPs or expression.
        GENO::read_fs(in_file, snp_mat_part, maf_arr, missing_rate_arr, nobs_geno_arr, start_snp, num_snp_read, index_vec, missing_in_geno_vec);
    }else{
        fatal_error("The input genotypic format is not supported");
    }
}


/**
 * @brief Read an arbitrary SNP subset from a BED file.
 *
 * The output columns follow `snp_index_vec` exactly, so the caller can request
 * SNPs in any order. Missing genotypes are imputed SNP by SNP using the allele
 * mean estimated from the same individuals that are returned in the output.
 * When `index_vec` is provided, allele frequency, missing rate, and observed
 * count are all computed on that requested sample subset.
 *
 * Used in:
 * - `utils/geno.cpp`
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
void GENO::read_bed_by_snp_indices(const string& in_file, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr, 
                const vector<std::int64_t>& snp_index_vec, const vector<std::int64_t>& index_vec){
    const std::int64_t num_byte_for_one_snp = bed_bytes_per_snp(_num_id);
    FILE* fin = fopen(in_file.c_str(), "rb");
    if (!fin) {
        fatal_error("Fail to open the plink bed file: " + in_file);
    }
    const bool use_all_ids = index_vec.empty();
    const std::int64_t num_used_id = use_all_ids ? _num_id : static_cast<std::int64_t>(index_vec.size());
    const std::int64_t num_snp_read = static_cast<std::int64_t>(snp_index_vec.size());
    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    nobs_geno_arr.setZero(num_snp_read);
    snp_mat_by_snp_index.resize(num_used_id, num_snp_read);
    vector<char> snp_bytes(num_byte_for_one_snp);

    for(std::int64_t isnp = 0; isnp < num_snp_read; isnp++) {
        const std::int64_t snp_idx = snp_index_vec[isnp];
        if(snp_idx < 0 || snp_idx >= _num_snp){
            fclose(fin);
            fatal_error("The requested SNP index {} exceeds the BED boundary", snp_idx);
        }

        if(fseek(fin, snp_idx * num_byte_for_one_snp + 3, SEEK_SET) != 0){
            fclose(fin);
            fatal_error("Fail to seek the predefined position");
        }
        const size_t read_count = fread(snp_bytes.data(), sizeof(char), num_byte_for_one_snp, fin);
        if(read_count != static_cast<size_t>(num_byte_for_one_snp) || ferror(fin)){
            fclose(fin);
            fatal_error("Failed to read the requested SNP block from: " + in_file);
        }

        double* output_col = snp_mat_by_snp_index.data() + isnp * num_used_id;
        double dosage_sum = 0.0;
        std::int64_t num_missing_geno = 0;
        double code_val;

        if(use_all_ids){
            for(std::int64_t i = 0; i < _num_id; i++){
                code_val = decode_bed_sample(snp_bytes.data(), i);
                if(is_missing_genotype(code_val)){
                    num_missing_geno++;
                }else{
                    dosage_sum += code_val;
                }
                output_col[i] = code_val;
            }
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes.data(), index_vec[i]);
                if(is_missing_genotype(code_val)){
                    num_missing_geno++;
                }else{
                    dosage_sum += code_val;
                }
            }
        }

        const double missing_rate = static_cast<double>(num_missing_geno) / num_used_id;
        missing_rate_arr(isnp) = missing_rate;
        nobs_geno_arr(isnp) = num_used_id - num_missing_geno;
        const double freq = (num_missing_geno == num_used_id) ? 0.0 : dosage_sum / (2 * (num_used_id - num_missing_geno));

        const double imputed_value = 2 * freq;
        if(use_all_ids){
            for(std::int64_t i = 0; i < _num_id; i++){
                if(is_missing_genotype(output_col[i])){
                    output_col[i] = imputed_value;
                }
            }
            maf_arr(isnp) = freq;
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                code_val = decode_bed_sample(snp_bytes.data(), index_vec[i]);
                output_col[i] = is_missing_genotype(code_val) ? imputed_value : code_val;
            }
            maf_arr(isnp) = freq;
        }
    }
    fclose(fin);
}


/**
 * @brief Read an arbitrary feature subset from a feature-by-sample text matrix.
 *
 * The matrix has no header row and no leading feature-name column. Output
 * columns follow `feature_index_vec` exactly, including repeated feature
 * indices. Missing values are imputed feature by feature using the row mean.
 * When `index_vec` is provided, the row mean, missing rate, and observed count
 * are all computed on the requested sample subset.
 *
 * Used in:
 * - `utils/geno.cpp`
 */
void GENO::read_fs_by_feature_indices(const string& in_file, MatrixXd& geno_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr, 
                const vector<std::int64_t>& feature_index_vec, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec){
    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        fatal_error("Fail to open the dosage genotype file: " + in_file);
    }

    string line;
    const bool use_all_ids = index_vec.empty();
    const std::int64_t num_used_id = use_all_ids ? _num_id : static_cast<std::int64_t>(index_vec.size());
    const std::int64_t num_geno_read = static_cast<std::int64_t>(feature_index_vec.size());
    maf_arr.setZero(num_geno_read);
    missing_rate_arr.setZero(num_geno_read);
    nobs_geno_arr.setZero(num_geno_read);
    geno_mat_by_snp_index.setZero(num_used_id, num_geno_read);

    std::unordered_map<std::int64_t, vector<std::int64_t>> feature_positions;
    feature_positions.reserve(feature_index_vec.size());
    for(std::int64_t pos = 0; pos < num_geno_read; pos++){
        feature_positions[feature_index_vec[pos]].push_back(pos);
    }

    vector<string> tmp;
    vector<double> parsed_values;
    vector<std::int64_t> missing_geno_index;
    vector<double> row_values(use_all_ids ? _num_id : 0);
    vector<double> selected_values(use_all_ids ? 0 : num_used_id);
    std::int64_t feature_idx = 0;
    std::int64_t num_features_found = 0;
    while(getline(fin, line)){
        auto pos_it = feature_positions.find(feature_idx);
        if(pos_it == feature_positions.end()){
            feature_idx++;
            continue;
        }

        process_line(line);
        tmp = split_string(line);
        if(tmp.size() != _num_id){
            fatal_error("The number of genotyped individuals for locus " + std::to_string(feature_idx + 1) + " should be " + std::to_string(_num_id));
        }

        missing_geno_index.clear();
        double observed_sum = 0.0;
        std::int64_t num_missing_geno = 0;
        parse_feature_row_or_throw(tmp, missing_in_geno_vec, parsed_values, observed_sum, num_missing_geno);

        double selected_sum = 0.0;
        std::int64_t selected_missing = 0;
        if(use_all_ids){
            selected_sum = observed_sum;
            selected_missing = num_missing_geno;
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                const double geno_value = parsed_values[index_vec[i]];
                if(is_missing_genotype(geno_value)){
                    selected_missing++;
                }else{
                    selected_sum += geno_value;
                }
            }
        }

        const std::int64_t nobs_geno = num_used_id - selected_missing;
        const double feature_mean = (nobs_geno == 0) ? 0.0 : selected_sum / nobs_geno;

        if(use_all_ids){
            for(std::int64_t i = 0; i < _num_id; i++){
                if(is_missing_genotype(parsed_values[i])){
                    missing_geno_index.push_back(i);
                }else{
                    row_values[i] = parsed_values[i];
                }
            }
        }

        const double missing_rate = static_cast<double>(selected_missing) / num_used_id;
        const vector<std::int64_t>& output_positions = pos_it->second;

        if(use_all_ids){
            for(std::int64_t idx : missing_geno_index){
                row_values[idx] = feature_mean;
            }
            for(std::int64_t output_pos : output_positions){
                double* output_col = geno_mat_by_snp_index.data() + output_pos * num_used_id;
                std::copy(row_values.begin(), row_values.end(), output_col);
                maf_arr(output_pos) = feature_mean;
                missing_rate_arr(output_pos) = missing_rate;
                nobs_geno_arr(output_pos) = nobs_geno;
            }
        }else{
            for(std::int64_t i = 0; i < num_used_id; i++){
                const double geno_value = parsed_values[index_vec[i]];
                selected_values[i] = is_missing_genotype(geno_value) ? feature_mean : geno_value;
            }
            for(std::int64_t output_pos : output_positions){
                double* output_col = geno_mat_by_snp_index.data() + output_pos * num_used_id;
                std::copy(selected_values.begin(), selected_values.end(), output_col);
                maf_arr(output_pos) = feature_mean;
                missing_rate_arr(output_pos) = missing_rate;
                nobs_geno_arr(output_pos) = nobs_geno;
            }
        }

        num_features_found += output_positions.size();
        if(num_features_found >= num_geno_read){
            break;
        }
        feature_idx++;
    }
    fin.close();
}


/**
 * @brief Dispatch arbitrary-index genotype reading by input format.
 *
 * Used in:
 * - currently no in-tree call sites were found
 */
void GENO::read_geno_by_index(int input_geno_fmt, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr,
                const vector<std::int64_t>& snp_index_vec, const vector<std::int64_t>& index_vec, const vector<string>& missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::read_bed_by_snp_indices(in_file, snp_mat_by_snp_index, maf_arr, missing_rate_arr, nobs_geno_arr,
                snp_index_vec, index_vec);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".fmat";
        GENO::read_fs_by_feature_indices(in_file, snp_mat_by_snp_index, maf_arr, missing_rate_arr, nobs_geno_arr,
                snp_index_vec, index_vec, missing_in_geno_vec);
    }else{
        fatal_error("The input genotypic format is not supported");
    }
}
