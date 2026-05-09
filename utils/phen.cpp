#include <cstdint>

#define EIGEN_USE_MKL_ALL  // must be before include Eigen

#include <algorithm>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <spdlog/spdlog.h>

#include "omp.h"

#include "EigenMatrix_utils.hpp"
#include "iterator_utils.hpp"
#include "phen.hpp"
#include "string_utils.hpp"
#include "fatal_error.hpp"

using std::ifstream;
using std::vector;
using Eigen::HouseholderQR;

namespace {

size_t count_remaining_lines(ifstream& fin) {
    const std::streampos data_start = fin.tellg();
    if (data_start == std::streampos(-1)) {
        return 0;
    }

    size_t num_rows = 0;
    string line;
    while (getline(fin, line)) {
        ++num_rows;
    }

    fin.clear();
    fin.seekg(data_start);
    return num_rows;
}

const vector<string>& get_checked_column_or_exit(const vector<vector<string>>& data_vec,
                                                 std::int64_t column_no) {
    if (column_no < 0 || static_cast<size_t>(column_no) >= data_vec.size()) {
        fatal_error("Phenotype column index {} is out of range for data with {} columns",
                      column_no, data_vec.size());
    }

    return data_vec[static_cast<size_t>(column_no)];
}

const string& get_checked_header_name_or_exit(const vector<string>& header_vec,
                                              std::int64_t column_no) {
    if (column_no < 0 || static_cast<size_t>(column_no) >= header_vec.size()) {
        fatal_error("Phenotype header index {} is out of range for {} stored headers",
                      column_no, header_vec.size());
    }

    return header_vec[static_cast<size_t>(column_no)];
}

}  // namespace

/**
 * @brief Read the full phenotype table into a column-major string buffer.
 *
 * The input file is expected to contain a header row followed by whitespace-
 * separated records. An extra column is appended to `data_vec` to store the
 * original 0-based row number for each record.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
void PHEN::read_all_columns_raw(const string& pheno_file, vector<vector<string>>& data_vec){
    ifstream fin(pheno_file);
    if (!fin) {
        fatal_error("Fail to open the phenotypic file");
    }

    string line;

    // ===================== Read header =====================
    getline(fin, line);
    m_head_column_vec = split_string(line);
    const size_t num_col = m_head_column_vec.size();

    // +1 for row index column
    data_vec.assign(num_col + 1, {});
    const size_t num_rows = count_remaining_lines(fin);
    for (auto& column : data_vec) {
        column.reserve(num_rows);
    }

    size_t row_id = 0;

    // ===================== Read data =====================
    while (getline(fin, line)) {
        process_line(line);

        auto fields = split_string(line);

        if (fields.size() != num_col) {
            fatal_error("Column mismatch at row {}: {}", row_id, line);
        }

        // Move strings directly into column storage (no copy)
        for (size_t i = 0; i < num_col; ++i)
            data_vec[i].push_back(std::move(fields[i]));

        data_vec[num_col].push_back(std::to_string(row_id));
        ++row_id;
    }
}

/**
 * @brief Read only selected phenotype columns plus the original row index.
 *
 * This is a lighter-weight alternative to `read_all_columns_raw` when the caller only
 * needs a subset of columns.
 *
 * Used in:
 * - currently no in-tree call sites were found
 */
void PHEN::read_selected_columns_raw(const string& pheno_file,
                                     const vector<std::int64_t>& target_cols,
                                     vector<vector<string>>& data_vec) {

    // Open phenotype file
    ifstream fin(pheno_file);
    if (!fin) {
        fatal_error("Fail to open the phenotypic file");
    }

    string line;

    // ===================== Read header =====================
    getline(fin, line);
    m_head_column_vec = split_string(line);
    const size_t num_col = m_head_column_vec.size();

    // Validate requested column indices once and convert to storage indices.
    vector<size_t> selected_indices;
    selected_indices.reserve(target_cols.size());
    for (auto c : target_cols) {
        if (c < 0 || static_cast<size_t>(c) >= num_col) {
            fatal_error("Requested column {} out of range (0..{})", c, num_col - 1);
        }
        selected_indices.push_back(static_cast<size_t>(c));
    }

    const size_t K = selected_indices.size();

    // Allocate output: K selected columns + 1 column for row index
    data_vec.assign(K + 1, {});
    const size_t num_rows = count_remaining_lines(fin);
    for (auto& column : data_vec) {
        column.reserve(num_rows);
    }

    size_t row_id = 0;

    // ===================== Read data =====================
    while (getline(fin, line)) {
        process_line(line);

        auto fields = split_string(line);

        // Check column consistency
        if (fields.size() != num_col) {
            fatal_error("Column mismatch at row {}: {}", row_id, line);
        }

        // Move only the requested columns into the output buffer.
        for (size_t j = 0; j < K; ++j) {
            const size_t col = selected_indices[j];
            data_vec[j].push_back(std::move(fields[col]));
        }

        // Store row index
        data_vec[K].push_back(std::to_string(row_id));
        ++row_id;
    }
}

/**
 * @brief Keep only records whose IDs appear in every requested GRM ID list.
 *
 * `random_grm_vec` stores the phenotype-column indices that hold the sample ID
 * for each random effect, and `grm_id_used_vec_vec` stores the matched GRM IDs.
 * The result is accumulated into `m_index_keep_vec` using the project's 0-based
 * row-index convention.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
void PHEN::update_keep_indices_by_grm_ids(const vector<std::int64_t>& random_grm_vec,
                                          const vector<vector<string>>& grm_id_used_vec_vec){
    if (random_grm_vec.size() != grm_id_used_vec_vec.size()) {
        fatal_error("random_grm_vec size ({}) does not match grm_id_used_vec_vec size ({})",
                      random_grm_vec.size(), grm_id_used_vec_vec.size());
    }

    for (size_t i = 0; i < random_grm_vec.size(); ++i) {
        const std::int64_t column_index = random_grm_vec[i];
        if (column_index < 0 || static_cast<size_t>(column_index) >= m_data_vec.size()) {
            fatal_error("GRM ID column index {} is out of range for phenotype data with {} columns",
                          column_index, m_data_vec.size());
        }

        const std::unordered_set<string> grm_id_set(grm_id_used_vec_vec[i].begin(),
                                                    grm_id_used_vec_vec[i].end());
        const auto& data_id_vec = m_data_vec[static_cast<size_t>(column_index)];

        vector<std::int64_t> tmp_index_keep_vec;
        tmp_index_keep_vec.reserve(data_id_vec.size());

        for (size_t row = 0; row < data_id_vec.size(); ++row) {
            if (grm_id_set.find(data_id_vec[row]) != grm_id_set.end()) {
                tmp_index_keep_vec.push_back(static_cast<std::int64_t>(row));
            }
        }

        m_index_keep_vec = set_intersection_(m_index_keep_vec, tmp_index_keep_vec);
        if (m_index_keep_vec.empty()) {
            break;
        }
    }
}

/**
 * @brief Remove records with missing values in any of the requested columns.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
void PHEN::update_keep_indices_by_non_missing_values(const vector<std::int64_t>& col_used_vec,
                                                     const vector<string>& missing_in_data_vec){
    if (m_index_keep_vec.empty() || col_used_vec.empty() || missing_in_data_vec.empty()) {
        return;
    }

    const std::unordered_set<string> missing_tags(missing_in_data_vec.begin(),
                                                  missing_in_data_vec.end());

    vector<std::int64_t> next_keep_indices;
    next_keep_indices.reserve(m_index_keep_vec.size());

    for (std::int64_t col : col_used_vec) {
        if (col < 0 || static_cast<size_t>(col) >= m_data_vec.size()) {
            fatal_error("Phenotype column index {} is out of range for data with {} columns",
                          col, m_data_vec.size());
        }

        const auto& column = m_data_vec[static_cast<size_t>(col)];
        next_keep_indices.clear();

        for (std::int64_t row_index_zero_based : m_index_keep_vec) {
            const size_t row_index = static_cast<size_t>(row_index_zero_based);
            if (missing_tags.find(column[row_index]) == missing_tags.end()) {
                next_keep_indices.push_back(row_index_zero_based);
            }
        }

        m_index_keep_vec.swap(next_keep_indices);
        if (m_index_keep_vec.empty()) {
            break;
        }
    }
}

/**
 * @brief Materialize the filtered row set into `m_data_vec`.
 *
 * After `m_index_keep_vec` is finalized, this helper rewrites every stored
 * column so that only retained rows remain, preserving their current order.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
void PHEN::apply_keep_indices_to_data(){
    if (m_data_vec.empty()) {
        return;
    }

    const size_t num_rows = m_data_vec[0].size();
    const size_t kept_count = m_index_keep_vec.size();

    if (kept_count == 0) {
        for (auto& column : m_data_vec) {
            column.clear();
        }
        return;
    }

    vector<size_t> keep_indices_zero_based(kept_count);
    bool keeps_all_rows_in_order = kept_count == num_rows;

    for (size_t i = 0; i < kept_count; ++i) {
        const std::int64_t row_index_zero_based = m_index_keep_vec[i];
        if (row_index_zero_based < 0 || static_cast<size_t>(row_index_zero_based) >= num_rows) {
            fatal_error("Keep index {} is out of range for phenotype data with {} rows",
                          row_index_zero_based, num_rows);
        }

        keep_indices_zero_based[i] = static_cast<size_t>(row_index_zero_based);
        if (keeps_all_rows_in_order && keep_indices_zero_based[i] != i) {
            keeps_all_rows_in_order = false;
        }
    }

    if (keeps_all_rows_in_order) {
        return;
    }

    #pragma omp parallel for schedule(static)
    for (size_t column_index = 0; column_index < m_data_vec.size(); ++column_index) {
        auto& column = m_data_vec[column_index];

        vector<string> filtered_column;
        filtered_column.reserve(kept_count);

        for (size_t row_index : keep_indices_zero_based) {
            filtered_column.push_back(std::move(column[row_index]));
        }

        column.swap(filtered_column);
    }
}

/**
 * @brief Read phenotype data and apply record-level filtering.
 * 
 * @param pheno_file phenotypic file
 * @param random_grm_vec random effects linked with grm files
 * @param grm_id_used_vec_vec iids vector for grm files
 * @param col_used_vec used columns
 * @param missing_in_data_vec missing data
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
void PHEN::read_and_filter(const string& pheno_file,
                           const vector<std::int64_t>& random_grm_vec,
                           const vector<vector<string>>& grm_id_used_vec_vec,
                           const vector<std::int64_t>& col_used_vec,
                           const vector<string>& missing_in_data_vec){
    // Read all phenotype data into the class-owned storage.
    m_data_vec.clear();
    PHEN::read_all_columns_raw(pheno_file, m_data_vec);

    if (m_data_vec.empty()) {
        fatal_error("No phenotype columns were loaded from {}", pheno_file);
    }

    const size_t num_raw_records = m_data_vec[0].size();
    if (num_raw_records == 0) {
        fatal_error("No phenotype records were found in {}", pheno_file);
    }

    spdlog::info("Loaded {} phenotype records", num_raw_records);

    m_index_keep_vec.resize(num_raw_records);
    std::iota(m_index_keep_vec.begin(), m_index_keep_vec.end(), 0LL);

    PHEN::update_keep_indices_by_grm_ids(random_grm_vec, grm_id_used_vec_vec);
    spdlog::info("{} records remain after GRM-id filtering", m_index_keep_vec.size());
    
    PHEN::update_keep_indices_by_non_missing_values(col_used_vec, missing_in_data_vec);
    spdlog::info("{} records remain after non-missing filtering", m_index_keep_vec.size());
    
    const size_t num_records_kept = m_index_keep_vec.size();
    const double ratio = static_cast<double>(num_records_kept) /
                         static_cast<double>(num_raw_records);
    if (num_records_kept == 0) {
        fatal_error("No records remain after phenotype filtering.");
    }
    if (ratio < 0.05) {
        spdlog::warn("{} records remain after phenotype filtering ({:.2f}% of {}).",
                     num_records_kept, ratio * 100.0, num_raw_records);
    }
    
    PHEN::apply_keep_indices_to_data();
}

/**
 * @brief Reorder all retained phenotype records by a caller-provided sample-ID order.
 *
 * `sample_id_order_vec` must be a permutation of the current sample-ID column
 * (column 0). Every stored column, including the original-row bookkeeping
 * column, is rewritten in that order.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 */
void PHEN::reorder_records_by_sample_id_order(const vector<string>& sample_id_order_vec){
    if (m_data_vec.empty()) {
        if (!sample_id_order_vec.empty()) {
            fatal_error("Cannot reorder phenotype records: no phenotype data are loaded.");
        }
        return;
    }

    const auto& current_sample_ids = m_data_vec[0];
    const size_t num_records = current_sample_ids.size();
    if (sample_id_order_vec.size() != num_records) {
        fatal_error("Sample-ID order size ({}) does not match the number of retained records ({})",
                      sample_id_order_vec.size(), num_records);
    }

    std::unordered_map<string, size_t> sample_id_to_index;
    sample_id_to_index.reserve(num_records);
    for (size_t row = 0; row < num_records; ++row) {
        const auto [it, inserted] = sample_id_to_index.emplace(current_sample_ids[row], row);
        if (!inserted) {
            fatal_error("Duplicate sample ID '{}' found in phenotype data while reordering.",
                          current_sample_ids[row]);
        }
    }

    vector<size_t> reordered_indices;
    reordered_indices.reserve(num_records);
    vector<char> used_index(num_records, 0);
    bool keeps_current_order = true;

    for (size_t target_row = 0; target_row < sample_id_order_vec.size(); ++target_row) {
        const auto& sample_id = sample_id_order_vec[target_row];
        const auto it = sample_id_to_index.find(sample_id);
        if (it == sample_id_to_index.end()) {
            fatal_error("Sample ID '{}' was requested for reordering but is absent from phenotype data.",
                          sample_id);
        }

        const size_t source_row = it->second;
        if (used_index[source_row]) {
            fatal_error("Sample ID '{}' appears multiple times in the requested sample-ID order.",
                          sample_id);
        }

        used_index[source_row] = 1;
        reordered_indices.push_back(source_row);
        if (keeps_current_order && source_row != target_row) {
            keeps_current_order = false;
        }
    }

    if (keeps_current_order) {
        return;
    }

    #pragma omp parallel for schedule(static)
    for (size_t column_index = 0; column_index < m_data_vec.size(); ++column_index) {
        auto& column = m_data_vec[column_index];

        vector<string> reordered_column;
        reordered_column.reserve(num_records);

        for (size_t source_row : reordered_indices) {
            reordered_column.push_back(std::move(column[source_row]));
        }

        column.swap(reordered_column);
    }
}

/**
 * @brief Copy stored column names, including the appended row index column when present.
 *
 * Used in:
 * - `test/phen_workflow_example.cpp`
 */
void PHEN::get_stored_column_names(vector<string>& column_names) const{
    column_names = m_head_column_vec;
    if (!m_data_vec.empty() && m_data_vec.size() == m_head_column_vec.size() + 1) {
        column_names.push_back("row_id");
    }
}

/**
 * @brief Copy one stored column as strings.
 * 
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
void PHEN::get_given_column(std::int64_t column_index, vector<string>& values) const{
    values = get_checked_column_or_exit(m_data_vec, column_index);
}

/**
 * @brief Copy one stored column into an Eigen vector of doubles.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
void PHEN::get_given_column(std::int64_t column_index, VectorXd& values) const{
    const auto& source_column = get_checked_column_or_exit(m_data_vec, column_index);
    values.resize(static_cast<Eigen::Index>(source_column.size()));

    for (Eigen::Index row_index = 0; row_index < values.size(); ++row_index) {
        values(row_index) = string_to_double(source_column[static_cast<size_t>(row_index)]);
    }
}

/**
 * @brief Gather several phenotype columns into a dense numeric matrix.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
void PHEN::get_columns_as_matrix(const vector<std::int64_t>& column_indices, MatrixXd& matrix) const{
    if (m_data_vec.empty()) {
        fatal_error("Cannot access phenotype columns: no phenotype data are loaded.");
    }

    const Eigen::Index num_rows = static_cast<Eigen::Index>(m_data_vec[0].size());
    const Eigen::Index num_columns = static_cast<Eigen::Index>(column_indices.size());
    matrix.resize(num_rows, num_columns);

    for (Eigen::Index column_index = 0; column_index < num_columns; ++column_index) {
        const auto& source_column = get_checked_column_or_exit(
            m_data_vec, column_indices[static_cast<size_t>(column_index)]);
        for (Eigen::Index row_index = 0; row_index < num_rows; ++row_index) {
            matrix(row_index, column_index) =
                string_to_double(source_column[static_cast<size_t>(row_index)]);
        }
    }
}

/**
 * @brief Build the fixed-effect design matrix and aligned response matrix.
 *
 * @param covariate_columns Phenotype columns treated as numeric covariates.
 * @param class_columns Phenotype columns treated as categorical effects.
 * @param fixed_effect_column_offsets Cumulative column offsets for each effect group.
 * @param fixed_effect_names Descriptive labels for the resulting fixed-effect columns.
 * @param trait_columns Phenotype columns to extract into the response matrix.
 * @param response_matrix Output matrix for the trait columns.
 * @return MatrixXd Fixed-effect design matrix.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
MatrixXd PHEN::build_fixed_effect_design_matrix(const vector<std::int64_t>& covariate_columns,
                                                const vector<std::int64_t>& class_columns,
                                                vector<std::int64_t>& fixed_effect_column_offsets,
                                                vector<string>& fixed_effect_names,
                                                const vector<std::int64_t>& trait_columns,
                                                MatrixXd& response_matrix) const{
    if (m_data_vec.empty()) {
        fatal_error("Cannot build a fixed-effect design matrix: no phenotype data are loaded.");
    }

    const Eigen::Index num_rows = static_cast<Eigen::Index>(m_data_vec[0].size());
    if (!trait_columns.empty()) {
        PHEN::get_columns_as_matrix(trait_columns, response_matrix);
    } else {
        response_matrix.resize(num_rows, 0);
    }

    MatrixXd covariate_matrix;
    if (!covariate_columns.empty()) {
        PHEN::get_columns_as_matrix(covariate_columns, covariate_matrix);
    } else {
        covariate_matrix.resize(num_rows, 0);
    }

    vector<MatrixXd> class_design_blocks;
    vector<vector<string>> class_level_names;
    class_design_blocks.reserve(class_columns.size());
    class_level_names.reserve(class_columns.size());

    Eigen::Index total_num_design_columns = 1 + covariate_matrix.cols();
    size_t total_num_fixed_effect_names = 1 + static_cast<size_t>(covariate_matrix.cols());
    for (std::int64_t class_column : class_columns) {
        vector<string> class_levels;
        MatrixXd class_block = PHEN::build_design_block_for_class(class_column, class_levels);
        total_num_design_columns += class_block.cols();
        total_num_fixed_effect_names += class_levels.size();
        class_design_blocks.push_back(std::move(class_block));
        class_level_names.push_back(std::move(class_levels));
    }

    fixed_effect_column_offsets.clear();
    fixed_effect_column_offsets.reserve(2 + covariate_columns.size() + class_columns.size());
    fixed_effect_names.clear();
    fixed_effect_names.reserve(total_num_fixed_effect_names);

    MatrixXd design_matrix(num_rows, total_num_design_columns);
    Eigen::Index current_column = 0;

    design_matrix.col(current_column).setOnes();
    fixed_effect_column_offsets.push_back(0);
    ++current_column;
    fixed_effect_column_offsets.push_back(static_cast<std::int64_t>(current_column));
    fixed_effect_names.push_back("0 mu 1");

    if (covariate_matrix.cols() > 0) {
        design_matrix.middleCols(current_column, covariate_matrix.cols()) = covariate_matrix;
        for (std::int64_t covariate_column : covariate_columns) {
            fixed_effect_column_offsets.push_back(static_cast<std::int64_t>(current_column + 1));
            fixed_effect_names.push_back("1 " + get_checked_header_name_or_exit(m_head_column_vec, covariate_column) + " 1");
            ++current_column;
        }
    }

    for (size_t class_idx = 0; class_idx < class_columns.size(); ++class_idx) {
        const auto& class_block = class_design_blocks[class_idx];
        design_matrix.middleCols(current_column, class_block.cols()) = class_block;
        current_column += class_block.cols();
        fixed_effect_column_offsets.push_back(static_cast<std::int64_t>(current_column));

        const auto& class_levels = class_level_names[class_idx];
        const string& class_name =
            get_checked_header_name_or_exit(m_head_column_vec, class_columns[class_idx]);
        for (const auto& class_level : class_levels) {
            fixed_effect_names.push_back("2 " + class_name + " " + class_level);
        }
    }

    return design_matrix;
}


/**
 * @brief Build the dummy-coded design block for one categorical effect.
 *
 * If `class_level_vec` is empty, levels are inferred from the current data in
 * first-seen order.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
MatrixXd PHEN::build_design_block_for_class(std::int64_t effect_column,
                                            vector<string>& class_levels) const{
    const auto& class_values = get_checked_column_or_exit(m_data_vec, effect_column);
    if (class_levels.empty()) {
        class_levels = vector_unique(class_values);
    }

    std::unordered_map<string, Eigen::Index> class_level_to_column;
    class_level_to_column.reserve(class_levels.size());
    for (Eigen::Index level_index = 0; level_index < static_cast<Eigen::Index>(class_levels.size()); ++level_index) {
        const auto [it, inserted] =
            class_level_to_column.emplace(class_levels[static_cast<size_t>(level_index)], level_index);
        if (!inserted) {
            fatal_error("Duplicated class level '{}' is not allowed for fixed-effect design construction.",
                          class_levels[static_cast<size_t>(level_index)]);
        }
    }

    const Eigen::Index num_rows = static_cast<Eigen::Index>(class_values.size());
    MatrixXd design_block = MatrixXd::Zero(num_rows, static_cast<Eigen::Index>(class_levels.size()));
    for (Eigen::Index row_index = 0; row_index < num_rows; ++row_index) {
        const auto it = class_level_to_column.find(class_values[static_cast<size_t>(row_index)]);
        if (it == class_level_to_column.end()) {
            fatal_error("Categorical level '{}' in column {} is missing from the provided class level set.",
                          class_values[static_cast<size_t>(row_index)], effect_column);
        }

        design_block(row_index, it->second) = 1.0;
    }

    return design_block;
}



/**
 * @brief Append sparse design-matrix triplets for one categorical random effect.
 *
 * The phenotype column `effect_column` provides the row-side class labels. The
 * column-side level set comes from `class_levels`; if it is empty, levels are
 * inferred from the data in first-seen order.
 *
 * Used in:
 * - currently no in-tree call sites were found
 */
void PHEN::append_random_effect_design_triplets(std::int64_t effect_column,
                                                vector<string>& class_levels,
                                                std::int64_t row_offset,
                                                std::int64_t column_offset,
                                                vector<Eigen::Triplet<double>>& triplets) const{
    const auto& class_values = get_checked_column_or_exit(m_data_vec, effect_column);
    if (class_levels.empty()) {
        class_levels = vector_unique(class_values);
    }

    std::unordered_map<string, std::int64_t> class_level_to_index;
    class_level_to_index.reserve(class_levels.size());
    for (std::int64_t level_index = 0; level_index < static_cast<std::int64_t>(class_levels.size()); ++level_index) {
        const auto [it, inserted] = class_level_to_index.emplace(class_levels[static_cast<size_t>(level_index)],
                                                                 level_index);
        if (!inserted) {
            fatal_error("Duplicated class level '{}' is not allowed for random-effect design construction.",
                          class_levels[static_cast<size_t>(level_index)]);
        }
    }

    triplets.reserve(triplets.size() + class_values.size());
    for (size_t row_index = 0; row_index < class_values.size(); ++row_index) {
        const auto it = class_level_to_index.find(class_values[row_index]);
        if (it != class_level_to_index.end()) {
            triplets.emplace_back(row_offset + static_cast<std::int64_t>(row_index),
                                  column_offset + it->second,
                                  1.0);
        }
    }
}



/**
 * @brief Remove dependent columns in the xmat
 * 
 * @param design_matrix Design matrix for fixed effects.
 * @param fixed_effect_column_offsets Updated fixed-effect column offsets.
 * @param fixed_effect_names Fixed-effect labels; removed columns are marked.
 * @param removed_column_indices Zero-based design-matrix columns removed to make the matrix full rank.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
void PHEN::remove_dependent_design_columns(MatrixXd& design_matrix,
                                           vector<std::int64_t>& fixed_effect_column_offsets,
                                           vector<string>& fixed_effect_names,
                                           vector<std::int64_t>& removed_column_indices) const{
    removed_column_indices.clear();
    if (design_matrix.cols() == 0) {
        return;
    }

    HouseholderQR<MatrixXd> qr(design_matrix);
    const Eigen::Index num_columns = design_matrix.cols();
    const Eigen::Index num_rows = design_matrix.rows();
    const Eigen::Index diagonal_count = std::min(num_rows, num_columns);
    constexpr double kTolerance = 1.0e-6;

    removed_column_indices.reserve(static_cast<size_t>(num_columns));
    const auto& qr_storage = qr.matrixQR();
    for (Eigen::Index diag_index = 0; diag_index < diagonal_count; ++diag_index) {
        if (std::fabs(qr_storage(diag_index, diag_index)) < kTolerance) {
            removed_column_indices.push_back(static_cast<std::int64_t>(diag_index));
        }
    }
    for (Eigen::Index column_index = diagonal_count; column_index < num_columns; ++column_index) {
        removed_column_indices.push_back(static_cast<std::int64_t>(column_index));
    }

    if (removed_column_indices.empty()) {
        return;
    }

    remove_col(design_matrix, removed_column_indices);

    for (auto& column_offset : fixed_effect_column_offsets) {
        const auto removed_before_offset =
            std::lower_bound(removed_column_indices.begin(), removed_column_indices.end(), column_offset) -
            removed_column_indices.begin();
        column_offset -= static_cast<std::int64_t>(removed_before_offset);
    }

    for (std::int64_t removed_column : removed_column_indices) {
        if (removed_column < 0 || static_cast<size_t>(removed_column) >= fixed_effect_names.size()) {
            fatal_error("Removed design-matrix column {} is out of range for {} fixed-effect labels.",
                          removed_column, fixed_effect_names.size());
        }
        fixed_effect_names[static_cast<size_t>(removed_column)] += " 0 0";
    }
}
