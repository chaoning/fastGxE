#pragma once


#include <cstdint>

#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigen>

using std::string;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class PHEN {
public:
    PHEN() = default;

    void read_all_columns_raw(const string& pheno_file, vector<vector<string>>& data_vec);
    void read_selected_columns_raw(const string& pheno_file,
                                   const vector<std::int64_t>& target_cols,
                                   vector<vector<string>>& data_vec);
    void update_keep_indices_by_grm_ids(const vector<std::int64_t>& random_grm_vec,
                                        const vector<vector<string>>& grm_id_used_vec_vec);
    void update_keep_indices_by_non_missing_values(const vector<std::int64_t>& col_used_vec,
                                                   const vector<string>& missing_in_data_vec);
    void apply_keep_indices_to_data();
    void reorder_records_by_sample_id_order(const vector<string>& sample_id_order_vec);
    void read_and_filter(const string& pheno_file,
                         const vector<std::int64_t>& random_grm_vec,
                         const vector<vector<string>>& grm_id_used_vec_vec,
                         const vector<std::int64_t>& col_used_vec,
                         const vector<string>& missing_in_data_vec);

    MatrixXd build_fixed_effect_design_matrix(const vector<std::int64_t>& covariate_columns,
                                              const vector<std::int64_t>& class_columns,
                                              vector<std::int64_t>& fixed_effect_column_offsets,
                                              vector<string>& fixed_effect_names,
                                              const vector<std::int64_t>& trait_columns,
                                              MatrixXd& response_matrix) const;
    void remove_dependent_design_columns(MatrixXd& design_matrix,
                                         vector<std::int64_t>& fixed_effect_column_offsets,
                                         vector<string>& fixed_effect_names,
                                         vector<std::int64_t>& removed_column_indices) const;

    MatrixXd build_design_block_for_class(std::int64_t effect_column,
                                          vector<string>& class_levels) const;

    void append_random_effect_design_triplets(std::int64_t effect_column,
                                              vector<string>& class_levels,
                                              std::int64_t row_offset,
                                              std::int64_t column_offset,
                                              vector<Eigen::Triplet<double>>& triplets) const;

    void get_stored_column_names(vector<string>& column_names) const;
    void get_given_column(std::int64_t column_index, vector<string>& values) const;
    void get_given_column(std::int64_t column_index, VectorXd& values) const;
    void get_columns_as_matrix(const vector<std::int64_t>& column_indices, MatrixXd& matrix) const;

private:
    vector<string> m_head_column_vec;
    vector<vector<string>> m_data_vec;
    vector<std::int64_t> m_index_keep_vec;
};
