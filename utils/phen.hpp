/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-31 15:05:20
 * @LastEditTime: 2025-01-31 19:45:03
 * @LastEditors: Chao Ning
 */

#pragma once

#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigen>

using std::string;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

class PHEN{
    public:
        PHEN();
        void read_full(string pheno_file, vector<vector<string>>& data_vec);
        void index_keep_byGRM(vector<long long> random_grm_vec, vector<vector<string>> grm_id_used_vec_vec);
        void index_keep_deleteNA(vector<long long> col_used_vec, vector<string> missing_in_data_vec);
        void update_data_vec();
        void update_use_given_idOrder(vector<string> given_id_order_vec);
        void read(string pheno_file, vector<long long> random_grm_vec, vector<vector<string>> grm_id_used_vec_vec, vector<long long> col_used_vec, vector<string> missing_in_data_vec);
        
        
        MatrixXd dmat(vector<long long> covariate_vec, vector<long long> class_vec, vector <long long>& index_fixed_effect_vec, vector <string>& fixed_effect_name_vec,
            vector<long long> trait_vec, MatrixXd &y);
        void dmat_column_full_rank(MatrixXd& xmat, vector <long long>& index_fixed_effect_vec, vector <string>& fixed_effect_name_vec, vector<long long>& column_to_remove_vec);

        MatrixXd dmat_for_one_effect(int type, long long effect_col, vector<string>& class_vec);
        MatrixXd dmat_for_one_class(long long effect_col, vector<string>& class_level_vec);
        MatrixXd dmat_for_one_covariate(long long effect_col);

        void dmat_for_one_random(long long effect_col, vector<string>& element_uni_vec,
            long long start_row, long long start_column, vector < Eigen::Triplet < double > >& triplet_vec);
        
        

        void dmat_rowid_colid(vector<string>& row_id_vec, vector<string>& col_id_vec, long long row_start, long long col_start, 
            vector < Eigen::Triplet < double > >& triplet_vec);


        void get_given_column(long long columnNo, vector<string>& vec);
        void get_given_column(long long columnNo, vector<double>& vec);
        void get_given_column(long long columnNo, VectorXd& Vec);
        void get_given_columns_double(vector<long long> column_vec, MatrixXd& Mat);
        vector<string> get_used_record_number();

        long long get_num_used_records();

    private:
        vector<string> m_head_column_vec;
        vector<vector<string> > m_data_vec;
        vector<long long> m_index_keep_vec;
};
