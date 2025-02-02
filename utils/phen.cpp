/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-07-06 11:03:03
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-31 19:45:56
 */


#include<iostream>
#include<fstream>
#include<set>
#include<vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <spdlog/spdlog.h>

#include <boost/algorithm/string.hpp>

#include "phen.hpp"
#include "iterator_utils.hpp"
#include "string_utils.hpp"
#include "EigenMatrix_utils.hpp"


using std::cin;
using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;
using std::ifstream;
using std::ofstream;
using Eigen::HouseholderQR;


PHEN::PHEN(){
}


void PHEN::read_full(string pheno_file, vector<vector<string>>& data_vec){
    ifstream fin(pheno_file);
    if(!fin.is_open()) {
        spdlog::error("Fail to open the phenotypic file");
        exit(1);
    }
    
    // head line
    string one_line;
    getline(fin, one_line); // skip header line
    m_head_column_vec.clear(); // ! class member
    boost::split(m_head_column_vec, one_line, isspace, boost::token_compress_on);
    long long num_col_head = m_head_column_vec.size();

    // read
    data_vec.resize(num_col_head + 1);
    long long ith_line = 0;
    while (getline(fin, one_line)) {
        vector<string> tmp_vec;
        boost::split(tmp_vec, one_line, isspace, boost::token_compress_on);
        if(tmp_vec.size() != num_col_head){ // locate line with different columns compared with head lines.
            cout << "Check this row: " << one_line <<", this row has different number of columns with head lines" << endl;
            exit(1);
        }
        ith_line++;
        tmp_vec.push_back(std::to_string(ith_line)); // add line NO.
        for(int i = 0; i < tmp_vec.size(); i++){
            data_vec[i].push_back(tmp_vec[i]);
        }
    }
    fin.close();
}


void PHEN::index_keep_byGRM(vector<long long> random_grm_vec, vector<vector<string>> grm_id_used_vec_vec){
    
    for(int i = 0; i < random_grm_vec.size(); i++){
        set<string> grm_id_set(grm_id_used_vec_vec[i].begin(), grm_id_used_vec_vec[i].end());
        vector<string> data_id_vec = m_data_vec[random_grm_vec[i]];
        set<string> data_id_set(data_id_vec.begin(), data_id_vec.end());
        set<string> set_inter = set_intersection_(grm_id_set, data_id_set);
        set<string> set_diff = set_difference_(data_id_set, grm_id_set);
        
        map<string, long long> data_id_map;
        for(auto val:set_inter){
            data_id_map[val] = 1;
        }

        for(auto val:set_diff){
            data_id_map[val] = -1;
        }
        
        
        long long k = 0;
        vector<long long> tmp_index_keep_vec;
        for(auto val:data_id_vec){
            k++;
            if(data_id_map[val] == 1){
                tmp_index_keep_vec.push_back(k);
            }
        }

        m_index_keep_vec = set_intersection_(m_index_keep_vec, tmp_index_keep_vec);

    }
}

void PHEN::index_keep_deleteNA(vector<long long> col_used_vec, vector<string> missing_in_data_vec){
    for(auto val:col_used_vec){
        auto tmp_vec = m_data_vec[val];
        long long k = 0;
        vector<long long> tmp_index_keep_vec;
        for(auto val:tmp_vec){
            k++;
            if(!is_nan(val, missing_in_data_vec)){
                tmp_index_keep_vec.push_back(k);
            }
        }
        m_index_keep_vec = set_intersection_(m_index_keep_vec, tmp_index_keep_vec);
    }
}

void PHEN::update_data_vec(){
    for(auto& tmp_vec:m_data_vec){
        vector<string> tmp_vec2(m_index_keep_vec.size());
        long long ith = 0;
        for(auto val:m_index_keep_vec){
            tmp_vec2[ith] = tmp_vec[val-1];
            ith++;
        }
        tmp_vec = tmp_vec2;
    }
}

/**
 * @brief read the phenotypic file
 * 
 * @param pheno_file phenotypic file
 * @param random_grm_vec random effects linked with grm files
 * @param grm_id_used_vec_vec iids vector for grm files
 * @param col_used_vec used columns
 * @param missing_in_data_vec missing data
 */
void PHEN::read(string pheno_file, vector<long long> random_grm_vec, vector<vector<string>> grm_id_used_vec_vec, vector<long long> col_used_vec, vector<string> missing_in_data_vec){
    
    // read all data
    m_data_vec.clear(); // ! class member
    PHEN::read_full(pheno_file, m_data_vec);
    long long num_raw_records = m_data_vec[0].size();
    m_index_keep_vec.resize(num_raw_records); // ! class member
    long long ith = 0;
    for(auto& val:m_index_keep_vec){
        ith++;
        val = ith;
    }

    spdlog::info("Map the data with the GRM", 0);
    PHEN::index_keep_byGRM(random_grm_vec, grm_id_used_vec_vec);
    spdlog::info("Delete lines with NA values", 0);
    PHEN::index_keep_deleteNA(col_used_vec, missing_in_data_vec);
    
    double ratio = m_index_keep_vec.size() * 1.0 / num_raw_records;
    if(ratio < 0.1){
        spdlog::error("Only " + std::to_string(ratio * 100) + "% records are remained for analysis, which is less than 10%, Please check the data.");
        exit(1);
    }

    spdlog::info("Update the dataframe", 0);
    PHEN::update_data_vec();
}

void PHEN::update_use_given_idOrder(vector<string> given_id_order_vec){
    
    vector<string> id_in_data_keep_vec;
    this->get_given_column(0, id_in_data_keep_vec);
    map<string, long long> id_map;
    long long k = 0;
    for(auto tmp:id_in_data_keep_vec){
        id_map[tmp] = k++;
    }

    
    vector<long long> id_index;
    for(auto tmp:given_id_order_vec){
        id_index.push_back(id_map[tmp]);
    }

    for(auto& tmp_vec:m_data_vec){
        vector<string> tmp_vec2(id_index.size());
        long long ith = 0;
        for(auto val:id_index){
            tmp_vec2[ith] = tmp_vec[val];
            ith++;
        }
        tmp_vec = tmp_vec2;
    }
    
}

/**
 * @brief get the given column
 * 
 * @param columnNo 
 * @param vec 
 */
void PHEN::get_given_column(long long columnNo, vector<string>& vec){
    vec = m_data_vec[columnNo];
}


vector<string> PHEN::get_used_record_number(){
    vector<string> vec = m_data_vec.back();
    return vec;
}


void PHEN::get_given_column(long long columnNo, vector<double>& vec){
    vector<string> tmp_vec = m_data_vec[columnNo];
    char *endptr = new char[1000];
    for(auto tmp_val:tmp_vec){
        double val = string_to_double(tmp_val, endptr);
        vec.push_back(val);
    }
    delete[] endptr;
}

void PHEN::get_given_column(long long columnNo, VectorXd& Vec){
    vector<string> tmp_vec = m_data_vec[columnNo];
    Vec.resize(tmp_vec.size());

    char *endptr = new char[1000];
    long long i = 0;
    for(auto tmp_val:tmp_vec){
        double val = string_to_double(tmp_val, endptr);
        Vec(i++) = val;
    }
    delete[] endptr;
}

void PHEN::get_given_columns_double(vector<long long> column_vec, MatrixXd& Mat){
    long long nrows = m_data_vec[0].size();
    long long ncols = column_vec.size();
    Mat = MatrixXd::Zero(nrows, ncols);
    char *endptr = new char[1000];
    for(long long i = 0; i < ncols; i++){
        vector<string> tmp_vec = m_data_vec[column_vec[i]];
        for(long long j = 0; j < nrows; j++){
            Mat(j, i) = string_to_double(tmp_vec[j], endptr);
        }
    }
    delete[] endptr;
}

/**
 * @brief get the design matrix and index for fixed effects
 * 
 * @param covariate_vec covariate in the model
 * @param class_vec class variate in the model
 * @param index_fixed_effect_vec index vector for fixed effect
 * @param fixed_effect_name_vec a string contains three columns splited with space: mean(0), covariate(1) or class(2); term name; levels
 * @param trait_vec vector of column numbers for traits
 * @param y matrix to store the phenotypic values
 * @return MatrixXd design matrix
 */
MatrixXd PHEN::dmat(vector<long long> covariate_vec, vector<long long> class_vec, vector <long long>& index_fixed_effect_vec, vector <string>& fixed_effect_name_vec,
            vector<long long> trait_vec, MatrixXd& y){
    // y
    long long num_row_in_data = m_data_vec[0].size();
    long long num_trait = trait_vec.size();
    // cout << "The number of traits is: " << num_trait << endl;
    if(num_trait > 0){
        y.resize(num_row_in_data, num_trait);
        for(int i = 0; i < num_trait; i++){
            VectorXd tmp_Vec;
            PHEN::get_given_column(trait_vec[i], tmp_Vec);
            y.col(i) = tmp_Vec;
        }
    }

    
    // xmat
    index_fixed_effect_vec.clear();
    fixed_effect_name_vec.clear();
    MatrixXd xmat;
    index_fixed_effect_vec.push_back(0);
    // 1s
    xmat.setOnes(num_row_in_data, 1);
    index_fixed_effect_vec.push_back(xmat.cols());
    fixed_effect_name_vec.push_back("0 mu 1"); // 0: mean; 1, term; 1, levels


    // covariate
    long long num_covariate = covariate_vec.size();
    if(num_covariate > 0){
        int xmat_df = xmat.cols();
        xmat.conservativeResize(num_row_in_data, xmat_df + num_covariate);
        vector <string> element_uni_vec;
        for(int i = 0; i < num_covariate; i++){
            MatrixXd xmat_tmp = PHEN::dmat_for_one_effect(1, covariate_vec[i], element_uni_vec);
            xmat.col(xmat_df + i) = xmat_tmp.col(0);
            index_fixed_effect_vec.push_back(index_fixed_effect_vec.back() + 1);
            fixed_effect_name_vec.push_back("1 " + element_uni_vec[i] + " 1"); // 1 for covariate; term name; levels
        }
    }

    // class
    long long num_class = class_vec.size();
    if(num_class > 0){
        for(long long i = 0; i < num_class; i++){
            vector <string> element_uni_vec;
            MatrixXd xmat_tmp = PHEN::dmat_for_one_effect(2, class_vec[i], element_uni_vec);
            int xmat_df = xmat.cols();
            xmat.conservativeResize(num_row_in_data, xmat_df + xmat_tmp.cols());
            xmat.rightCols(xmat_tmp.cols()) = xmat_tmp;

            index_fixed_effect_vec.push_back(xmat.cols());
            for(vector<string>::iterator it = element_uni_vec.begin(); it != element_uni_vec.end(); it++){
                fixed_effect_name_vec.push_back("2 " + m_head_column_vec[class_vec[i]] + " " + *it); // 2 for class; i+1, order; levels
            }
        }
    }
    return xmat;
}


/**
 * @brief get the design matrix and index for one fixed effect
 * 
 * @param type covariate (1) or class(2)
 * @param effect_col column number 
 * @param element_uni_vec vector of unique elements, given or not
 * @return MatrixXd 
 */
MatrixXd PHEN::dmat_for_one_effect(int type, long long effect_col, vector<string>& element_uni_vec){
    MatrixXd mat;
    if(type == 1){
        mat = PHEN::dmat_for_one_covariate(effect_col);
        element_uni_vec.push_back(m_head_column_vec[effect_col]);
    }else if(type == 2){
        mat = PHEN::dmat_for_one_class(effect_col, element_uni_vec);
    }else{
        cout << "The parameter type must be 1 (covariate) or 2 (class)" << endl;
        exit(1);
    }
    return mat;
}



MatrixXd PHEN::dmat_for_one_class(long long effect_col, vector<string>& class_level_vec){
    vector<string> tmp_vec = m_data_vec[effect_col];
    if(class_level_vec.empty()){
        class_level_vec = vector_unique(tmp_vec);
    }

    map<string, long long> class_map;
    long long val = 0;
    for(auto tmp_val:class_level_vec){
        class_map[tmp_val] = val++;
    }

    long long nrow = m_data_vec[0].size();
    MatrixXd mat = MatrixXd::Zero(nrow, class_level_vec.size());
    for(long long j = 0; j < nrow; j++){
        mat(j, class_map[tmp_vec[j]]) = 1.0;
    }
    return mat;
}



MatrixXd PHEN::dmat_for_one_covariate(long long effect_col){
    VectorXd tmp_Vec;
    PHEN::get_given_column(effect_col, tmp_Vec);
    MatrixXd mat(tmp_Vec.size(), 1);
    mat.col(0) = tmp_Vec;
    return mat;
}



/**
 * @brief Obtain the design matrix for one random effect (class)
 * 
 * @param effect_col 
 * @param element_uni_vec
 * @return SparseMatrix < double > 
 */
void PHEN::dmat_for_one_random(long long effect_col, vector<string>& element_uni_vec,
            long long start_row, long long start_column, vector < Eigen::Triplet < double > >& triplet_vec){
    long long nrow = m_data_vec.size();
    vector<string> element_vec = m_data_vec[effect_col];

    if(element_uni_vec.empty()){
        element_uni_vec = vector_unique(element_vec);
    }
    
    PHEN::dmat_rowid_colid(element_vec, element_uni_vec, start_row, start_column, triplet_vec);
}


/**
 * @brief design matrix with row id and column id
 * 
 * @param row_id_vec row id, id in data
 * @param col_id_vec column id, effect id, may no in the data
 * @param row_start 
 * @param col_start 
 * @param triplet_vec 
 */
void PHEN::dmat_rowid_colid(vector<string>& row_id_vec, vector<string>& col_id_vec, long long row_start, long long col_start, 
        vector < Eigen::Triplet < double > >& triplet_vec){
    map<string, long long> id_map;
    long long val = 0;
    for(vector<string>::iterator it = col_id_vec.begin(); it != col_id_vec.end(); it++){
        id_map[*it] = val++;
    }
    eigen_assert(id_map.size() == col_id_vec.size() && "Duplicated iids are not allowed in the effect vector");


    for(long long i = 0; i < row_id_vec.size(); i++){
        auto it = id_map.find(row_id_vec[i]);
        if(it != id_map.end()){
            long long j = it->second;
            triplet_vec.emplace_back(row_start + i, col_start + j, 1.0);
        }
    }
}



/**
 * @brief Remove dependent columns in the xmat
 * 
 * @param xmat design matrix for fixed effects
 * @param index_fixed_effect_vec modified index vector
 * @param fixed_effect_name_vec for a removed effect, add "0 0"
 * @param column_to_remove_vec columns are removed to guarantee xmat is full rank
 */
void PHEN::dmat_column_full_rank(MatrixXd& xmat, vector <long long>& index_fixed_effect_vec, vector <string>& fixed_effect_name_vec, vector<long long>& column_to_remove_vec){
    // delete dependent columns
    HouseholderQR<MatrixXd> qr;
    qr.compute(xmat);
    column_to_remove_vec.clear();
    MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    for(long long i = 0; i < R.cols(); i++){ // locate dependent columns
        if(std::fabs(R(i, i)) < 1.0e-6){
            column_to_remove_vec.push_back(i);
        }
    }
    if(column_to_remove_vec.size() != 0){
        remove_col(xmat, column_to_remove_vec);
    }
    // Modify the the index for fixed effect
    for(long long i = column_to_remove_vec.size()-1; i >= 0; i--){
        for(long long j = 0; j < index_fixed_effect_vec.size(); j++){
            if(index_fixed_effect_vec[j] > column_to_remove_vec[i]){
                index_fixed_effect_vec[j] -= 1; 
            }
        }
    }
    // Dependent effect set effect and se to zero
    for(long long i = 0; i < column_to_remove_vec.size(); i++){
        fixed_effect_name_vec[column_to_remove_vec[i]] += " 0 0";
    }
}


/**
 * @brief get the number of used records
 * 
 * @return long long 
 */
long long PHEN::get_num_used_records(){
    return m_data_vec[0].size();
}
