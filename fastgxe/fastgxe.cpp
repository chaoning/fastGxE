/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-18 21:22:40
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-17 16:33:55
 */

#define EIGEN_USE_MKL_ALL  // must be before include Eigen

#include <getopt.h>
#include "mkl.h"
#include "omp.h"

#include <set>
#include <map>
#include <cmath>
#include <string>
#include <random>
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
#include "../utils/random_num.hpp"
#include "fastgxe.hpp"


using std::set;
using std::map;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::endl;
using Eigen::HouseholderQR;
using Eigen::ColPivHouseholderQR;


fastGxE::fastGxE(){

}

fastGxE::~fastGxE() {
}

/**
 * @Description: Read the data file and GRM file
 */
int fastGxE::pre_data(string out_file, string data_file, string agrm_file, vector<long long> covariate_arr, 
        vector<long long> class_arr, vector<long long> bye_arr,  vector<long long> trait, 
        vector<string> missing_in_data_vec){
    spdlog::info("Start analysis...");
    m_out_file = out_file;

    m_num_trait = trait.size();
    long long num_covariate = covariate_arr.size();
    long long num_class = class_arr.size();
    this->m_num_bye = bye_arr.size();
    spdlog::info("The number of analyzed traits: {}", m_num_trait);
    spdlog::info("The number of covariates: {}", num_covariate);
    spdlog::info("The number of classes: {}", num_class);
    spdlog::info("The number of interaction environments: {}", m_num_bye);

    // id in grm
    spdlog::info("Read iids in GRM");
    ProcessGRM ProcessGRMA;
    vector<string> id_in_gmat_vec = ProcessGRMA.read_grm_id(agrm_file + ".agrm");
    long long num_id_in_grm = id_in_gmat_vec.size();
    spdlog::info("The sample size in the GRM is: {}" + num_id_in_grm);

    // read data file
    spdlog::info("Read data file");
    PHEN phenoA;
    // vector<long long> tmp_bye_arr;
    // if(bye_arr.size() > 0) tmp_bye_arr = {bye_arr[0]};
    vector<long long> col_used_vec = vector_merged({covariate_arr, class_arr, bye_arr, trait});
    phenoA.read(data_file, {0}, {id_in_gmat_vec}, col_used_vec, missing_in_data_vec);
    
    // id in data
    this->m_id_in_data_vec.clear();
    phenoA.get_given_column(0, m_id_in_data_vec);
    set<string> id_in_data_st(m_id_in_data_vec.begin(), m_id_in_data_vec.end());
    if(m_id_in_data_vec.size() > id_in_data_st.size()){ // check whether duplicated ids exist
        spdlog::error("Duplicated sample IDs exist in the data file!");
        exit(1);
    }
    m_num_id = m_id_in_data_vec.size();
    spdlog::info("The number of used iids in data file: {}", m_num_id);


    // group
    spdlog::info("Read the group files: {}.agrm.group", agrm_file);
    map<string, long long> id_map;
    for(auto tmp:m_id_in_data_vec){ // ! Only keep the iids both in data and GRM
        id_map[tmp] = 1;
    }

    ifstream fin(agrm_file + ".agrm.group");
    if(!fin.is_open()){
        spdlog::error("Fail to open the {}.agrm.group", agrm_file);
        exit(1);
    }
    string line;
    vector<vector<long long>> group_vec;
    while(std::getline(fin, line)){
        process_line(line);
        vector <string> vec = split_string(line);
        if(id_map.find(vec[1]) != id_map.end())
            group_vec.push_back({std::atoll(vec[0].c_str()), std::atoll(vec[2].c_str()), std::atoll(vec[3].c_str())});
    }
    fin.close();


    // Re-calculate the group size. Some iids in the GRM are unused.
    spdlog::info("Re-calculate the group size");
    map<long long, long long> group_size_map;
    for(auto tmp_vec:group_vec){
        auto it = group_size_map.find(tmp_vec[1]);
        if(it == group_size_map.end()){
            group_size_map[tmp_vec[1]] = 1;
        }else{
            group_size_map[tmp_vec[1]] += 1;
        }
    }
    for(auto &tmp_vec:group_vec){
        tmp_vec[2] = group_size_map[tmp_vec[1]];
    }
    

    // Sort
    spdlog::info("Sort by group size and group index");
    sort(group_vec.begin(), group_vec.end(), [](vector<long long>a, vector<long long>b){
        if(a[2] != b[2])
            return a[2] < b[2];
        else{
            if(a[1] != b[1])
                return a[1] < b[1];
            else
                return a[0] < b[0];
        }
    });

    // update dataframe
    spdlog::info("Update the dataframe using the current sample id orders");
    vector<string> id_in_gmat_used_reorder_vec;
    for(auto tmp_vec:group_vec){
        id_in_gmat_used_reorder_vec.push_back(id_in_gmat_vec[tmp_vec[0] - 1]);
    }
    phenoA.update_use_given_idOrder(id_in_gmat_used_reorder_vec);
    this->m_id_in_data_vec.clear();
    phenoA.get_given_column(0, m_id_in_data_vec); // ! update id order

    spdlog::info("Design matrix for fixed effects");
    vector <long long> index_fixed_effect_vec;
    vector <string> fixed_effect_name_vec; // name for fixed effect
    m_xmat = phenoA.dmat(covariate_arr, class_arr, index_fixed_effect_vec, fixed_effect_name_vec, trait, m_y);
    vector<long long> column_to_remove_vec;
    phenoA.dmat_column_full_rank(m_xmat, index_fixed_effect_vec, fixed_effect_name_vec, column_to_remove_vec); // full rank of _xmat
    if (!column_to_remove_vec.empty()) {
        spdlog::info("Dependent columns of xmat are removed: ");
        std::string removed_columns;
        for (const auto& col : column_to_remove_vec) {
            removed_columns += std::to_string(col) + " ";
        }
        spdlog::info("{}", removed_columns);
    }

    spdlog::info("Get interaction environment covariates");
    if(!bye_arr.empty()){
        phenoA.get_given_columns_double(bye_arr, this->m_bye_mat);
    }
    
    spdlog::info("Store sample index within subgroup");
    vector<long long> group_one_index_vec; // Store all groups with only one individuals
    vector<vector<long long>> group_index_vec;
    vector<long long> tmp_group_vec = {0};
    long long igroup = 0;
    for(auto& tmp:group_vec){
        tmp_group_vec.push_back(tmp[1]);
        long long n0 = tmp_group_vec[tmp_group_vec.size() - 2];
        long long n1 = tmp_group_vec[tmp_group_vec.size() - 1];

        if(tmp[2] == 1){
            group_one_index_vec.push_back(tmp[0]);
        }else{
            if(n0 != n1){
                group_index_vec.push_back({tmp[0]});
            }else{
                group_index_vec[group_index_vec.size() - 1].push_back(tmp[0]);
            }
        }
    }


    map<string, vector<long long>> grm_index_group_map; // key: Pair-wise sample index within subgroup; val: {group order, MatrixXd index0, MatrixXd index1}
    m_grm_mat_group_vec = {MatrixXd::Zero(group_one_index_vec.size(), 1)};
    long long val = 0;
    for(auto tmp:group_one_index_vec){
        string key_str = std::to_string(tmp) + "_" + std::to_string(tmp);
        grm_index_group_map[key_str] = {0, val, 0};
        val++;
    }

    for(long long i = 0; i < group_index_vec.size(); i++){
        m_grm_mat_group_vec.push_back(MatrixXd::Zero(group_index_vec[i].size(), group_index_vec[i].size()));
        for(long long m = 0; m < group_index_vec[i].size(); m++){
            for(long long n = 0; n < group_index_vec[i].size(); n++){
                string key_str = std::to_string(group_index_vec[i][m]) + "_" + std::to_string(group_index_vec[i][n]);
                grm_index_group_map[key_str] = {i+1, m, n};
            }
        }
    }

    spdlog::info("Read the GRMs");
    
    std::ifstream fin2(agrm_file + ".agrm.sp.bin", std::ios::binary);
    if(!fin2.is_open()){
        spdlog::error("Fail to open: {}.agrm.sp.bin", agrm_file);
        exit(1);
    }

    double tmp_val;
    long long index0, index1;

    while(fin2.read((char*)&index0, sizeof(long long)) && 
            fin2.read((char*)&index1, sizeof(long long)) && fin2.read((char*)&tmp_val, sizeof(double))){
        string key_str = std::to_string(index0 + 1) + "_" + std::to_string(index1 + 1);
        if(grm_index_group_map.find(key_str) != grm_index_group_map.end()){
            vector<long long> tmp_vec = grm_index_group_map[key_str];
            if(tmp_vec[0] != 0){
                m_grm_mat_group_vec[tmp_vec[0]](tmp_vec[1], tmp_vec[2]) = m_grm_mat_group_vec[tmp_vec[0]](tmp_vec[2], tmp_vec[1]) = tmp_val;
            }else{
                m_grm_mat_group_vec[tmp_vec[0]](tmp_vec[1], tmp_vec[2]) = tmp_val; // tmp_vec[2] = 0 now
            }
        }
    }
    fin2.close();

    
    return 0;
}








/**
 * @Description: Sparse GRM or perform Eigen Decompostion of GRM
 * @param {bool} use_eigen
 */
void fastGxE::process_grm(bool use_eigen){
    m_grm_index_vec = {0}; // store start index in sparse GRM for each subgroup
    for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
        MatrixXd tmp_mat = m_grm_mat_group_vec[i];
        long long start_index = m_grm_index_vec.back();
        long long num_element = tmp_mat.rows();
        m_grm_index_vec.push_back(start_index + num_element);
    }

    if(!use_eigen){
        spdlog::info("Prepare sparse GRMs");
        std::vector<Eigen::Triplet<double>> tripletList;
        for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
            MatrixXd tmp_mat = m_grm_mat_group_vec[i];
            long long start_index = m_grm_index_vec[i];
            long long num_element = tmp_mat.rows();
            if(i == 0){
                for(long long m = 0; m < num_element; m++){
                    tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat(m, 0)));
                }
            }else{
                for(long long m = 0; m < num_element; m++){
                    for(long long n = 0; n < num_element; n++)
                        tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat(m, n)));
                }
            }
        }
        m_grm_mat.resize(m_num_id, m_num_id);
        m_grm_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    }else{
        spdlog::info("Eigen Decompostion of GRM");
        std::vector<Eigen::Triplet<double>> tripletList;
        m_grm_eigenvals.resize(m_num_id);
        for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
            MatrixXd tmp_mat = m_grm_mat_group_vec[i];
            long long start_index = m_grm_index_vec[i];
            long long num_element = tmp_mat.rows();
            if(i == 0){
                if(tmp_mat.size() != 0){
                    for(long long m = 0; m < num_element; m++){
                        m_grm_eigenvals(m) = tmp_mat(m, 0);
                        tripletList.push_back(Eigen::Triplet<double>(m, m, 1));
                    }
                }
            }else{
                Eigen::EigenSolver<Eigen::MatrixXd> eigensolver;
                eigensolver.compute(tmp_mat);
                Eigen::VectorXd tmp_eigen_values = eigensolver.eigenvalues().real();
                Eigen::MatrixXd tmp_eigen_vectors = eigensolver.eigenvectors().real();
                for(long long m = 0; m < num_element; m++){
                    m_grm_eigenvals(start_index + m) = tmp_eigen_values(m);
                    for(long long n = 0; n < num_element; n++)
                        tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_eigen_vectors(m, n)));
                }
            }
        }
        m_grm_mat_group_vec.shrink_to_fit();
        m_grm_eigenvecs.resize(m_num_id, m_num_id);
        m_grm_eigenvecs.setFromTriplets(tripletList.begin(), tripletList.end());

        m_y_trans = m_grm_eigenvecs.transpose() * m_y;
        m_xmat_trans = m_grm_eigenvecs.transpose() * m_xmat;
    }
}



/**
 * @Description: Scale the interaction environment covariates and add to xmat
 */
void fastGxE::pre_data_GxE(bool standardize_env, bool phen_correct){
    spdlog::info("Scale the interaction environment covariates");
    // rank of bye_mat
    ColPivHouseholderQR<MatrixXd> qr(this->m_bye_mat);
    if(qr.rank() < m_num_bye){
        spdlog::error("There are dependent columns in the interaction environment covariates!");
        exit(1);
    }

    if(standardize_env){
        // scale to (0, 1)
        VectorXd meanVec = this->m_bye_mat.colwise().mean();
        MatrixXd centered = this->m_bye_mat.rowwise() - meanVec.transpose();
        VectorXd stdVec = (centered.cwiseProduct(centered).colwise().sum().array() / centered.rows()).sqrt();
        this->m_bye_mat = centered.array().rowwise() * (1 / stdVec.transpose().array());
    }
    // std::cout << _bye_mat.transpose() << std::endl << std::endl;
    long long nrows = this->m_xmat.rows();
    long long ncols = this->m_xmat.cols();
    this->m_xmat.conservativeResize(nrows, ncols + this->m_bye_mat.cols());
    this->m_xmat.rightCols(this->m_bye_mat.cols()) = this->m_bye_mat; // interaction covariate
    qr.compute(this->m_xmat);
    if(qr.rank() < m_xmat.cols()){
        spdlog::error("There are dependent columns between interaction environment covariates and other covariates!");
        exit(1);
    }

    if(phen_correct){
        this->m_y = (this->m_y - m_xmat * ((m_xmat.transpose() * m_xmat).inverse() * (m_xmat.transpose() * m_y))).eval();
        this->m_xmat = MatrixXd::Ones(nrows, 1);
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
VectorXd fastGxE::varcom_GxE(VectorXd& init_varcom, int maxiter0, double cc_par0, double cc_gra0, double cc_logL0){
    spdlog::info("Estimate variances...");

    // * set the initial variances
    long long num_cov = 3;
    VectorXd varcom;
    if((init_varcom.size() == num_cov) && (init_varcom.minCoeff() > 0)){
        varcom = init_varcom;
    }else{
        varcom = VectorXd::Ones(num_cov);
    }
    
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "", "", "", "", "");
    std::ostringstream oss;
    oss << varcom.transpose().format(fmt);
    spdlog::info("Initial variances are: {}", oss.str()); 


    spdlog::info("Calculate GxE relationship matrix");
    // MM'/m ⊙ EE'/L
    std::vector<Eigen::Triplet<double>> tripletList;
    vector<MatrixXd> gxe_rm_group_vec(m_grm_mat_group_vec.size());
    for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
        MatrixXd tmp_mat = m_grm_mat_group_vec[i];
        gxe_rm_group_vec[i] = tmp_mat; // Elements index 0 will not be null
        long long start_index = m_grm_index_vec[i];
        long long num_element = tmp_mat.rows();
        MatrixXd bye_mat_part = this->m_bye_mat.middleRows(start_index, num_element);
        if(i == 0){
            if(tmp_mat.size() != 0){
                MatrixXd tmp_mat2 = (bye_mat_part.cwiseProduct(bye_mat_part)).rowwise().sum() / this->m_num_bye;
                tmp_mat2 = (tmp_mat2.cwiseProduct(tmp_mat)).eval();
                gxe_rm_group_vec[i] = tmp_mat2;
                for(long long m = 0; m < num_element; m++){
                    tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat2(m, 0)));
                }
            }
        }else{
            MatrixXd tmp_mat2 = bye_mat_part * bye_mat_part.transpose() / this->m_num_bye;
            tmp_mat2 = (tmp_mat2.cwiseProduct(tmp_mat)).eval();
            gxe_rm_group_vec[i] = tmp_mat2;
            for(long long m = 0; m < num_element; m++){
                for(long long n = 0; n < num_element; n++)
                    tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat2(m, n)));
            }
        }
    }
    tripletList.shrink_to_fit();
    SparseMatrix < double > gxe_rm(m_num_id, m_num_id);
    gxe_rm.setFromTriplets(tripletList.begin(), tripletList.end());
    tripletList.shrink_to_fit();

    spdlog::info("Start the iteration");

    double logL = -1e100, logL_curr = -1e100;
    bool isCC = true;
    for(int iter_count = 0; iter_count < maxiter0; iter_count++){
        spdlog::info("Iteration {}", iter_count + 1);

        // Vi log|V|
        double vmat_logdet = 0;
        std::vector<Eigen::Triplet<double>> tripletList;
        for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
            MatrixXd tmp_grm_mat = m_grm_mat_group_vec[i];
            MatrixXd tmp_gxe_mat = gxe_rm_group_vec[i];
            long long start_index = m_grm_index_vec[i];
            long long num_element = tmp_grm_mat.rows();
            if(i == 0){
                if(tmp_grm_mat.size() != 0){
                    Eigen::ArrayXXd tmp_mat = tmp_grm_mat.array() * varcom(0) + 
                                tmp_gxe_mat.array() * varcom(1) + varcom(2);
                    vmat_logdet += tmp_mat.log().sum();
                    tmp_mat = 1 / tmp_mat.array();
                    for(long long m = 0; m < num_element; m++){
                        tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat(m, 0)));
                    }
                }
            }else{
                MatrixXd tmp_mat = MatrixXd::Identity(num_element, num_element) * varcom(2);
                tmp_mat += tmp_grm_mat * varcom(0) + tmp_gxe_mat * varcom(1);

                Eigen::LDLT<Eigen::MatrixXd> ldlt(tmp_mat);
                Eigen::VectorXd diag = ldlt.vectorD();
                vmat_logdet += diag.array().log().sum();
                tmp_mat = tmp_mat.inverse();
                for(long long m = 0; m < num_element; m++){
                    for(long long n = 0; n < num_element; n++)
                        tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat(m, n)));
                }
            }
        }
        
        //spdlog::info("Allocate");
        SparseMatrix < double > vmat(m_num_id, m_num_id);
        vmat.setFromTriplets(tripletList.begin(), tripletList.end());
        // tripletList.shrink_to_fit();
       // spdlog::info("Free");
        std::vector<Eigen::Triplet<double>>().swap(tripletList);

        // FD, logL, KPy, PKPy
        // Gradient
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
        logL = vmat_logdet + XViX_logdet + (m_y.transpose() * Py).sum();
        spdlog::info("-2logL: {}", logL);
        MatrixXd KPy = m_grm_mat * Py;
        double tr_ViK = (vmat.cwiseProduct(m_grm_mat)).sum();
        double tr_XViXi_XViKViX = (XViXi.cwiseProduct(XViKViX)).sum();
        fd_mat(0) = tr_ViK - tr_XViXi_XViKViX - (Py.transpose() * KPy).sum();
        KPy_vec[0] = KPy;

        double tr_ViKe = (vmat.cwiseProduct(gxe_rm)).sum();
        MatrixXd XViKeViX = ViX.transpose() * gxe_rm * ViX;
        double tr_XViXi_XViKeViX = (XViXi.cwiseProduct(XViKeViX)).sum();
        MatrixXd KePy = gxe_rm * Py;
        fd_mat(1) = tr_ViKe - tr_XViXi_XViKeViX - (Py.transpose() * KePy).sum();
        KPy_vec[1] = KePy;

        fd_mat(num_cov - 1) = vmat.diagonal().sum() - (XViXi.cwiseProduct(ViX.transpose() * ViX)).sum() - (Py.transpose() * Py).sum();
        fd_mat *= -0.5;
        KPy_vec[num_cov - 1] = Py;

        vector <MatrixXd> PKPy_vec(num_cov);
        for(long long m = 0; m < num_cov; m++){
            MatrixXd KPy_tmp = KPy_vec[m];
            PKPy_vec[m] = vmat * KPy_tmp - ViX * (XViXi * (ViX.transpose() * KPy_tmp));
        }

        // AI
        MatrixXd ai_mat(num_cov, num_cov);
        for(long long m = 0; m < num_cov; m++){
            ai_mat(m, m) = (KPy_vec[m].transpose() * PKPy_vec[m]).sum();
            for(long long n = 0; n < m; n++){
                ai_mat(m, n) = ai_mat(n, m) = (KPy_vec[m].transpose() * PKPy_vec[n]).sum();
            }
        }
        ai_mat *= 0.5;

        // EM
        MatrixXd em_mat = MatrixXd::Zero(num_cov, num_cov);
        em_mat.diagonal() = m_num_id / (2 * varcom.array() * varcom.array());
        
        // Update
        VectorXd delta, varcom_update;
        for (int j = 0; j <= 100; j++) {
            double gamma = j * 0.01;
            MatrixXd wemai_mat = (1 - gamma) * ai_mat + gamma * em_mat;
            delta = wemai_mat.inverse() * fd_mat;
            varcom_update = varcom + delta;
            if (varcom_update.minCoeff() > 0) {
                spdlog::info("EM weight value: {}", gamma);
                break;
            }
        }

        std::ostringstream oss;
        oss << varcom_update.transpose().format(fmt);
        spdlog::info("Updated variances: {}", oss.str()); 

        varcom = varcom_update;
        
        double cc_par_val = (delta.cwiseProduct(delta)).sum() / (varcom_update.cwiseProduct(varcom_update)).sum();
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
    
    // Vi
    tripletList.clear();
    logL = 0;
    for(long long i = 0; i < m_grm_mat_group_vec.size(); i++){
        MatrixXd tmp_grm_mat = m_grm_mat_group_vec[i];
        MatrixXd tmp_gxe_mat = gxe_rm_group_vec[i];
        long long start_index = m_grm_index_vec[i];
        long long num_element = tmp_grm_mat.rows();
        if(i == 0){
            if(tmp_grm_mat.size() != 0){
                Eigen::ArrayXXd tmp_mat = tmp_grm_mat.array() * varcom(0) + 
                            tmp_gxe_mat.array() * varcom(1) + varcom(2);
                logL += tmp_mat.log().sum();
                tmp_mat = 1 / tmp_mat.array();
                for(long long m = 0; m < num_element; m++){
                    tripletList.push_back(Eigen::Triplet<double>(m, m, tmp_mat(m, 0)));
                }
            }
        }else{
            MatrixXd tmp_mat = MatrixXd::Identity(num_element, num_element) * varcom(2);
            tmp_mat += tmp_grm_mat * varcom(0) + tmp_gxe_mat * varcom(1);

            Eigen::LDLT<Eigen::MatrixXd> ldlt(tmp_mat);
            Eigen::VectorXd diag = ldlt.vectorD();
            logL += diag.array().log().sum();
            tmp_mat = tmp_mat.inverse();

            for(long long m = 0; m < num_element; m++){
                for(long long n = 0; n < num_element; n++)
                    tripletList.push_back(Eigen::Triplet<double>(start_index + m, start_index + n, tmp_mat(m, n)));
            }
        }
    }
    this->m_V0_logdet = logL;
    
    this->m_Vi0.resize(m_num_id, m_num_id);
    this->m_Vi0.setZero();
    this->m_Vi0.setFromTriplets(tripletList.begin(), tripletList.end());
    tripletList.shrink_to_fit();

    // logL
    MatrixXd ViX = m_Vi0 * m_xmat;
    MatrixXd XViX = m_xmat.transpose() * ViX;

    Eigen::LDLT<Eigen::MatrixXd> ldlt(XViX);
    Eigen::VectorXd diag = ldlt.vectorD();
    double XViX_logdet = diag.array().log().sum();
    logL += XViX_logdet;
    MatrixXd XViXi = XViX.inverse();
    
    MatrixXd Py = this->m_Vi0 * m_y - ViX * (XViXi * (ViX.transpose() * m_y));
    logL += (m_y.transpose() * Py).sum();
    spdlog::info("-2logL in null model: {}", logL);
    m_logL_null = logL;
    
    m_varcom_null = varcom;

    ofstream fout;
    fout.open(m_out_file + ".var");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file: {}.var", m_out_file);
        exit(1);
    }

    fout << varcom << std::endl;

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

    // * set the initial variances
    init = VectorXd::Ones(2);
    VectorXd varcom = init;

    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "", "", "", "", "");
    std::ostringstream oss;
    oss << varcom.transpose().format(fmt);
    spdlog::info("Initial variances are: {}", oss.str()); 

    double logL = -1e100, logL_curr = -1e100;
    bool isCC = true;

    for (int i = 0; i < maxiter0; i++) {

        spdlog::info("Iteration {}", i + 1);
        
        VectorXd vmat = (1.0 / (m_grm_eigenvals.array() * varcom(0) + varcom(1))).matrix();
        double logL_curr = - (vmat.array().log()).sum(); // !reml log|v|
        Eigen::ArrayXXd vxmat_arr = m_xmat_trans.array();
        vxmat_arr.colwise() *= vmat.array();
        MatrixXd vxmat = vxmat_arr.matrix();
        MatrixXd xvxmat = m_xmat_trans.transpose() * vxmat;
        CustomLLT llt_solver;
        llt_solver.compute(xvxmat);
        logL_curr += llt_solver.logDeterminant(); // !reml log|X'VX|
        xvxmat = llt_solver.inverse();
        VectorXd xvy = vxmat.transpose() * m_y_trans;
        VectorXd py = vmat.cwiseProduct(m_y_trans - m_xmat_trans * (xvxmat * xvy));
        logL_curr += (m_y_trans.transpose() * py).sum(); // !reml y'py
        VectorXd dpy = m_grm_eigenvals.cwiseProduct(py);
        xvy = vxmat.transpose() * dpy;
        VectorXd pdpy = vmat.cwiseProduct(dpy - m_xmat_trans * (xvxmat * xvy));
        xvy = vxmat.transpose() * py;
        VectorXd ppy = vmat.cwiseProduct(py - m_xmat_trans * (xvxmat * xvy));
        // *first-order partial derivative
        VectorXd fd_mat(2);
        fd_mat(0) = (vmat.cwiseProduct(m_grm_eigenvals)).sum();
        vxmat_arr.colwise() *= m_grm_eigenvals.array();
        fd_mat(0) -= (xvxmat.cwiseProduct(vxmat_arr.matrix().transpose() * vxmat)).sum();
        fd_mat(0) -= (py.transpose() * dpy).sum();
        fd_mat(0) *= -0.5;
        fd_mat(1) = vmat.sum() - (xvxmat.cwiseProduct(vxmat.transpose() * vxmat)).sum();
        fd_mat(1) -= (py.transpose() * py).sum();
        fd_mat(1) *= -0.5;
        // *AI matrix
        MatrixXd ai_mat(2,2);
        ai_mat(0, 0) = 0.5*dpy.dot(pdpy);
        ai_mat(1, 1) = 0.5*py.dot(ppy);
        ai_mat(0, 1) = 0.5*dpy.dot(ppy);
        ai_mat(1, 0) = ai_mat(0, 1);
        // *EM matrix
        MatrixXd em_mat = MatrixXd::Zero(2, 2);
        em_mat(0, 0) = m_num_id / (2 * varcom(0) * varcom(0));
        em_mat(1, 1) = m_num_id / (2 * varcom(1) * varcom(1));

        VectorXd delta, varcom_update;
        for (int j = 0; j <= 100; j++) {
            double gamma = j * 0.01;
            MatrixXd wemai_mat = (1 - gamma) * ai_mat + gamma * em_mat;
            delta = wemai_mat.inverse() * fd_mat;
            varcom_update = varcom + delta;
            if (varcom_update.minCoeff() > 0) {
                spdlog::info("EM weight value: {}", gamma);
                break;
            }
        }

        std::ostringstream oss;
        oss << varcom_update.transpose().format(fmt);
        spdlog::info("Updated variances: {}", oss.str()); 
        
        varcom = varcom_update;

        double cc_par_val = 1000.0, cc_gra_val = 1000.0;
        isCC = convergence_criteria(logL_curr, logL, cc_logL0,
                    delta, cc_par_val, cc_par0,
                    fd_mat, cc_gra_val, cc_gra0);
        logL = logL_curr;
        if (isCC)
            break;
    }
    
    if (isCC)
        spdlog::info("Variances Converged");
    else
        spdlog::info("Variances Not Converged");
    m_varcom_null = varcom;

    ofstream fout;
    fout.open(m_out_file + ".var");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file: {}.var", m_out_file);
        exit(1);
    }

    fout << varcom << std::endl;
    
    return m_varcom_null;
}



/**
 * @Description: Test SNPs with one random-effect model
 */
void fastGxE::test_main(string bed_file, long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut, 
                int maxiter, double cc_par, double cc_gra){
    spdlog::info("Assocation...");
    GENO GenoA(bed_file);
    vector<long long> id_index_in_bed_vec = GenoA.find_fam_index_pro(m_id_in_data_vec); 
    long long num_id_used = m_id_in_data_vec.size();

    long long num_snp = GenoA.get_num_sid(); // the number of SNP in the plink file.
    long long num_id_in_bed = GenoA.get_num_iid(); // the number of IDs in the plink file.
    spdlog::info("The number of individuals and SNPs in the plink file are: {} {}", num_id_in_bed, num_snp);
    GenoA.check_bed_numChar();

    // *xmat
    this->m_xmat_trans = this->m_grm_eigenvecs.transpose() * this->m_xmat;
    int xmat_df0 = m_xmat_trans.cols();
    int xmat_df1 = xmat_df0 + 1; // with snp

    VectorXd vmat0 = 1.0 / (this->m_grm_eigenvals.array() * this->m_varcom_null(0) + this->m_varcom_null(1));
    // std::cout << "varcom in null: " << _varcom_null.transpose() << std::endl;
    Eigen::ArrayXXd vxmat_arr0 = this->m_xmat_trans.array();
    vxmat_arr0.colwise() *= vmat0.array();
    MatrixXd vxmat0 = vxmat_arr0.matrix();
    MatrixXd xvxmat0 = m_xmat_trans.transpose() * vxmat0;
    MatrixXd xvxmat_inv0 = xvxmat0.inverse();
    MatrixXd xvxmat_inv_xvmat0 = xvxmat_inv0 * vxmat0.transpose();
    
    spdlog::info("Randomly select SNPs and calculate the gamma");
    vector <long long> random_snp_index_vec = GenerateDiffNumber(0, num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr;
    GenoA.read_bed_by_snp_index(bed_file + ".bed", snp_mat_part_random, freq_arr, missing_rate_arr, random_snp_index_vec, id_index_in_bed_vec);
    snp_mat_part_random.rowwise() -= 2*freq_arr.transpose();
    snp_mat_part_random = (m_grm_eigenvecs.transpose() * snp_mat_part_random).eval();

    VectorXd gamma_correct_Vec(num_random_snp);
    VectorXd p_Vec(num_random_snp);
    #pragma omp parallel for schedule(dynamic)
    for(long long k = 0; k < num_random_snp; k++){
        VectorXd isnp_arr = snp_mat_part_random.col(k);
        // !score test
        VectorXd psnp = vmat0.cwiseProduct(isnp_arr - m_xmat_trans * (xvxmat_inv_xvmat0 * isnp_arr));
        double p_score = gsl_cdf_chisq_Q(std::pow((psnp.transpose() * m_y_trans).sum(), 2) / (isnp_arr.transpose() * psnp).sum(), 1);

        double geno_var_isnp = isnp_arr.transpose() * isnp_arr;
        gamma_correct_Vec(k) = (isnp_arr.transpose() * psnp).sum() / geno_var_isnp;
        p_Vec(k) = p_score;
    }
    snp_mat_part_random.resize(0, 0);

    long long num_random_snp_used = 0;
    double gamma_correct = 0;
    ofstream fout;
    fout.open(m_out_file + ".random.res");
    if(!fout.is_open()){
        spdlog::error("Fail to open the file: {}.random.res", m_out_file);
        exit(1);
    }
    fout << "order af missing gamma p" << std::endl;
    for(long long k = 0; k < num_random_snp; k++){
        if (isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
            fout << random_snp_index_vec[k] << " " << freq_arr(k) << " " << 
                    missing_rate_arr(k) << " " << gamma_correct_Vec(k) << " " << p_Vec(k) << std::endl;
            num_random_snp_used++;
            gamma_correct += gamma_correct_Vec[k];
        }
    }
    fout.close();

    if(num_random_snp_used * 1.0 / num_random_snp < 0.1){
        spdlog::error("The number of used random SNPs to calculate gamma is: {}, which is less than 10% of random SNPs", 
                     num_random_snp_used);
        exit(1);
    }
    gamma_correct /= num_random_snp_used;
    spdlog::info("Gamma factor: {}", gamma_correct);
    spdlog::info("The number of used random SNPs is: {}", num_random_snp_used);

    // pymat of raw data (not rotated)
    spdlog::info("Start testing SNPs...");
    fout.open(m_out_file + ".res");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file: {}.res", m_out_file);
        exit(1);
    }
    fout << "order chrom SNP cm base allele1 allele2 af missing beta se p" << std::endl;
    
    VectorXd pyarr0 = this->m_grm_eigenvecs * (vmat0.cwiseProduct(this->m_y_trans) - vxmat0 * (xvxmat_inv0 * (vxmat0.transpose() * this->m_y_trans)));
    double *pyarr0_pt = pyarr0.data();
    long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
    double *snp_mat_part = new double[num_snp_read * num_id_used];
    double* geno_var_arr = new double[num_snp_read];
    double* se_arr = new double[num_snp_read];
    double* prod_arr = new double[num_snp_read];
    double* eff_arr = new double[num_snp_read];
    VectorXd p_arr = VectorXd::Constant(num_snp_read, std::numeric_limits<double>::quiet_NaN());
    
    for(long long i = 0; i < npart_snp; i++){
        long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
        long long start_snp = i*num_snp_read + start_pos;
        if(i == npart_snp - 1 && (end_pos - start_pos) % npart_snp != 0){
            num_snp_read += (end_pos - start_pos) % npart_snp;
            delete[] snp_mat_part;
            delete[] geno_var_arr;
            delete[] se_arr;
            delete[] prod_arr;
            delete[] eff_arr;
            snp_mat_part = new double[num_snp_read * num_id_used];
            geno_var_arr = new double[num_snp_read];
            se_arr = new double[num_snp_read];
            prod_arr = new double[num_snp_read];
            eff_arr = new double[num_snp_read];
            p_arr = VectorXd::Constant(num_snp_read, std::numeric_limits<double>::quiet_NaN());
        }

        spdlog::info("Part {}/{}: Start SNP, {}; the number of read SNP, {}", 
             i + 1, npart_snp, start_snp, num_snp_read);
        
        // Read
        GenoA.read_bed_omp_spMVLMM(bed_file + ".bed", snp_mat_part, freq_arr, missing_rate_arr, 
                start_snp, num_snp_read, id_index_in_bed_vec);
        
        // Test
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < num_snp_read; ++j) {
            geno_var_arr[j] = cblas_ddot(num_id_used, snp_mat_part + j * num_id_used, 1, snp_mat_part + j * num_id_used, 1);
        }

        // VectorXd se_arr = 1 / (gamma_correct * geno_var_arr).array();
        // Compute the inverse square root of each element
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < num_snp_read; ++j) {
            se_arr[j] = 1.0 / (gamma_correct * geno_var_arr[j]);
        }

        // VectorXd eff_arr = se_arr.cwiseProduct(snp_mat_part.transpose() * pyarr0);
        // Compute the product of the transpose of snp_mat_part and pyarr0"
        cblas_dgemv(CblasRowMajor, CblasNoTrans, num_snp_read, num_id_used, 1.0, snp_mat_part, 
                num_id_used, pyarr0_pt, 1, 0.0, prod_arr, 1);
        
        // Compute the element-wise product of se_arr and prod_arr
        vdMul(num_snp_read, se_arr, prod_arr, eff_arr);

        // se_arr = se_arr.cwiseSqrt()
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < num_snp_read; ++j) {
            se_arr[j] = sqrt(se_arr[j]);
        }
        
        // spdlog::info("Give P");
        #pragma omp parallel for schedule(dynamic)
        for (long long k = 0; k < num_snp_read; k++) {
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
        vector<string> snp_info_vec = GenoA.snp_anno(start_snp, num_snp_read);
        fastGxE::output(fout, snp_info_vec, freq_arr, missing_rate_arr,  
            eff_arr, se_arr, p_arr, start_snp, num_snp_read, p_cut);
    }
    fout.close();
    // free memory
    delete[] snp_mat_part;
    delete[] geno_var_arr;
    delete[] se_arr;
    delete[] prod_arr;
    delete[] eff_arr;
}


void fastGxE::output(ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr, 
                double* eff_arr, double* se_arr, VectorXd& p_arr,
                long long start_snp, long long num_snp_read, double p_cut){
    for(long long i = 0; i < num_snp_read; i++){
        fout << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " "  << missing_rate_arr(i) << " "
                << " " << eff_arr[i] << " " << se_arr[i] 
                << " "<< p_arr(i) << std::endl;
    }
}

void fastGxE::reset_output_prefix(string out_file, long long first, long long second){
    m_out_file = out_file + "." + std::to_string(first) + "_"  + std::to_string(second);
}

bool fastGxE::isValidSnp(double af, double missing_rate, double maf_cut, double missing_rate_cut) {
    return af > maf_cut && af < 1 - maf_cut && missing_rate < missing_rate_cut;
}



/**
 * @Description: Assocation of GxE using an improved structLMM model
 */
void fastGxE::test_GxE_multi(string bed_file, long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut){
    spdlog::info("Assocation...");
    GENO GenoA(bed_file);
    vector<long long> id_index_in_bed_vec = GenoA.find_fam_index_pro(m_id_in_data_vec);
    long long num_id_used = id_index_in_bed_vec.size();

    long long num_snp = GenoA.get_num_sid(); // the number of SNP in the plink file.
    long long num_id_in_bed = GenoA.get_num_iid(); // the number of IDs in the plink file.
    spdlog::info("The number of individuals and SNPs in the plink file are: {} {}", num_id_in_bed, num_snp);
    GenoA.check_bed_numChar();

    VectorXd Vi0y = this->m_Vi0 * m_y;
    double *bye_mat_pt = this->m_bye_mat.data();
    double *Vi0y_pt = Vi0y.data();

    spdlog::info("Randomly select SNPs and calculate the gamma");
    vector <long long> random_snp_index_vec = GenerateDiffNumber(0, num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr;
    GenoA.read_bed_by_snp_index(bed_file + ".bed", snp_mat_part_random, freq_arr, missing_rate_arr, random_snp_index_vec, id_index_in_bed_vec);
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
    for(long long i = 0; i < num_random_snp; i++){
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

    long long num_random_snp_used = 0;
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
    for(long long i = 0; i < num_random_snp; i++){
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
    long long _omp_max_threads = omp_get_max_threads();

    double **EP0y_part = new double*[_omp_max_threads];
    for(int i = 0; i < _omp_max_threads; i++){
        EP0y_part[i] = new double[num_id_used];
    }

    double **WEP0y_part = new double*[_omp_max_threads];
    for(int i = 0; i < _omp_max_threads; i++){
        WEP0y_part[i] = new double[m_num_bye];
    }

    long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
    double *snp_mat_part = new double[num_snp_read * num_id_used];
    beta_main_Vec.resize(num_snp_read);
    se_main_Vec.resize(num_snp_read);
    p_main_Vec.resize(num_snp_read);
    VectorXd score_Vec(num_snp_read);
    VectorXd p_gxe_Vec(num_snp_read);
    
    for(long long ipart = 0; ipart < npart_snp; ipart++){
        long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
        long long start_snp = ipart*num_snp_read + start_pos;

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
        
        GenoA.read_bed_omp_spMVLMM(bed_file + ".bed", snp_mat_part, freq_arr, missing_rate_arr, 
                start_snp, num_snp_read, id_index_in_bed_vec);

        #pragma omp parallel for schedule(dynamic)
        for(long long k = 0; k < num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                long long thread_id = omp_get_thread_num();

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

void fastGxE::output_GxE_multi(std::ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr,
            VectorXd& beta_main_Vec, VectorXd& se_main_Vec, VectorXd& p_main_Vec, 
            VectorXd& score_Vec, VectorXd& p_Vec, long long start_snp, long long num_snp_read){
    for(long long i = 0; i < num_snp_read; i++){
        fout << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) << " "
            << beta_main_Vec(i) << " " << se_main_Vec(i) << " " << p_main_Vec(i) << " "
            << score_Vec(i) << " " << p_Vec(i) << std::endl;
    }
}

void fastGxE::test_GxE(string bed_file, 
                long long start_pos, long long end_pos, int npart_snp,
                int speed, long long num_random_snp,
                double p_approx_cut, double p_cut, double maf_cut, double missing_rate_cut){
    spdlog::info("GxE analysis for individual environments");
    GENO GenoA(bed_file);
    vector<long long> id_index_in_bed_vec = GenoA.find_fam_index_pro(m_id_in_data_vec); 
    long long num_id_used = id_index_in_bed_vec.size();

    long long num_snp = GenoA.get_num_sid(); // the number of SNP in the plink file.
    long long num_id_in_bed = GenoA.get_num_iid(); // the number of IDs in the plink file.
    spdlog::info("The number of individuals and SNPs in the plink file are: {} {}", num_id_in_bed, num_snp);
    spdlog::info("The number of used iids: {}", num_id_used);
    GenoA.check_bed_numChar();
    VectorXd Vi0y = this->m_Vi0 * m_y;
    
    spdlog::info("Randomly select SNPs and calculate the gamma");
    vector <long long> random_snp_index_vec = GenerateDiffNumber(0, num_snp, num_random_snp);
    MatrixXd snp_mat_part_random;
    VectorXd freq_arr, missing_rate_arr;
    GenoA.read_bed_by_snp_index(bed_file + ".bed", snp_mat_part_random, freq_arr, missing_rate_arr, random_snp_index_vec, id_index_in_bed_vec);
    snp_mat_part_random.rowwise() -= 2*freq_arr.transpose();
    VectorXd snp_squaredNorm_Vec = snp_mat_part_random.colwise().squaredNorm();

    spdlog::info("main");
    VectorXd p_main_Vec = VectorXd::Ones(num_random_snp) * NAN;
    VectorXd beta_main_var_Vec = 1 / ((snp_mat_part_random).cwiseProduct(this->m_Vi0 * snp_mat_part_random)).colwise().sum().array();
    VectorXd se_main_Vec = beta_main_var_Vec.array().sqrt();
    VectorXd beta_main_Vec = beta_main_var_Vec.cwiseProduct(snp_mat_part_random.transpose() * Vi0y);
    VectorXd gamma_main_Vec = 1 / (beta_main_var_Vec.array() * snp_squaredNorm_Vec.array());  // gamma = x'V(-1)x / x'x

    long long num_random_snp_used = 0;
    double gamma_main = 0;
    for(long long i = 0; i < num_random_snp; i++){
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
    gamma_main = gamma_main / num_random_snp_used;

    if(num_random_snp_used * 1.0 / num_random_snp < 0.5){
        spdlog::error("Used random SNPs: {} (<50%)", num_random_snp_used);
        exit(1);
    }

    ofstream fout;
    fout.open(m_out_file + ".main.random.res");
    if(!fout.is_open()){
        spdlog::error("Fail to open output file: {}.main.random.res", m_out_file);
        exit(1);
    }
    
    fout << "order af missing gamma beta se p" << std::endl;
    for(long long i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i) << " "
              << gamma_main_Vec[i] << " " << beta_main_Vec[i] << " " << se_main_Vec[i] << " " << p_main_Vec[i]
               << std::endl;
    }
    fout.close();
    gamma_main_Vec.resize(0); beta_main_Vec.resize(0); se_main_Vec.resize(0); p_main_Vec.resize(0); beta_main_var_Vec.resize(0);


    spdlog::info("GxE without main");
    VectorXd gamma_GxEnoMain_Vec = VectorXd::Zero(this->m_num_bye);
    MatrixXd beta_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    MatrixXd se_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    MatrixXd p_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    MatrixXd gamma_GxEnoMain_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    for(long long j = 0; j < this->m_num_bye; j++){
        MatrixXd snp_mat_part_random_xE = snp_mat_part_random.array().colwise() * this->m_bye_mat.col(j).array();
        VectorXd beta_GxEnoMain_var_Vec = 1 / ((snp_mat_part_random_xE).cwiseProduct(this->m_Vi0 * snp_mat_part_random_xE)).colwise().sum().array();
        VectorXd se_GxEnoMain_Vec = beta_GxEnoMain_var_Vec.array().sqrt();
        VectorXd beta_GxEnoMain_Vec = beta_GxEnoMain_var_Vec.cwiseProduct(snp_mat_part_random_xE.transpose() * Vi0y);
        VectorXd gamma_GxEnoMain_j = 1 / (beta_GxEnoMain_var_Vec.array() * snp_squaredNorm_Vec.array());
        VectorXd p_GxEnoMain_Vec = VectorXd::Ones(num_random_snp) * NAN;
        for(long long i = 0; i < num_random_snp; i++){
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

    fout.open(m_out_file + ".GxEnoMain.random.res");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file: " + m_out_file + ".GxEnoMain.random.res");
        exit(1);
    }
    fout << "order af missing";
    for(long long i = 0; i < m_num_bye; i++){
        fout << " gamma" << i + 1 << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << std::endl;

    for(long long i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        for(long long j = 0; j < m_num_bye; j++){
            fout <<  " " << gamma_GxEnoMain_Mat(i, j) << " " << beta_GxEnoMain_Mat(i, j) << " " 
            << se_GxEnoMain_Mat(i, j) << " " << p_GxEnoMain_Mat(i, j);
        }
        fout << std::endl;
    }
    fout.close();
    gamma_GxEnoMain_Mat.resize(0, 0); beta_GxEnoMain_Mat.resize(0, 0); se_GxEnoMain_Mat.resize(0, 0); p_GxEnoMain_Mat.resize(0, 0);

    spdlog::info("GxE with main");
    MatrixXd beta_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    MatrixXd se_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    MatrixXd p_Mat = MatrixXd::Ones(num_random_snp, this->m_num_bye) * NAN;
    vector<vector<MatrixXd>> gamma_Mat_2Dvec;
    for(long long j = 0; j < this->m_num_bye; j++){
        vector<MatrixXd> gamma_Mat_vec;
        for(long long i = 0; i < num_random_snp; i++){
            gamma_Mat_vec.push_back(MatrixXd::Zero(2, 2));
        }
        gamma_Mat_2Dvec.push_back(gamma_Mat_vec);
    }

    #pragma omp parallel for schedule(dynamic)
    for(long long i = 0; i < num_random_snp; i++){
        if(isValidSnp(freq_arr(i), missing_rate_arr(i), maf_cut, missing_rate_cut)){
            VectorXd snpi = snp_mat_part_random.col(i);
            double snp_norm2 = snp_squaredNorm_Vec(i);
            MatrixXd snpiME = this->m_bye_mat.array().colwise() * snpi.array();
            MatrixXd xmati = MatrixXd::Zero(num_id_used, 2);
            xmati.col(0) = snpi;
            VectorXd beta_bye_Vec(m_num_bye);
            VectorXd p_bye_Vec(m_num_bye);
            for(long long j = 0; j < this->m_num_bye; j++){
                xmati.col(1) = snpiME.col(j);
                MatrixXd XVi0X = xmati.transpose() * this->m_Vi0 * xmati;
                MatrixXd XVi0Xi = XVi0X.inverse();
                VectorXd beta_Vec = XVi0Xi * (xmati.transpose() * Vi0y);
                double beta_var = XVi0Xi(1, 1);
                double beta = beta_Vec(1);
                double p = gsl_cdf_chisq_Q(beta * beta / beta_var, 1);
                beta_Mat(i, j) = beta;
                beta_bye_Vec(j) = beta;
                se_Mat(i, j) = sqrt(beta_var);
                p_Mat(i, j) = p;
                p_bye_Vec(j) = p;
                gamma_Mat_2Dvec[j][i] = XVi0X / snp_norm2;
            }
        }
    }

    vector<MatrixXd> gamma_Mat_vec;
    for(long long j = 0; j < this->m_num_bye; j++){
        MatrixXd gamma_Mat = MatrixXd::Zero(2, 2);
        for(long long i = 0; i < num_random_snp; i++){
            gamma_Mat += gamma_Mat_2Dvec[j][i];
        }
        gamma_Mat_vec.push_back(gamma_Mat / num_random_snp_used);
    }

    
    fout.open(m_out_file + ".GxE.random.res");
    if(!fout.is_open()){
        spdlog::error("Fail to open output file: {}.GxE.random.res", m_out_file);
        exit(1);
    }
    fout << "order af missing";
    for(long long i = 0; i < m_num_bye; i++){
        fout << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << std::endl;
    for(long long i = 0; i < num_random_snp; i++){
        fout << random_snp_index_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        for(long long j = 0; j < m_num_bye; j++){
            fout << " " << beta_Mat(i, j) << " " << se_Mat(i, j) << " " << p_Mat(i, j);
        }
        fout << std::endl;
    }
    fout.close();

    // the correlation among the environments
    MatrixXd corrE = this->m_bye_mat.transpose() * this->m_bye_mat / this->m_bye_mat.rows();
    corrE.diagonal() += Eigen::VectorXd::Constant(corrE.rows(), 0.001);

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(corrE);
    if (eigensolver.info() != Eigen::Success) {
        spdlog::error("Fail to eigen decomposition of GRM");
        exit(1);
    }
    
    VectorXd corrE_eigenvalues = eigensolver.eigenvalues();
    
    spdlog::info("Genome-wide association");
    fout.open(m_out_file + ".res");
    if(!fout.is_open()){
        spdlog::error("Fail to open output file: {}.res", m_out_file);
        exit(1);
    }
    fout << "order chrom SNP cm base allele1 allele2 af missing beta_main se_main p_main";
    for(long long i = 0; i < m_num_bye; i++){
        fout << " beta" << i + 1 << " se" << i + 1 << " p" << i + 1;
    }
    fout << " p_single p_multi p_gxe" << std::endl;

    long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
    double *snp_mat_part = new double[num_snp_read * num_id_used];
    VectorXd geno_var_Vec = VectorXd::Zero(num_snp_read);
    beta_Mat.resize(num_snp_read, this->m_num_bye + 1);
    se_Mat.resize(num_snp_read, this->m_num_bye + 1);
    p_Mat.resize(num_snp_read, m_num_bye + 1);
    MatrixXd p_combined_Mat = MatrixXd::Zero(num_snp_read, 3);
    VectorXd EVy = VectorXd::Zero(num_id_used);
    VectorXd prod_Vec = VectorXd::Zero(num_snp_read);
    
    for(long long ipart = 0; ipart < npart_snp; ipart++){
        long long num_snp_read = (long long)((end_pos - start_pos)/npart_snp);
        long long start_snp = ipart*num_snp_read + start_pos;
        if(ipart == npart_snp - 1 && (end_pos - start_pos) % npart_snp != 0){
            num_snp_read = num_snp_read + (end_pos - start_pos) % npart_snp;
            delete [] snp_mat_part; snp_mat_part = new double[num_snp_read * num_id_used];
            geno_var_Vec.resize(num_snp_read);
            beta_Mat.resize(num_snp_read, this->m_num_bye + 1);
            se_Mat.resize(num_snp_read, this->m_num_bye + 1);
            p_Mat.resize(num_snp_read, m_num_bye + 1);
            p_combined_Mat.resize(num_snp_read, 3);
            prod_Vec.resize(num_snp_read);
        }
        
        spdlog::info("Part {}/{}: Start SNP {}, read SNPs {}", ipart + 1, npart_snp, start_snp, num_snp_read);

        GenoA.read_bed_omp_spMVLMM(bed_file + ".bed", snp_mat_part, freq_arr, missing_rate_arr, start_snp, num_snp_read, id_index_in_bed_vec);

        // norm2
        #pragma omp parallel for schedule(dynamic)
        for (long long j = 0; j < num_snp_read; ++j) {
            geno_var_Vec(j) = cblas_ddot(num_id_used, snp_mat_part + j * num_id_used, 1, snp_mat_part + j * num_id_used, 1);
        }

        // main
        se_Mat.col(0) = 1.0 / (gamma_main * geno_var_Vec.array());

        cblas_dgemv(CblasRowMajor, CblasNoTrans, num_snp_read, num_id_used, 1.0, snp_mat_part, 
                num_id_used, Vi0y.data(), 1, 0.0, prod_Vec.data(), 1);
        vdMul(num_snp_read, se_Mat.col(0).data(), prod_Vec.data(), beta_Mat.col(0).data());
        
        // GxE
        for(long long k = 0; k < this->m_num_bye; k++){
            se_Mat.col(k+1) = 1.0 / (gamma_GxEnoMain_Vec(k) * geno_var_Vec.array());
            vdMul(num_id_used, this->m_bye_mat.col(k).data(), Vi0y.data(), EVy.data());
            cblas_dgemv(CblasRowMajor, CblasNoTrans, num_snp_read, num_id_used, 1.0, snp_mat_part, 
                num_id_used, EVy.data(), 1, 0.0, prod_Vec.data(), 1);
            vdMul(num_snp_read, se_Mat.col(k+1).data(), prod_Vec.data(), beta_Mat.col(k+1).data());
        }
        se_Mat = se_Mat.array().sqrt();
        MatrixXd z_mat = beta_Mat.array() / se_Mat.array();

        #pragma omp parallel for schedule(dynamic)
        for(long long k = 0; k < num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                VectorXd p_Vec = VectorXd::Ones(m_num_bye + 1) * NAN;

                for(long long m = 0; m < m_num_bye + 1; m++){
                    double z = z_mat(k, m);
                    double chi2_val = z * z;
                    p_Vec(m) = gsl_cdf_chisq_Q(chi2_val, 1);
                }

                if(p_Vec.minCoeff() < p_approx_cut){
                    Eigen::Map<Eigen::VectorXd> snpk(snp_mat_part + k * num_id_used, num_id_used);
                    double snpk_norm2 = geno_var_Vec(k);
                    MatrixXd snpkME = this->m_bye_mat.array().colwise() * snpk.array();
                    MatrixXd xmatk = MatrixXd::Zero(num_id_used, 2);
                    xmatk.col(0) = snpk;

                    for(long long m = 0; m < m_num_bye; m++){
                        if(p_Vec(m+1) < p_approx_cut){
                            xmatk.col(1) = snpkME.col(m);
                            if(speed == 0){
                                // exact test
                                MatrixXd XVi0X = xmatk.transpose() * this->m_Vi0 * xmatk;
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

                }

                p_Mat.row(k) = p_Vec;
                Eigen::VectorXd result = combine_pvalues_cauchy(p_Vec.tail(this->m_num_bye));
                p_combined_Mat(k, 0) = result(1);

            }else{
                beta_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                se_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                p_Mat.row(k) = VectorXd::Ones(m_num_bye+1) * NAN;
                p_combined_Mat(k, 0) = NAN;
            }
        }

        VectorXd score_Vec = z_mat.rightCols(this->m_num_bye).array().square().matrix().rowwise().sum();
        
        #pragma omp parallel for schedule(dynamic)
        for(long long k = 0; k < num_snp_read; k++){
            if(isValidSnp(freq_arr(k), missing_rate_arr(k), maf_cut, missing_rate_cut)){
                double p1 = saddle(score_Vec(k), corrE_eigenvalues);
                p_combined_Mat(k, 1) = p1;

                VectorXd p_Vec(2);
                p_Vec(0) = p_combined_Mat(k, 0);
                p_Vec(1) = p1;
                p_combined_Mat(k, 2) = combine_pvalues_cauchy(p_Vec)(1);
            }else{
                p_combined_Mat(k, 1) = p_combined_Mat(k, 2) = NAN;
            }
        }

        vector <std::string> snp_info_vec = GenoA.snp_anno(start_snp, num_snp_read);
        fastGxE::output_GxE(fout, snp_info_vec, freq_arr, missing_rate_arr,
            beta_Mat, se_Mat, p_Mat, p_combined_Mat,
            start_snp, num_snp_read);
    }
    delete [] snp_mat_part;
    fout.close();
}

void fastGxE::output_GxE(ofstream& fout, vector<string> snp_info_vec, VectorXd& freq_arr, VectorXd& missing_rate_arr,
            MatrixXd& beta_mat, MatrixXd& se_mat, MatrixXd& p_mat, MatrixXd& p_combined_Mat,
            long long start_snp, long long num_snp_read){
    
    std::stringstream ss;
    ss.precision(8);

    for (long long i = 0; i < num_snp_read; i++) {
        ss << start_snp + i << " " << snp_info_vec[i] << " " << freq_arr(i) << " " << missing_rate_arr(i);
        
        for (long long m = 0; m < beta_mat.cols(); m++) {
            ss << " " << beta_mat(i, m) << " " << se_mat(i, m) << " " << p_mat(i, m);
        }

        for (long long m = 0; m < p_combined_Mat.cols(); m++) {
            ss << " " << p_combined_Mat(i, m);
        }

        ss << std::endl;
    }

    fout << ss.str();
}


int fastGxE::run(int argc, char **argv) {
    CLI::App app{"fastGxE - Scalable and fast multivariate GxE analysis"};

    app.description(R"(
    Quick Start:

      Test GxE:
        fastgxe --test-gxe --grm grm_file --bfile bed_file --data data_file --trait BMI --env-int age smoking:alcohol --out test_gxe

      Test SNP Main Effects:
        fastgxe --test-main --grm grm_file --bfile bed_file --data data_file --trait BMI --out test_main
    )");

    int threads = 10;
    bool test_main = false;
    bool test_gxe = false;
    bool phen_correct = false;
    bool standardize_env = true;
    
    std::string data_file, agrm_file, bed_file, out_file;
    std::vector<std::string> trait_vec, covariate_vec, class_vec, bye_vec, snp_range_vec;
    std::vector<std::string> missing_in_data_vec = {"NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan", "<NA>", "<na>", "N/A", "n/a"};
    std::vector<std::string> missing_in_geno_vec = missing_in_data_vec;
    std::vector<long long> split_task_vec;  // Default: {1, 1} (1 total part, executing part 1)

    long long num_random_snp = 2000;
    int npart_snp = 20, speed = 2, maxiter = 100;
    int continue_varcom = 0;
    double p_approx_cut = 1.0e-3, p_cut = 2, maf_cut = 0.01, missing_rate_cut = 0.05;
    double cc_par = 1.0e-7, cc_gra = 1.0e-6, cc_logL = 5.0e-5;

    // Register CLI11 options
    app.add_option("-p,--threads", threads, 
    "Number of threads to use (default: 10).")
    ->default_val(10);

    // Define both flags before referencing them in 'excludes()'
    auto* test_main_flag = app.add_flag("--test-main", test_main, 
        "Enable testing of SNP main effects (mutually exclusive with --test-gxe).");

    auto* test_gxe_flag = app.add_flag("--test-gxe", test_gxe, 
        "Enable testing of SNP-environment interactions (mutually exclusive with --test-main).");

    // Now use 'excludes()' safely
    test_main_flag->excludes(test_gxe_flag);
    test_gxe_flag->excludes(test_main_flag);

    app.add_option("--data", data_file, 
        "Path to input data file (required).")
        ->required();

    app.add_option("--trait", trait_vec, 
        "Trait to analyze (required, expects 1 value).")
        ->expected(1)->required();

    app.add_option("--env-int", bye_vec, 
        "List of interacting environmental covariates (required).\n"
        "  - Supports multiple values (e.g., --env-int age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --env-int age:BMI).\n"
        "  - Expands to include all covariates in the specified range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1)->required();
    
    app.add_option("--out", out_file, 
        "Path to output file (required).")
        ->required();

    app.add_option("--grm", agrm_file, 
        "Path to genetic relationship matrix (GRM) file (required).")
        ->required();

    app.add_option("--bfile", bed_file, 
        "Path to binary PLINK BED file (optional).");

    app.add_option("--covar", covariate_vec, 
        "List of covariates for the analysis (required).\n"
        "  - Supports multiple values (e.g., --covar age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --covar age:BMI).\n"
        "  - Expands to include all covariates in the range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);

    app.add_option("--class", class_vec, 
        "List of categorical class variables (required).\n"
        "  - Supports multiple values (e.g., --class gender ethnicity region).\n"
        "  - Use ':' to specify a range (e.g., --class gender:region).\n"
        "  - Expands to include all variables in the range.\n"
        "  - Example order: {gender, ethnicity, region, education, income}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1);

    app.add_flag("--no-standardize-env", [&standardize_env](int count) {
        if (count > 0) standardize_env = false;  // Flip standardize_env to false when flag is used
    }, "Disable standardization of interacting environmental covariates.\n"
       "  - Standardization is ENABLED by default (mean 0, std 1).\n"
       "  - Use --no-standardize-env to turn it OFF.");

    app.add_option("--missing-data", missing_in_data_vec, 
        "List of missing value indicators for phenotype/covariates.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-data . -999 \"?\").")
        ->expected(-1);

    app.add_option("--missing-geno", missing_in_geno_vec, 
        "List of missing value indicators for genotype data.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-geno 0 ./ .).")
        ->expected(-1);

    app.add_option("--maxiter", maxiter, 
        "Maximum number of optimization iterations (default: 200).")
        ->default_val(200);

    app.add_option("--cc-par", cc_par, 
        "Convergence threshold for parameter updates (default: 1e-7).")
        ->default_val(1e-7);

    app.add_option("--cc-gra", cc_gra, 
        "Convergence threshold for gradient norm (default: 1e-6).")
        ->default_val(1e-6);

    app.add_option("--cc-logL", cc_logL, 
        "Convergence threshold for log-likelihood change (default: 5e-5).")
        ->default_val(5e-5);

    app.add_option("--snp-range", snp_range_vec, 
        "Specify start and end SNPs (expects 2 values).")
        ->expected(2);

    app.add_option("--npart-SNP", npart_snp, 
        "Number of SNP partitions (default: 1000).")
        ->default_val(1000);

    app.add_option("--speed", speed, 
        "Computation speed level (default: 2, higher is faster).")
        ->default_val(2);

    app.add_option("--p-approx-cut", p_approx_cut, 
        "P-value approximation threshold (default: 1e-3).")
        ->default_val(1e-3);

    app.add_option("--p-cut", p_cut, 
        "P-value significance cutoff for output (default: 2).")
        ->default_val(2);

    app.add_option("--maf", maf_cut, 
        "Minor allele frequency (MAF) cutoff (default: 0.01).")
        ->default_val(0.01);

    app.add_option("--missing-rate", missing_rate_cut, 
        "Missing genotype rate cutoff (default: 0.05).")
        ->default_val(0.05);

    app.add_option("--num-random-snp", num_random_snp, 
        "Number of randomly selected SNPs (default: 2000).")
        ->default_val(2000);

    app.add_option("--split-task", split_task_vec, 
        "Partition the task for parallel execution (expects 2 values: total parts, current part).")
        ->expected(2);


    // Parse command-line arguments
    CLI11_PARSE(app, argc, argv);

    // Log parsed arguments
    spdlog::info("=== Parsed Arguments ===");
    spdlog::info("Threads: {}", threads);
    spdlog::info("Input Data File: {}", data_file);
    spdlog::info("Output File: {}", out_file);
    spdlog::info("GRM File: {}", agrm_file);
    spdlog::info("BED File: {}", bed_file.empty() ? "Not Provided" : bed_file);
    spdlog::info("Trait: {}", trait_vec[0]);
    spdlog::info("Interacting environmental covariates: {}", join_string(bye_vec, ", "));
    spdlog::info("Standardization of Interacting environments: {}", standardize_env ? "ENABLED" : "DISABLED");
    
    
    spdlog::info("Covariates: {}", covariate_vec.empty() ? "None" : join_string(covariate_vec, ", "));
    spdlog::info("Class Variables: {}", class_vec.empty() ? "None" : join_string(class_vec, ", "));
    spdlog::info("Missing Data Indicators: {}", join_string(missing_in_data_vec, ", "));

    spdlog::info("Max Iterations: {}", maxiter);
    spdlog::info("CC Par: {}, CC Gra: {}, CC LogL: {}", cc_par, cc_gra, cc_logL);

    if (snp_range_vec.size() == 2) {
        spdlog::info("SNP Range: {} to {}", snp_range_vec[0], snp_range_vec[1]);
    }

    spdlog::info("SNP Partitions: {}", npart_snp);
    spdlog::info("P-value Approx Cutoff: {}", p_approx_cut);
    spdlog::info("P-value Output Cutoff: {}", p_cut);
    spdlog::info("Minor Allele Frequency Cutoff: {}", maf_cut);
    spdlog::info("Missing Rate Cutoff: {}", missing_rate_cut);
    spdlog::info("Number of Random SNPs: {}", num_random_snp);
    
    if(!split_task_vec.empty()) spdlog::info("Task Partitioning: {} parts, running part {}", split_task_vec[0], split_task_vec[1]);
    
    spdlog::info("========================");
    
    // Set MKL and OpenMP threads
    if (threads <= 0) threads = 10;
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);

    // Obtain head row
    std::ifstream fin(data_file);
    if (!fin.is_open()) {
        spdlog::error("Failed to open data file: {}", data_file);
        exit(1);  // Ensure early exit to prevent further errors
    }

    std::string line;
    if (std::getline(fin, line)) {
        process_line(line);
        if (line.empty()) {
            spdlog::error("Head line is empty or starts with #");
            exit(1);
        }
    } else {
        spdlog::error("Failed to read the head line from file: {}", data_file);
        exit(1);
    }

    std::vector<std::string> head_vec = split_string(line);
    fin.close();


    if (test_main + test_gxe != 1) {
        spdlog::error("Exactly one of --test-main or --test-gxe must be specified.");
        exit(1);
    }

    // Trait index
    
    std::vector<std::string> strNoFound_vec;
    std::vector<long long> trait_index_vec = find_index(head_vec, trait_vec, strNoFound_vec);

    if (!strNoFound_vec.empty()) {
        spdlog::error("Trait names not found in the header: {}", join_string(strNoFound_vec));
        exit(1);
    }
    strNoFound_vec.clear();

    // covariate index
    covariate_vec = expand_variable_ranges(covariate_vec, head_vec);
    if(!covariate_vec.empty()){
        spdlog::info("Number of covariates: {}, Covariates: {}", 
             covariate_vec.size(), join_string(covariate_vec, ", "));
    }
    vector <long long> covariate_index_vec = find_index(head_vec, covariate_vec, strNoFound_vec);
    strNoFound_vec.clear();

    // class index
    class_vec = expand_variable_ranges(class_vec, head_vec);
    if(!class_vec.empty()){
        spdlog::info("Number of class variables: {}, class variables: {}", 
             class_vec.size(), join_string(class_vec, ", "));
    }
    vector <long long> class_index_vec = find_index(head_vec, class_vec, strNoFound_vec);
    strNoFound_vec.clear();


    // interacting environment
    bye_vec = expand_variable_ranges(bye_vec, head_vec);
    spdlog::info("Number of interacting environmental covariates: {}, Interacting environmental covariates: {}", 
             bye_vec.size(), join_string(bye_vec, ", "));
    vector <long long> bye_index_vec = find_index(head_vec, bye_vec, strNoFound_vec);
    if(strNoFound_vec.size() != 0){
        spdlog::error("Interacting environments not found in the header: {}", join_string(strNoFound_vec));
        exit(1);
    }
    strNoFound_vec.clear();
    
    if (!bye_index_vec.empty()) {
        std::vector<long long> tmp1 = find_index(covariate_vec, bye_vec, strNoFound_vec);
        std::vector<long long> tmp2 = find_index(class_vec, bye_vec, strNoFound_vec);
        
        if (!tmp1.empty() || !tmp2.empty()) {
            spdlog::error("Interacting environments should not be included in --covar or --class. They are treated as covariates by default.");
            exit(1);
        }
    }
    strNoFound_vec.clear();



    // prepare data
    this->pre_data(out_file, data_file, agrm_file, covariate_index_vec, class_index_vec, 
            bye_index_vec, trait_index_vec, missing_in_data_vec);
    
    VectorXd init_varcom;
    

    this->process_grm(test_main);
    Eigen::VectorXd varcom;
    if(test_main){
        varcom = this->varcom_main(init_varcom, maxiter, cc_gra, cc_gra, cc_logL);
    }else if(test_gxe){
        phen_correct = true; // must be true, save resources
        this->pre_data_GxE(standardize_env, phen_correct);
        varcom = this->varcom_GxE(init_varcom, maxiter, cc_par, cc_gra, cc_logL);
    }
    
    
    if(bed_file.empty()){
        spdlog::info("If you want to perform GWAS, please give --bfile.");
        return 0;
    }
    
    vector<long long> tmp_vec = this->get_start_end_pos(bed_file, snp_range_vec, split_task_vec);
    long long start_pos = tmp_vec[0], end_pos = tmp_vec[1];
    if(end_pos - start_pos <= npart_snp) npart_snp = 1;

    if(test_main){
        this->test_main(bed_file, start_pos, end_pos, npart_snp,
                speed, num_random_snp, p_approx_cut, p_cut, maf_cut, missing_rate_cut, 
                maxiter, cc_par, cc_gra);
    }else if(test_gxe){
        this->test_GxE(bed_file, start_pos, end_pos, npart_snp,
                                 speed, num_random_snp, p_approx_cut, p_cut, 
                                 maf_cut, missing_rate_cut);
    }
    return 0;
}


vector<long long> fastGxE::get_start_end_pos(string bed_file, 
            vector<string> snp_range_vec, vector<long long> split_task_vec){
    GENO GenoA(bed_file);
    long long num_snp = GenoA.get_num_sid();
    long long start_pos = 0, end_pos = num_snp;
    if((!snp_range_vec.empty()) && (!split_task_vec.empty())){
        spdlog::error("Options --snp-range and --split-task cannot be used together.");
        exit(1);
    }else if(!snp_range_vec.empty()){
        vector<long long> tmp = GenoA.find_bim_index(snp_range_vec);
        start_pos = tmp[0];
        end_pos = tmp[1] + 1;
        this->reset_output_prefix(this->m_out_file, start_pos, end_pos);
    }else if(!split_task_vec.empty()){
        if(split_task_vec[0] < split_task_vec[1] || split_task_vec[1] <= 0){
            spdlog::error("The second number in --split-task must be greater than 0 and not exceed the first number.");
            exit(1);
        }
        long long num_snp_part = num_snp / split_task_vec[0];
        start_pos = num_snp_part * (split_task_vec[1] - 1);
        end_pos = num_snp_part * split_task_vec[1];
        if(split_task_vec[0] == split_task_vec[1]) end_pos = num_snp;
        this->reset_output_prefix(this->m_out_file, split_task_vec[0], split_task_vec[1]);
    }
    if(start_pos >= end_pos){
        spdlog::error("The start SNP position must not exceed the end SNP position.");
        exit(1);
    }
    vector<long long> vec = {start_pos, end_pos};
    return vec;
}
