#define EIGEN_USE_MKL_ALL  // must be before include Eigen



#include <getopt.h>
#include "mkl.h"
#include "omp.h"
#include <set>
#include <map>
#include <cmath>
#include <string>
#include <random>
#include <cstdlib>
#include <gsl/gsl_cdf.h>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>  // Include CLI11


#include "mom.hpp"
#include "../utils/geno.hpp"
#include "../utils/phen.hpp"
#include "../utils/iterator_utils.hpp"
#include "../utils/string_utils.hpp"




using std::set;
using std::map;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::HouseholderQR;
using Eigen::ColPivHouseholderQR;

MoM::MoM(){
    m_maf_cut = 0.01;
    m_missing_rate_cut = 0.05;
}

MoM::~MoM() {
}


double cal_std(Eigen::VectorXd data) {
    double mean = data.mean();
    double variance = (data.array() - mean).square().sum() / (data.size() - 1);
    double stddev = std::sqrt(variance);
    return stddev;
}


void read_data(const string &data_file, const vector<long long> &covariate_arr,
               const vector<long long> &class_arr, const vector<long long> &bye_arr,
               const vector<long long> &trait, const vector<string> &missing_in_data_vec,
               PHEN &phenoA, const vector<string> &in_in_fam_vec,
               vector<string> &id_in_data_vec) {
    vector<long long> col_used_vec = vector_merged({covariate_arr, class_arr, bye_arr, trait});
    phenoA.read(data_file, {0}, {in_in_fam_vec}, col_used_vec, missing_in_data_vec);

    phenoA.get_given_column(0, id_in_data_vec);
    set<string> id_in_data_st = set<string>(id_in_data_vec.begin(), id_in_data_vec.end());

    if (id_in_data_vec.size() > id_in_data_st.size()) {
        spdlog::error("Duplicated sample IDs exist in the data file!");
        exit(1);
    }
}


MatrixXd adjust_for_fixed_effects(PHEN &phenoA, const vector<long long> &covariate_arr,
                                  const vector<long long> &class_arr, const vector<long long> &bye_arr,
                                  const vector<long long> &trait, MatrixXd &y) {
    vector<long long> index_fixed_effect_vec;
    vector<string> fixed_effect_name_vec;
    MatrixXd xmat = phenoA.dmat(covariate_arr, class_arr, index_fixed_effect_vec, fixed_effect_name_vec, trait, y);
    vector<long long> column_to_remove_vec;
    phenoA.dmat_column_full_rank(xmat, index_fixed_effect_vec, fixed_effect_name_vec, column_to_remove_vec);

    if (!column_to_remove_vec.empty()) {
        spdlog::info("Dependent columns of xmat are removed: ");
    
        std::string removed_columns;
        for (auto col : column_to_remove_vec) {
            removed_columns += std::to_string(col) + " ";
        }

        spdlog::info("{}", removed_columns);
    }

    MatrixXd bye_mat;
    if (!bye_arr.empty()) {
        phenoA.get_given_columns_double(bye_arr, bye_mat);
    }

    long long nrows = xmat.rows();
    long long ncols = xmat.cols();
    xmat.conservativeResize(nrows, ncols + bye_mat.cols());
    xmat.rightCols(bye_mat.cols()) = bye_mat;
    ColPivHouseholderQR<MatrixXd> qr(xmat);

    if (qr.rank() < xmat.cols()) {
        spdlog::error("There are dependent columns between interaction environment covariates and other covariates!");
        exit(1);
    }

    y = (y - xmat * ((xmat.transpose() * xmat).inverse() * (xmat.transpose() * y))).eval();

    return bye_mat;
}

void read_plink_bed(const string &bed_file, long long start_snp, long long num_snp_read,
                    long long num_iid_bed, char* &bytes_vec) {
    string in_file = bed_file + ".bed";
    FILE* fin = fopen(in_file.c_str(), "rb");
    if (!fin) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }

    long long num_byte_for_one_snp = (num_iid_bed + 3) / 4;
    std::streamoff startPos = start_snp * num_byte_for_one_snp + 3;
    if (fseek(fin, startPos, SEEK_SET) != 0) {
        spdlog::error("Fail to seek the predefined position");
        exit(1);
    }

    bytes_vec = new char[num_byte_for_one_snp * num_snp_read];
    long long read_count = fread(bytes_vec, sizeof(char), num_byte_for_one_snp * num_snp_read, fin);
    if (read_count != num_byte_for_one_snp * num_snp_read) {
        delete[] bytes_vec;
        spdlog::error("fread failed to read the expected amount");
        exit(1);
    }
    if (ferror(fin)) {
        delete[] bytes_vec;
        spdlog::error("Failed to read from file: " + in_file);
        exit(1);
    }
    fclose(fin);
}

MatrixXd generate_random_Bs(long long num_iid_data, long long num_randomB, const MatrixXd &y) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);
    MatrixXd B = MatrixXd::Zero(num_iid_data, num_randomB + 1);
    B.col(0) = y.col(0);

    for (long long i = 0; i < num_iid_data; i++) {
        for (long long j = 1; j < num_randomB + 1; j++) {
            B(i, j) = d(gen);
        }
    }

    return B;
}

void MoM::process_snps(long long num_snp_read, long long num_iid_bed, long long num_used_id,
                  long long num_randomB, const char* bytes_vec, const vector<long long> &index_vec,
                  MatrixXd &GB, MatrixXd &VB, const MatrixXd &B, const MatrixXd &bye_mat) {
    long long num_byte_for_one_snp = (num_iid_bed + 3) / 4;

    #pragma omp parallel for schedule(dynamic)
    for (long long isnp = 0; isnp < num_snp_read; isnp++) {
        Eigen::VectorXd snp_Vec(num_iid_bed);
        for (long long i = 0; i < num_iid_bed; i++) {
            unsigned char byte = bytes_vec[i / 4 + isnp * num_byte_for_one_snp];
            unsigned char genotype = (byte >> (2 * (i % 4))) & 3;
            double code_val;
            switch (genotype) {
                case 0: code_val = 2.0; break;
                case 1: code_val = -9; break;
                case 2: code_val = 1.0; break;
                case 3: code_val = 0.0; break;
                default: throw std::logic_error("Invalid genotype value");
            }
            snp_Vec[i] = code_val;
        }

        Eigen::VectorXd snp_used_Vec = Eigen::VectorXd::Zero(num_used_id);
        vector<long long> missing_index_vec;
        double sum = 0;
        for (long long i = 0; i < num_used_id; i++) {
            long long code_int = snp_Vec[index_vec[i]];
            if (code_int == -9) {
                missing_index_vec.push_back(i);
            } else {
                snp_used_Vec(i) = code_int;
                sum += code_int;
            }
        }
        long long num_missing = missing_index_vec.size();
        double missing_rate = num_missing * 1.0 / num_used_id;
        if(missing_rate < m_missing_rate_cut){
            double freq = sum / (2 * (num_used_id - num_missing));
            snp_used_Vec = snp_used_Vec.array() - 2 * freq;
            for(auto index:missing_index_vec){
                snp_used_Vec(index) = 0; // replace missing with mean value
            }
            if(freq > m_maf_cut && freq < 1 - m_maf_cut) {
                snp_used_Vec = snp_used_Vec / cal_std(snp_used_Vec);
            }else{
                snp_used_Vec.setZero();
            }
        }else{
            snp_used_Vec.setZero();
        }

        Eigen::MatrixXd local_GB = snp_used_Vec * (snp_used_Vec.transpose() * B) / num_snp_read;
        Eigen::MatrixXd snpByeE = bye_mat.array().colwise() * snp_used_Vec.array();
        Eigen::MatrixXd local_VB = snpByeE * (snpByeE.transpose() * B) / (num_snp_read * bye_mat.cols());

        #pragma omp critical
        {
            GB += local_GB;
            VB += local_VB;
        }
    }
}


void process_snps_MV(long long num_snp_read, long long num_iid_bed, long long num_used_id,
                  long long num_randomB, const char* bytes_vec, const vector<long long> &index_vec,
                  vector<MatrixXd> &VB_vec, const MatrixXd &B, const MatrixXd &bye_mat) {
    long long num_byte_for_one_snp = (num_iid_bed + 3) / 4;

    #pragma omp parallel for schedule(dynamic)
    for (long long isnp = 0; isnp < num_snp_read; isnp++) {
        Eigen::VectorXd snp_Vec(num_iid_bed);
        for (long long i = 0; i < num_iid_bed; i++) {
            unsigned char byte = bytes_vec[i / 4 + isnp * num_byte_for_one_snp];
            unsigned char genotype = (byte >> (2 * (i % 4))) & 3;
            double code_val;
            switch (genotype) {
                case 0: code_val = 2.0; break;
                case 1: code_val = -9; break;
                case 2: code_val = 1.0; break;
                case 3: code_val = 0.0; break;
                default: throw std::logic_error("Invalid genotype value");
            }
            snp_Vec[i] = code_val;
        }

        Eigen::VectorXd snp_used_Vec = Eigen::VectorXd::Zero(num_used_id);
        vector<long long> missing_index_vec;
        double sum = 0;
        for (long long i = 0; i < num_used_id; i++) {
            long long code_int = snp_Vec[index_vec[i]];
            if (code_int == -9) {
                missing_index_vec.push_back(i);
            } else {
                snp_used_Vec(i) = code_int;
                sum += code_int;
            }
        }
        long long num_missing = missing_index_vec.size();
        if(num_missing != num_used_id){
            double freq = sum / (2 * (num_used_id - num_missing));
            snp_used_Vec = snp_used_Vec.array() - 2 * freq;
            for(auto index:missing_index_vec){
                snp_used_Vec(index) = 0; // replace missing with mean value
            }
            if(freq > 0 && freq < 1) {
                snp_used_Vec = snp_used_Vec / cal_std(snp_used_Vec);
            }else{
                snp_used_Vec.setZero();
            }
        }else{
            snp_used_Vec.setZero();
        }

        Eigen::MatrixXd local_GB = snp_used_Vec * (snp_used_Vec.transpose() * B) / num_snp_read;
        Eigen::MatrixXd snpByeE = bye_mat.array().colwise() * snp_used_Vec.array();

        vector <Eigen::MatrixXd> local_VB_vec;
        for(long long i = 0; i < snpByeE.cols(); i++){
            Eigen::MatrixXd local_VB = snpByeE.col(i) * (snpByeE.col(i).transpose() * B) / (num_snp_read * bye_mat.cols());
            local_VB_vec.push_back(local_VB);
        }
        

        #pragma omp critical
        {
            VB_vec[0] += local_GB;
            for(long long i = 0; i < snpByeE.cols(); i++){
                VB_vec[i+1] += local_VB_vec[i];
            }
            
        }
    }
}


VectorXd calculate_variance_components(const MatrixXd &GB, const MatrixXd &VB, const MatrixXd &B, const MatrixXd &y) {
    long long num_randomB = B.cols() - 1;
    long long num_used_id = GB.rows();

    MatrixXd leftMat(3, 3);
    VectorXd rightVec(3);

    leftMat(0, 0) = (GB.rightCols(num_randomB).cwiseProduct(GB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(1, 1) = (VB.rightCols(num_randomB).cwiseProduct(VB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(2, 2) = num_used_id;
    leftMat(0, 1) = leftMat(1, 0) = (GB.rightCols(num_randomB).cwiseProduct(VB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(0, 2) = leftMat(2, 0) = (GB.rightCols(num_randomB).cwiseProduct(B.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(1, 2) = leftMat(2, 1) = (VB.rightCols(num_randomB).cwiseProduct(B.rightCols(num_randomB))).sum() / num_randomB;

    rightVec(0) = (GB.col(0).transpose() * y).sum();
    rightVec(1) = (VB.col(0).transpose() * y).sum();
    rightVec(2) = (y.transpose() * y).sum();

    VectorXd varcom = leftMat.inverse() * rightVec;
    return varcom;
}


VectorXd calculate_variance_components_withNxE(const MatrixXd &GB, const MatrixXd &VB, const VectorXd &nxe_vec, const MatrixXd &B, const MatrixXd &y) {
    long long num_randomB = B.cols() - 1;
    long long num_used_id = GB.rows();

    MatrixXd leftMat(4, 4);
    VectorXd rightVec(4);
    MatrixXd nxeB = B.array().colwise() * nxe_vec.array();

    leftMat(0, 0) = (GB.rightCols(num_randomB).cwiseProduct(GB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(1, 1) = (VB.rightCols(num_randomB).cwiseProduct(VB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(2, 2) = (nxe_vec.cwiseProduct(nxe_vec)).sum();
    leftMat(3, 3) = num_used_id;
    leftMat(0, 1) = leftMat(1, 0) = (GB.rightCols(num_randomB).cwiseProduct(VB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(0, 2) = leftMat(2, 0) = (GB.rightCols(num_randomB).cwiseProduct(nxeB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(0, 3) = leftMat(3, 0) = (GB.rightCols(num_randomB).cwiseProduct(B.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(1, 2) = leftMat(2, 1) = (VB.rightCols(num_randomB).cwiseProduct(nxeB.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(1, 3) = leftMat(3, 1) = (VB.rightCols(num_randomB).cwiseProduct(B.rightCols(num_randomB))).sum() / num_randomB;
    leftMat(2, 3) = leftMat(3, 2) = nxe_vec.sum();

    rightVec(0) = (GB.col(0).transpose() * y).sum();
    rightVec(1) = (VB.col(0).transpose() * y).sum();
    rightVec(2) = (nxeB.col(0).transpose() * y).sum();
    rightVec(3) = (y.transpose() * y).sum();

    VectorXd varcom = leftMat.inverse() * rightVec;
    return varcom;
}

VectorXd calculate_multi_variance_components(const vector<MatrixXd> &VB_vec, const MatrixXd &B, const MatrixXd &y) {
    long long num_randomB = B.cols() - 1;
    long long num_used_id = VB_vec[0].rows();
    long long num_varcom = VB_vec.size() + 1;

    MatrixXd leftMat(num_varcom, num_varcom);
    VectorXd rightVec(num_varcom);

    for(long long i = 0; i < num_varcom - 1; i++){
        leftMat(i, i) = (VB_vec[i].rightCols(num_randomB).cwiseProduct(VB_vec[i].rightCols(num_randomB))).sum() / num_randomB;
        for(long long j = i + 1; j < num_varcom - 1; j ++){
            leftMat(i, j) = leftMat(j, i) = (VB_vec[i].rightCols(num_randomB).cwiseProduct(VB_vec[j].rightCols(num_randomB))).sum() / num_randomB;
        }

        leftMat(i, num_varcom - 1) = leftMat(num_varcom - 1, i) = (VB_vec[i].rightCols(num_randomB).cwiseProduct(B.rightCols(num_randomB))).sum() / num_randomB;
    }
    leftMat(num_varcom - 1, num_varcom - 1) = num_used_id;

    for(long long i = 0; i < num_varcom - 1; i++){
        rightVec(i) = (VB_vec[i].col(0).transpose() * y).sum();
    }
    rightVec(num_varcom - 1) = (y.transpose() * y).sum();

    VectorXd varcom = leftMat.inverse() * rightVec;
    return varcom;
}


void MoM::MoMV(bool no_noisebye, const string &out_file, const string &data_file, const vector<long long> &covariate_arr,
          const vector<long long> &class_arr, const vector<long long> &bye_arr, const vector<long long> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, long long start_snp, long long num_snp_read,
          long long num_randomB) {
    spdlog::info("Start analysis...");

    long long num_trait = trait.size();
    long long num_covariate = covariate_arr.size();
    long long num_class = class_arr.size();
    long long num_bye = bye_arr.size();

    spdlog::info("Read IDs from Plink .fam file");
    GENO genoA(bed_file);
    long long num_iid_bed = genoA.get_num_iid();
    long long num_sid_bed = genoA.get_num_sid();
    vector<string> in_in_fam_vec = genoA.iid_vec();

    spdlog::info("Read data file");
    PHEN phenoA;
    vector<string> id_in_data_vec;
    read_data(data_file, covariate_arr, class_arr, bye_arr, trait, missing_in_data_vec, phenoA, in_in_fam_vec, id_in_data_vec);

    long long num_iid_data = id_in_data_vec.size();
    spdlog::info("The number of used iids in data file: {}", num_iid_data);

    spdlog::info("Design matrix for fixed effects");
    MatrixXd y;
    MatrixXd bye_mat = adjust_for_fixed_effects(phenoA, covariate_arr, class_arr, bye_arr, trait, y);

    spdlog::info("Read the plink bed file");
    char* bytes_vec = nullptr;
    read_plink_bed(bed_file, start_snp, num_snp_read, num_iid_bed, bytes_vec);

    spdlog::info("Random Bs");
    MatrixXd B = generate_random_Bs(num_iid_data, num_randomB, y);

    vector<long long> index_vec = genoA.find_fam_index_pro(id_in_data_vec);
    long long num_used_id = index_vec.size();

    MatrixXd GB = MatrixXd::Zero(num_used_id, num_randomB + 1);
    MatrixXd VB = MatrixXd::Zero(num_used_id, num_randomB + 1);
    VectorXd nxe_vec = bye_mat.rowwise().squaredNorm() / bye_mat.cols(); // nxe for each individual

    process_snps(num_snp_read, num_iid_bed, num_used_id, num_randomB, bytes_vec, index_vec, GB, VB, B, bye_mat);

    VectorXd varcom;
    if(no_noisebye){
        varcom = calculate_variance_components(GB, VB, B, y);
    }else{
        varcom = calculate_variance_components_withNxE(GB, VB, nxe_vec, B, y);
    }

    delete[] bytes_vec;

    for (auto val : varcom) {
        spdlog::info(std::to_string(val), 0, 0, " ");
    }
    spdlog::info("");

    ofstream fout(out_file + ".var");
    if(!fout.is_open()){
        spdlog::error("Fail to open the output file");
        exit(1);
    }
    for (auto val : varcom) {
        fout << val << std::endl;
    }
    fout.close();
}



void MoM::MoMMV(const string &out_file, const string &data_file, const vector<long long> &covariate_arr,
          const vector<long long> &class_arr, const vector<long long> &bye_arr, const vector<long long> &trait,
          const vector<string> &missing_in_data_vec, const string &bed_file, long long start_snp, long long num_snp_read,
          long long num_randomB) {
    spdlog::info("Start analysis...");

    long long num_trait = trait.size();
    long long num_covariate = covariate_arr.size();
    long long num_class = class_arr.size();
    long long num_bye = bye_arr.size();

    spdlog::info("Read IDs from Plink .fam file");
    GENO genoA(bed_file);
    long long num_iid_bed = genoA.get_num_iid();
    long long num_sid_bed = genoA.get_num_sid();
    vector<string> in_in_fam_vec = genoA.iid_vec();

    spdlog::info("Read data file");
    PHEN phenoA;
    vector<string> id_in_data_vec;
    read_data(data_file, covariate_arr, class_arr, bye_arr, trait, missing_in_data_vec, phenoA, in_in_fam_vec, id_in_data_vec);

    long long num_iid_data = id_in_data_vec.size();
    spdlog::info("The number of used iids in data file: {}", num_iid_data);


    spdlog::info("Design matrix for fixed effects");
    MatrixXd y;
    MatrixXd bye_mat = adjust_for_fixed_effects(phenoA, covariate_arr, class_arr, bye_arr, trait, y);

    spdlog::info("Read the plink bed file");
    char* bytes_vec = nullptr;
    read_plink_bed(bed_file, start_snp, num_snp_read, num_iid_bed, bytes_vec);

    spdlog::info("Random Bs");
    MatrixXd B = generate_random_Bs(num_iid_data, num_randomB, y);

    vector<long long> index_vec = genoA.find_fam_index_pro(id_in_data_vec);
    long long num_used_id = index_vec.size();

    vector<MatrixXd> VB_vec;
    for(long long i = 0; i < num_bye+1; i++){
        VB_vec.push_back(MatrixXd::Zero(num_used_id, num_randomB + 1));
    }

    process_snps_MV(num_snp_read, num_iid_bed, num_used_id, num_randomB, bytes_vec, index_vec,VB_vec, B, bye_mat);

    VectorXd varcom = calculate_multi_variance_components(VB_vec, B, y);

    delete[] bytes_vec;

    for (auto val : varcom) {
        spdlog::info(std::to_string(val));
    }
    spdlog::info("");
}



int MoM::run(int argc, char* argv[]) {
    CLI::App app{"MoM - GxE Heritability Estimation Using the Method of Moments"};

    app.description(R"(
    Quick Start:
        fastgxe --mom --bfile bed_file --data data_file --trait BMI --env-int age smoking:alcohol --out test_mom
    )");


    int threads = 10;
    std::string out_file, data_file, bed_file;
    std::vector<long long> block_vec;
    long long num_randomB = 10;
    std::vector<std::string> covariate_vec, class_vec, bye_vec, trait;
    std::vector<std::string> missing_in_data_vec = {"NA", "Na", "na", "NAN", "NaN", "nan", "-NAN", "-NaN", "-nan", "<NA>", "<na>", "N/A", "n/a"};

    // Define command-line options
    bool mom=false;
    bool no_noisebye = false;
    app.add_flag("--mom", mom, "GxE Heritability Estimation Using MoM");

    app.add_flag("--no-noisebye", [&no_noisebye](int count) {
        if (count > 0) no_noisebye = true;  // Flip no_noisebye to true when flag is used
    }, "Disable noise-by-environment interaction terms.\n"
       "  - Noise-by-environment interactions are ENABLED by default.\n"
       "  - Use --no-noisebye to turn it OFF.");

    app.add_option("--threads", threads, "Number of threads (default: 10)")->default_val(10);
    app.add_option("--out", out_file, "Output file")->required();
    app.add_option("--data", data_file, "Input data file")->required();
    app.add_option("--bfile", bed_file, "Binary PLINK file")->required();

    app.add_option("--env-int", bye_vec, 
        "List of interacting environmental covariates (required).\n"
        "  - Supports multiple values (e.g., --env-int age BMI smoking).\n"
        "  - Use ':' to specify a range (e.g., --env-int age:BMI).\n"
        "  - Expands to include all covariates in the specified range.\n"
        "  - Example order: {age, BMI, smoking, exercise, alcohol}.\n"
        "  - Invalid or duplicate values will cause an error.")
        ->expected(-1)->required();
    app.add_option("--trait", trait, "Trait columns")->required();

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
    
    app.add_option("--maf", m_maf_cut, 
        "Minor allele frequency (MAF) cutoff (default: 0.01).")
        ->default_val(0.01);

    app.add_option("--missing-rate", m_missing_rate_cut, 
        "Missing genotype rate cutoff (default: 0.05).")
        ->default_val(0.05);
    
    app.add_option("--block", block_vec, "Block size (2 values required)")->expected(2);
    app.add_option("--num-random-vector", num_randomB, "Number of random vector")->default_val(10);
    app.add_option("--missing-data", missing_in_data_vec, 
        "List of missing value indicators for phenotype/covariates.\n"
        "  - Default: {NA, Na, na, NAN, NaN, nan, -NAN, -NaN, -nan, <NA>, <na>, N/A, n/a}.\n"
        "  - Customize with space-separated values (e.g., --missing-data . -999 \"?\").")
        ->expected(-1);
    
    // Parse command-line arguments
    CLI11_PARSE(app, argc, argv);

    // Set number of threads for MKL and OpenMP
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

    // Trait index
    
    std::vector<std::string> strNoFound_vec;
    std::vector<long long> trait_index_vec = find_index(head_vec, trait, strNoFound_vec);

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

    // Load genotype data
    GENO genoA(bed_file);
    long long num_sid_bed = genoA.get_num_sid();
    long long start_pos = 0, num_snp_read = num_sid_bed;

    if (!block_vec.empty()) {
        if (block_vec[0] < block_vec[1] || block_vec[1] <= 0) {
            spdlog::error("Invalid block parameters: second number must be >0 and <= first number");
            exit(1);
        }
        long long num_snp_part = num_sid_bed / block_vec[0];
        start_pos = num_snp_part * (block_vec[1] - 1);
        long long end_pos = num_snp_part * block_vec[1];
        if (block_vec[0] == block_vec[1]) end_pos = num_sid_bed;
        num_snp_read = end_pos - start_pos;
    }

    // Run the appropriate MoM analysis method

    this->MoMV(no_noisebye, out_file, data_file, covariate_index_vec, class_index_vec, bye_index_vec, trait_index_vec, missing_in_data_vec, bed_file, start_pos, num_snp_read, num_randomB);

    return 0;
}
