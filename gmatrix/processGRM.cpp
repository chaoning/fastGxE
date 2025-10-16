/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-08-12 17:23:31
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-31 12:43:53
 */

#define EIGEN_USE_MKL_ALL  // !must be before include Eigen

#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <string>
#include <spdlog/spdlog.h>
#include <CLI/CLI.hpp>
#include "mkl.h"
#include "omp.h"
#include <Eigen/Eigen>

#include "processGRM.hpp"
#include "string_utils.hpp"
#include "iterator_utils.hpp"
#include "EigenMatrix_utils.hpp"


using std::ofstream; using std::ifstream;
using std::cout; using std::endl;
using std::getline;

using std::string;
using std::map;
using std::vector;
using std::set;

// typedef Eigen::SparseMatrix<double> SpMat; 
// typedef Eigen::Triplet<double> T;

ProcessGRM::ProcessGRM(){
}

ProcessGRM::~ProcessGRM(){
}


vector<string> ProcessGRM::read_grm_id(string grm_file){
    ifstream fin;
    fin.open(grm_file + ".id");
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}.id", grm_file);
        exit(EXIT_FAILURE);
    }

    vector<string> vec;
    string line;
    while(getline(fin, line)){
        process_line(line);
        vec.push_back(line);
    }
    fin.close();

    if(is_duplicated(vec)){
        spdlog::error("Duplicated iids exist in the grm file");
        exit(EXIT_FAILURE);
    }
    return vec;
}


vector<string> ProcessGRM::read_gcta_grm_id(string grm_file){
    ifstream fin;
    fin.open(grm_file + ".id");
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}.id", grm_file);
        exit(EXIT_FAILURE);
    }

    vector<string> vec;
    string line;
    while(getline(fin, line)){
        process_line(line);
        vector<string> tmp_vec = split_string(line);
        vec.push_back(tmp_vec[1]);
    }
    fin.close();

    if(is_duplicated(vec)){
        spdlog::error("Duplicated iids exist in the grm file");
        exit(EXIT_FAILURE);
    }
    return vec;
}



void ProcessGRM::out_grm_id(std::string out_file, std::vector<std::string> grm_id_vec){
    ofstream fout;
    fout.open(out_file + ".id");
    if(!fout.is_open()){
        spdlog::error("Fail to open: {}", out_file);
        exit(EXIT_FAILURE);
    }

    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        fout << *it << endl;
    }
    fout.close();
}


// Read the grm in the given order
void ProcessGRM::read_grm_bin(string grm_file, map<string, long long> grm_id_map, Eigen::MatrixXd& mat, double val){
    
    vector<string> grm_id_vec = ProcessGRM::read_grm_id(grm_file);
    vector<long long> grm_id_index_vec;
    long long num_id_used = 0;
    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        auto it2 = grm_id_map.find(*it);
        if(it2 != grm_id_map.end()){
            grm_id_index_vec.push_back(it2->second);
            num_id_used++;
        }else{
            grm_id_index_vec.push_back(-1);
        }
    }

    if(grm_id_map.empty()){
        // read in the original order
        grm_id_index_vec.clear(); // must clear
        num_id_used = grm_id_vec.size();
        for(long long i = 0; i < num_id_used; i++){
            grm_id_map[grm_id_vec[i]] = i;
            grm_id_index_vec.push_back(i);
        }
    }


    if(grm_id_map.size() != num_id_used){
        spdlog::error("iids in the grm do not match the data");
        exit(1);
    }

    ifstream fin(grm_file + ".bin", std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: " + grm_file + ".bin");
        exit(1);
    }

    mat.resize(num_id_used, num_id_used);
    long long num_id = grm_id_index_vec.size();
    double tmp_val;
    long long index0, index1;
    for(long long i = 0; i < num_id; i++){
        index0 = grm_id_index_vec[i];
        for(long long j = 0; j <= i; j++){
            index1 = grm_id_index_vec[j];
            fin.read((char*)&tmp_val, sizeof(double));
            if(index0 != -1 && index1 != -1){
                if(index0 == index1) tmp_val += val;
                mat(index0, index1) = mat(index1, index0) = tmp_val;
            }
        }
    }

    fin.close();
}


void ProcessGRM::read_gcta_grm_bin(string grm_file, map<string, long long> grm_id_map, Eigen::MatrixXd& mat, double val){
    
    vector<string> grm_id_vec = ProcessGRM::read_gcta_grm_id(grm_file);
    vector<long long> grm_id_index_vec;
    long long num_id_used = 0;
    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        auto it2 = grm_id_map.find(*it);
        if(it2 != grm_id_map.end()){
            grm_id_index_vec.push_back(it2->second);
            num_id_used++;
        }else{
            grm_id_index_vec.push_back(-1);
        }
    }

    if(grm_id_map.empty()){
        // read in the original order
        grm_id_index_vec.clear(); // must clear
        num_id_used = grm_id_vec.size();
        for(long long i = 0; i < num_id_used; i++){
            grm_id_map[grm_id_vec[i]] = i;
            grm_id_index_vec.push_back(i);
        }
    }


    if(grm_id_map.size() != num_id_used){
        spdlog::error("iids in the grm do not match the data");
        exit(1);
    }

    ifstream fin(grm_file + ".bin", std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}.bin", grm_file);
        exit(1);
    }

    mat.resize(num_id_used, num_id_used);
    long long num_id = grm_id_index_vec.size();
    float tmp_val;
    long long index0, index1;
    for(long long i = 0; i < num_id; i++){
        index0 = grm_id_index_vec[i];
        for(long long j = 0; j <= i; j++){
            index1 = grm_id_index_vec[j];
            fin.read((char*)&tmp_val, sizeof(float));
            if(index0 != -1 && index1 != -1){
                if(index0 == index1) tmp_val += val;
                mat(index0, index1) = mat(index1, index0) = tmp_val;
            }
        }
    }

    fin.close();
}



void ProcessGRM::read_grm_sp_bin(std::string grm_sparse_file, std::map<std::string, long long> grm_id_map, Eigen::SparseMatrix<double>& mat, double val){
    vector<string> grm_id_vec = ProcessGRM::read_grm_id(grm_sparse_file);
    vector<long long> grm_id_index_vec;
    long long num_id_used = 0;
    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        auto it2 = grm_id_map.find(*it);
        if(it2 != grm_id_map.end()){
            grm_id_index_vec.push_back(it2->second);
            num_id_used++;
        }else{
            grm_id_index_vec.push_back(-1);
        }
    }

    if(grm_id_map.empty()){
        // read in the original order
        grm_id_index_vec.clear(); // must clear
        num_id_used = grm_id_vec.size();
        for(long long i = 0; i < num_id_used; i++){
            grm_id_map[grm_id_vec[i]] = i;
            grm_id_index_vec.push_back(i);
        }
    }
    
    if(grm_id_map.size() != num_id_used){
        spdlog::error("iids in the grm do not match the data");
        exit(1);
    }

    ifstream fin(grm_sparse_file + ".sp.bin", std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}.sp.bin", grm_sparse_file);
        exit(1);
    }


    std::vector<Eigen::Triplet<double>> tripletList;
    double tmp_val;
    long long index0, index1, index_curr0, index_curr1;

    while(fin.read((char*)&index0, sizeof(long long)) && fin.read((char*)&index1, sizeof(long long)) && fin.read((char*)&tmp_val, sizeof(double))){
        // cout << index0 + 1 << " " << index1 + 1 << " " << tmp_val << endl;
        index_curr0 = grm_id_index_vec[index0];
        index_curr1 = grm_id_index_vec[index1];
        if(index_curr0 != -1 && index_curr1 != -1){
            if(index_curr0 == index_curr1){
                tmp_val += val;
                tripletList.push_back(Eigen::Triplet<double>(index_curr0, index_curr1, tmp_val));
            }else{
                tripletList.push_back(Eigen::Triplet<double>(index_curr0, index_curr1, tmp_val));
                tripletList.push_back(Eigen::Triplet<double>(index_curr1, index_curr0, tmp_val));
            }
        }
    }
    fin.close();
    
    mat.resize(num_id_used, num_id_used);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
}


/**
 * @brief merge the grm files from n parts produced by --make-grm --npart n i
 * 
 * @param grm_file 
 * @param sparse 
 * @param npart 
 */
void ProcessGRM::merge_grm(std::string grm_file, int sparse, int npart, string out_file){

    // id
    spdlog::info(" *.id file");
    
    vector<string> file_vec;
    for(int i = 1; i <= npart; i++){
        string tmp_file = grm_file +"." + std::to_string(npart) + "_" + std::to_string(i) + ".id";
        file_vec.push_back(tmp_file);
    }
    ProcessGRM::merge_file(file_vec, out_file + ".id", false);

    // N
    spdlog::info(" *.N file");
    file_vec.clear();
    for(int i = 1; i <= npart; i++){
        string tmp_file = grm_file +"." + std::to_string(npart) + "_" + std::to_string(i) + ".N";
        file_vec.push_back(tmp_file);
    }
    ProcessGRM::merge_file({file_vec[0]}, out_file + ".N", false);


    // bin
    spdlog::info(" *.bin file");
    file_vec.clear();
    for(int i = 1; i <= npart; i++){
        if(sparse == 0){
            string tmp_file = grm_file +"." + std::to_string(npart) + "_" + std::to_string(i) + ".bin";
            file_vec.push_back(tmp_file);
        }else{
            string tmp_file = grm_file +"." + std::to_string(npart) + "_" + std::to_string(i) + ".sp.bin";
            file_vec.push_back(tmp_file);
        }
    }

    if(sparse == 0){
        ProcessGRM::merge_file(file_vec, out_file + ".bin", true);
    }else{
        ProcessGRM::merge_file(file_vec, out_file + ".sp.bin", true);
    }
}

/**
 * @brief merge files
 * 
 * @param file_vec 
 * @param out_file 
 * @param bin 
 */
void ProcessGRM::merge_file(vector<string> file_vec, string out_file, bool bin=false){

    if(bin == true){

        ofstream fout(out_file, std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the: {}", out_file);
            exit(1);
        }

        for(auto it = file_vec.begin(); it != file_vec.end(); it++){
            ifstream fin(*it, std::ios::binary);
            if(!fin.is_open()){
                spdlog::error("Fail to open the: {}", *it);
                exit(1);
            }
            spdlog::info("Merging file: {}", *it);

            char tmp;
            while(fin.read(&tmp, sizeof(char))){
                fout.write(&tmp, sizeof(char));
            }
            fin.close();
        }
        fout.close();

    }else{
        ofstream fout(out_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the: {}", out_file);
            exit(1);
        }

        for(auto it = file_vec.begin(); it != file_vec.end(); it++){
            ifstream fin(*it);
            if(!fin.is_open()){
                spdlog::error("Fail to open the: {}", *it);
                exit(1);
            }
            spdlog::info("Merging file: {}", *it);

            string line;
            while(std::getline(fin, line)){
                fout << line << endl;
            }
            fin.close();
        }
        fout.close();
    }
    
}



void ProcessGRM::out_gmat(Eigen::MatrixXd& gmat, vector<string> grm_id_vec, int out_fmt, string out_file){
    
    long long num_id = grm_id_vec.size();
    string out_gmat_file = out_file;
    if(out_fmt == 0){ // matrix
        out_gmat_file += ".mat_fmt";
        ofstream fout;
        fout.open(out_gmat_file, std::ios::out);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        for(long long i = 0; i < num_id; i++){
            for(long long j = 0; j < num_id; j++){
                fout << gmat(i, j) << " ";
            }
            fout << endl;
        }
        fout.close();
    }else if(out_fmt == 1){ //row col val
        out_gmat_file += ".ind_fmt";
        ofstream fout;
        fout.open(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        for(long long i = 0; i < num_id; i++){
            for(long long j = 0; j <= i; j++){
                fout << i+1 << " " << j+1 << " " << gmat(i, j) << endl;
            }
        }
        fout.close();
    }else if(out_fmt == 2){ //id1 id2 val
        out_gmat_file += ".id_fmt";
        ofstream fout;
        fout.open(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        for(long long i = 0; i < num_id; i++){
            for(long long j = 0; j <= i; j++){
                fout << grm_id_vec[i] << " " << grm_id_vec[j] << " " << gmat(i, j) << endl;
            }
        }
        fout.close();
    }else if(out_fmt == 3){ // .bin
        std::ofstream fout(out_file + ".bin", std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}.bin", out_file);
            exit(1);
        }

        for(long long i = 0; i < num_id; i++){
            for(long long j = 0; j <= i; j++){
                fout.write((char*)&gmat(i, j), sizeof(double));
            }
        }
        fout.close();
    }else{
        spdlog::error("Unsupported output format: --out-fmt {}", out_fmt);
        exit(1);
    }
}


void ProcessGRM::out_gmat(Eigen::SparseMatrix<double>& gmat, vector<string> grm_id_vec, int out_fmt, std::string out_file){
    long long num_id = grm_id_vec.size();
    string out_gmat_file = out_file;
    if(out_fmt == 0){ // matrix
        spdlog::warn("Output sparse relationship matrix in dense format");
        out_gmat_file += ".mat_fmt";
        ofstream fout;
        fout.open(out_gmat_file, std::ios::out);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        fout << Eigen::MatrixXd(gmat) << endl;
        fout.close();
    }else if(out_fmt == 1){ //row col val
        out_gmat_file += ".ind_fmt";
        ofstream fout;
        fout.open(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        gmat = gmat.triangularView<Eigen::Upper>();  // keep the upper triangular part
        for (long long k = 0; k < gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(gmat, k); it; ++it){
                fout << it.col() + 1 << " " << it.row() + 1 << " " << it.value() << endl; // output the transpose of upper triangular part
            }
        }
        fout.close();
    }else if(out_fmt == 2){ //id1 id2 val
        out_gmat_file += ".id_fmt";
        ofstream fout;
        fout.open(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        gmat = gmat.triangularView<Eigen::Upper>();  // keep the upper triangular part
        for (long long k = 0; k < gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(gmat, k); it; ++it){
                fout << grm_id_vec[it.col()] << " " << grm_id_vec[it.row()] << " " << it.value() << endl; // output the transpose of upper triangular part
            }
        }

        fout.close();
    }else if(out_fmt == 3){
        std::ofstream fout(out_file + ".sp.bin", std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}.sp.bin", out_file);
            exit(1);
        }
        long long index0, index1;
        double val;
        gmat = gmat.triangularView<Eigen::Upper>();  // keep the upper triangular part
        for (long long k = 0; k < gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(gmat, k); it; ++it){
                //output the transpose of upper triangular part -> low triangular part
                index0 = it.col(), index1 = it.row();
                val = it.value();
                fout.write((char*)&index0, sizeof(long long));
                fout.write((char*)&index1, sizeof(long long));
                fout.write((char*)&val, sizeof(double));
            }
        }
    }else{
        spdlog::error("Unsupported output format: --out-fmt " + std::to_string(out_fmt));
        exit(1);
    }
}



void ProcessGRM::pca(Eigen::MatrixXd& gmat, std::string out_file){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(gmat);
    if (eigensolver.info() != Eigen::Success) {
        spdlog::error("Fail to compute eigenvalues and eigenvectors");
        exit(1);
    }
    
    Eigen::VectorXd gmat_eigenvals = eigensolver.eigenvalues();
    Eigen::MatrixXd gmat_eigenvecs = eigensolver.eigenvectors();
    long long num_element = gmat_eigenvals.size();
    gmat_eigenvals = (gmat_eigenvals(Eigen::lastN(num_element).reverse())).eval();
    gmat_eigenvecs = (gmat_eigenvecs(Eigen::all, Eigen::lastN(num_element).reverse())).eval();
    ofstream fout1, fout2;
    fout1.open(out_file + ".eigenvals");
    fout2.open(out_file + ".eigenvecs");
    if((!fout1.is_open()) || (!fout2.is_open())){
        spdlog::error("Fail to open the output file: {}.eigenvals, {}.eigenvecs", out_file , out_file);
        exit(1);
    }
    for(long long i = 0; i < gmat_eigenvals.size(); i++){
        fout1 << gmat_eigenvals(i) << endl;
    }
    fout1.close();
    for(long long i = 0; i < gmat_eigenvecs.rows(); i++){
        for(long long j = 0; j < gmat_eigenvecs.cols() - 1; j++){
            fout2 << gmat_eigenvecs(i, j) << " ";
        }
        fout2 << gmat_eigenvecs(i, gmat_eigenvecs.cols() - 1) << endl;
    }
    fout2.close();
}


void ProcessGRM::epistatic_grm(std::string grm_file, int sparse, std::string grm_type, std::string out_file){

    string in_file = grm_file;

    if(grm_type == "aagrm" || grm_type == "ddgrm_as" || grm_type == "ddgrm_gs"){

        if(grm_type == "aagrm"){
            in_file += ".agrm";
        }else if(grm_type == "ddgrm_as"){
            in_file += ".dgrm_as";
        }else{
            in_file += ".dgrm_gs";
        }
        
        ProcessGRM::merge_file({in_file + ".id"}, out_file + "." + grm_type + ".id", false);
        vector<string> grm_id_vec = ProcessGRM::read_grm_id(out_file + "." + grm_type);

        ProcessGRM::merge_file({in_file + ".N", in_file + ".N"}, out_file + "." + grm_type + ".N", false);

        if(sparse == 0){
            std::map<string, long long> grm_id_map;
            Eigen::MatrixXd mat;
            ProcessGRM::read_grm_bin(in_file, grm_id_map, mat);
            mat = mat.cwiseProduct(mat);
            ProcessGRM::out_gmat(mat, grm_id_vec, 3, out_file + "." + grm_type);
        }else{
            std::map<string, long long> grm_id_map;
            Eigen::SparseMatrix<double> mat;
            ProcessGRM::read_grm_sp_bin(in_file, grm_id_map, mat);
            mat = mat.cwiseProduct(mat);
            ProcessGRM::out_gmat(mat, grm_id_vec, 3, out_file + "." + grm_type);
        }
    }else if(grm_type == "adgrm_as" || grm_type == "adgrm_gs"){
        string in_file1;
        string in_file2;
        if(grm_type == "adgrm_as"){
            in_file1 = in_file + ".agrm";
            in_file2 = in_file + ".dgrm_as";
        }else if(grm_type == "adgrm_gs"){
            in_file1 = in_file + ".agrm";
            in_file2 = in_file + ".dgrm_gs";
        }
        
        ProcessGRM::merge_file({in_file1 + ".id"}, out_file + "." + grm_type + ".id", false);
        vector<string> grm_id_vec = ProcessGRM::read_grm_id(out_file + "." + grm_type);

        // N
        ProcessGRM::merge_file({in_file1 + ".N", in_file2 + ".N"}, out_file + "." + grm_type + ".N", false);

        // bin
        if(sparse == 0){
            std::map<string, long long> grm_id_map;
            Eigen::MatrixXd mat1;
            ProcessGRM::read_grm_bin(in_file1, grm_id_map, mat1);
            Eigen::MatrixXd mat2;
            ProcessGRM::read_grm_bin(in_file2, grm_id_map, mat2);
            mat1 = mat1.cwiseProduct(mat2);
            ProcessGRM::out_gmat(mat1, grm_id_vec, 3, out_file + "." + grm_type);
        }else{
            std::map<string, long long> grm_id_map;
            Eigen::SparseMatrix<double> mat1;
            ProcessGRM::read_grm_sp_bin(in_file1, grm_id_map, mat1);
            Eigen::SparseMatrix<double> mat2;
            ProcessGRM::read_grm_sp_bin(in_file2, grm_id_map, mat2);
            mat1 = mat1.cwiseProduct(mat2);
            ProcessGRM::out_gmat(mat1, grm_id_vec, 3, out_file + "." + grm_type);
        }
    }else{
        spdlog::error("Unsupported epistatic GRM type: --make-epis " + grm_type);
        exit(1);
    }
}


void ProcessGRM::remove_id(std::string grm_file, int sparse, std::string id_file, string out_file){
    ifstream fin(id_file);
    if(!fin.is_open()){
        spdlog::error("Fail to open: " + id_file);
        exit(1);
    }

    string line;
    vector<string> id_remove_vec;
    while (getline(fin, line))
    {
        process_line(line);
        id_remove_vec.push_back(line);
    }
    fin.close();
    
    std::vector<std::string> grm_id_vec = ProcessGRM::read_grm_id(grm_file);
    std::map<std::string, long long> grm_id_map;
    vector<string> id_keep_vec;
    long long k = 0;
    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        if(std::find(id_remove_vec.begin(), id_remove_vec.end(), *it) == id_remove_vec.end()){
            grm_id_map[*it] = k;
            id_keep_vec.push_back(*it);
            k++;
        }
    }

    ProcessGRM::out_grm_id(out_file, id_keep_vec);
    ProcessGRM::merge_file({grm_file + ".N"}, out_file + ".N", false);

    if(sparse==0){
        Eigen::MatrixXd mat;
        ProcessGRM::read_grm_bin(grm_file, grm_id_map, mat);
        ProcessGRM::out_gmat(mat, id_keep_vec, 3, out_file);
    }else{
        Eigen::SparseMatrix<double> mat;
        ProcessGRM::read_grm_sp_bin(grm_file, grm_id_map, mat);
        ProcessGRM::out_gmat(mat, id_keep_vec, 3, out_file);
    }
}


void ProcessGRM::keep_id(std::string grm_file, int sparse, std::string id_file, string out_file){
    ifstream fin(id_file);
    if(!fin.is_open()){
        spdlog::error("Fail to open: " + id_file);
        exit(1);
    }

    string line;
    vector<string> id_keep_all_vec;
    while (getline(fin, line))
    {
        process_line(line);
        id_keep_all_vec.push_back(line);
    }
    fin.close();
    
    std::vector<std::string> grm_id_vec = ProcessGRM::read_grm_id(grm_file);
    std::map<std::string, long long> grm_id_map;
    vector<string> id_keep_vec;
    long long k = 0;
    for(auto it = grm_id_vec.begin(); it != grm_id_vec.end(); it++){
        if(std::find(id_keep_all_vec.begin(), id_keep_all_vec.end(), *it) != id_keep_all_vec.end()){
            grm_id_map[*it] = k;
            id_keep_vec.push_back(*it);
            k++;
        }
    }

    ProcessGRM::out_grm_id(out_file, id_keep_vec);
    ProcessGRM::merge_file({grm_file + ".N"}, out_file + ".N", false);

    if(sparse==0){
        Eigen::MatrixXd mat;
        ProcessGRM::read_grm_bin(grm_file, grm_id_map, mat);
        ProcessGRM::out_gmat(mat, id_keep_vec, 3, out_file);
    }else{
        Eigen::SparseMatrix<double> mat;
        ProcessGRM::read_grm_sp_bin(grm_file, grm_id_map, mat);
        ProcessGRM::out_gmat(mat, id_keep_vec, 3, out_file);
    }
}


void ProcessGRM::remove_close_id(std::string grm_file, int sparse, double cutoff, std::string out_file){
    vector<string> id_in_grm_vec = ProcessGRM::read_grm_id(grm_file);
    map<string, set<string>> close_id_map;
    long long num_close = 0;
    for(auto it = id_in_grm_vec.begin(); it != id_in_grm_vec.end(); it++){
        close_id_map[*it] = {};
    }
    if(sparse==0){
        ifstream fin(grm_file + ".bin", std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".bin");
            exit(1);
        }

        long long num_id = id_in_grm_vec.size();
        double tmp_val;
        for(long long i = 0; i < num_id; i++){
            for(long long j = 0; j <= i; j++){
                fin.read((char*)&tmp_val, sizeof(double));
                if(tmp_val > cutoff && i != j){
                    num_close++;
                    close_id_map[id_in_grm_vec[i]].insert(id_in_grm_vec[j]);
                    close_id_map[id_in_grm_vec[j]].insert(id_in_grm_vec[i]);
                }
            }
        }
        fin.close();
    }else{
        ifstream fin(grm_file + ".sp.bin", std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".sp.bin");
            exit(1);
        }

        double tmp_val;
        long long index0, index1;
        while(fin.read((char*)&index0, sizeof(long long)) && fin.read((char*)&index1, sizeof(long long)) && fin.read((char*)&tmp_val, sizeof(double))){
            if(tmp_val > cutoff && index0 != index1){
                num_close++;
                close_id_map[id_in_grm_vec[index0]].insert(id_in_grm_vec[index1]);
                close_id_map[id_in_grm_vec[index1]].insert(id_in_grm_vec[index0]);
            }
        }
        fin.close();
    }

    spdlog::info("There are " + std::to_string(num_close) + " pair of close individuals");

    vector <string> key_vec;
    for(auto it = close_id_map.begin(); it != close_id_map.end(); it++){
        key_vec.push_back(it->first);
    }

    for(auto it = key_vec.begin(); it != key_vec.end(); it++){
        if(close_id_map[*it].size() == 0){
            close_id_map.erase(*it);
        }
    }

    vector <string> remove_id_vec;
    while(close_id_map.size() != 0){

        long long max_num_close = 0;
        string closest_id;

        for(auto it = close_id_map.begin(); it != close_id_map.end(); it++){
            if(it->second.size() > max_num_close){
                max_num_close = it->second.size();
                closest_id = it->first;
            }
        }

        remove_id_vec.push_back(closest_id);

        set<string> removed_set = close_id_map[closest_id];
        close_id_map.erase(closest_id);

        for(auto it = removed_set.begin(); it != removed_set.end(); it++){
            close_id_map[*it].erase(closest_id);
        }

        vector <string> key_vec;
        for(auto it = close_id_map.begin(); it != close_id_map.end(); it++){
            key_vec.push_back(it->first);
        }

        for(auto it = key_vec.begin(); it != key_vec.end(); it++){
            if(close_id_map[*it].size() == 0){
                close_id_map.erase(*it);
            }
        }
    }

    ofstream fout(out_file + ".removed_id");
    if(!fout.is_open()){
        spdlog::error("Fail to open: " + out_file + ".removed_id");
        exit(1);
    }

    for(auto it = remove_id_vec.begin(); it != remove_id_vec.end(); it++){
        fout << *it << endl;
    }
    fout.close();

}


void ProcessGRM::group(std::string grm_file, int sparse, double cutoff, std::string out_file){
    vector<string> id_in_grm_vec = ProcessGRM::read_grm_id(grm_file);
    long long num_id_in_grm = id_in_grm_vec.size();
    set<long long> id_index_in_grm_set;
    map<long long, set<long long>> group_map;

    spdlog::info("Read");
    for(long long i = 0; i < id_in_grm_vec.size(); i++){
        set<long long> vec;
        group_map[i] = vec;
        id_index_in_grm_set.insert(i);
    }

    if(sparse == 0){
        ifstream fin(grm_file + ".bin", std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".bin");
            exit(1);
        }

        double tmp_val;
        for(long long i = 0; i < id_in_grm_vec.size(); i++){
            for(long long j = 0; j <= i; j++){
                fin.read((char*)&tmp_val, sizeof(double));
                if(tmp_val > cutoff){
                    group_map[i].insert(j);
                    group_map[j].insert(i);
                }
            }
        }
        fin.close();
    }else{
        ifstream fin(grm_file + ".sp.bin", std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".sp.bin");
            exit(1);
        }

        double tmp_val;
        long long index0, index1;
        while(fin.read((char*)&index0, sizeof(long long)) && fin.read((char*)&index1, sizeof(long long)) && fin.read((char*)&tmp_val, sizeof(double))){
            if(index0 != index1){
                if(tmp_val > cutoff){
                    group_map[index0].insert(index1);
                    group_map[index1].insert(index0);
                }
            }
        }
        fin.close();
    }

    spdlog::info("Group...");
    vector<set<long long>> group_vec;
    set<long long> id_index_in_grm_left_set = id_index_in_grm_set;
    while(id_index_in_grm_left_set.size() != 0){

        spdlog::info("  The proportion of Un-group samples: " + 
            std::to_string(id_index_in_grm_left_set.size() * 100.0 / num_id_in_grm) + "%");
        set<long long> id_index_in_group_set;
        id_index_in_group_set.insert(*id_index_in_grm_left_set.begin());
        vector<long long> id_index_in_group_added_vec = {*id_index_in_grm_left_set.begin()};

        while(1){
            long long size0 = id_index_in_group_set.size();
            set<long long> id_index_in_group_added_set;
            for(auto tmp:id_index_in_group_added_vec){
                id_index_in_group_added_set = set_union_(id_index_in_group_added_set, group_map[tmp]);
            }

            id_index_in_group_added_set = set_difference_(id_index_in_group_added_set, id_index_in_group_set);
            id_index_in_group_set = set_union_(id_index_in_group_added_set, id_index_in_group_set);

            id_index_in_group_added_vec.clear();
            id_index_in_group_added_vec.assign(id_index_in_group_added_set.begin(), id_index_in_group_added_set.end());
            long long size1 = id_index_in_group_set.size();
            spdlog::info("    Group size (increasing): " + std::to_string(id_index_in_group_set.size()));
            if(size0 == size1){
                break;
            }
        }
        id_index_in_grm_left_set = set_difference_(id_index_in_grm_left_set, id_index_in_group_set);
        group_vec.push_back(id_index_in_group_set);
    }

    spdlog::info("All iids are divided into " + std::to_string(group_vec.size()), 0);
    // string cutoff_str = std::to_string(cutoff);
    // cutoff_str.erase(cutoff_str.find_last_not_of('0') + 1, string::npos);
    // if(cutoff_str == "") cutoff_str = "0.0";
    ofstream fout;
    fout.open(grm_file + ".group.size" );
    if(!fout.is_open()){
        spdlog::error("Fail to open: " + out_file + ".group.size");
        exit(1);
    }
    map<long long, long long> id_index_group_map;
    map<long long, long long> group_size_map;
    for(long long i = 0; i < group_vec.size(); i++){
        set<long long> tmp_set = group_vec[i];
        // LOGGER.i("Group " + std::to_string(i+1) + ": Sample size is " + std::to_string(tmp_set.size()), 0);
        fout << i + 1 << " " << tmp_set.size() << endl;
        for(auto tmp:tmp_set){
            id_index_group_map[tmp] = i + 1;
            group_size_map[tmp] = tmp_set.size();
        }
    }
    fout.close();
    
    fout.open(grm_file + ".group");
    if(!fout.is_open()){
        spdlog::error("Fail to open: " + out_file + ".group");
        exit(1);
    }

    for(long long i = 0; i < id_in_grm_vec.size(); i++){
        fout << i + 1 << " " << id_in_grm_vec[i] << " " << id_index_group_map[i] << " " << group_size_map[i] << endl;
    }

    fout.close();

    
    if(sparse == 0){
        spdlog::info("Keep within-group elements and output the sparse GRM");
        ofstream fout;
        // fout.open("group.test.txt");
        map<string, int> group_index0_index1_map;
        for(long long i = 0; i < group_vec.size(); i++){
            set<long long> tmp_set = group_vec[i];
            for(auto tmp0:tmp_set){
                for(auto tmp1:tmp_set){
                    // fout << tmp0 << "_" << tmp1 << " ";
                    group_index0_index1_map[std::to_string(tmp0) + "_" + std::to_string(tmp1)] = 1;
                }
                fout << std::endl;
            }
        }
        // fout.close();

        ifstream fin(grm_file + ".bin", std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".bin");
            exit(1);
        }

        fout.open(grm_file + ".sp.bin", std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open: " + grm_file + ".sp.bin");
            exit(1);
        }

        double tmp_val;
        for(long long i = 0; i < id_in_grm_vec.size(); i++){
            for(long long j = 0; j <= i; j++){
                fin.read((char*)&tmp_val, sizeof(double));
                string tmp_str = std::to_string(i) + "_" + std::to_string(j);
                if(group_index0_index1_map.find(tmp_str) != group_index0_index1_map.end()){
                    fout.write((char*)&i, sizeof(long long));
                    fout.write((char*)&j, sizeof(long long));
                    fout.write((char*)&tmp_val, sizeof(double));
                }
            }
        }
        fin.close();
        fout.close();
    }
    
}


int ProcessGRM::run(int argc, char* argv[]){
    CLI::App app{"GMAT: GRM Processing Tool"};

    std::string description = R"(Quick start: 
    Merge GRM files
        gmat --process-grm --merge --grm test.agrm --npart 5
        gmat --process-grm --sparse --merge --grm test.agrm --npart 5

    PCA
        gmat --process-grm --pca --grm test.agrm --out test

    Calculate epistatic GRM
        gmat --process-grm --make-epis --grm test --epis-type aagrm --out test
        gmat --process-grm --make-epis --sparse --grm test --epis-type aagrm --out test

    Reformat GRM
        gmat --process-grm --reformat --grm test.agrm --out-fmt 0/1/2/3 --out test
        gmat --process-grm --reformat --sparse --grm test.agrm --out-fmt 0/1/2/3 --out test

    Reformat GRM from GCTA
        gmat --process-grm --reformat-gcta --grm test.grm --out test.agrm

    Check if the GRM is positive
        gmat --process-grm --check-positive --grm test.agrm
        gmat --process-grm --check-positive --grm test.agrm --sparse

    Calculate the inverse of GRM
        gmat --process-grm --make-inv --grm test.agrm --out-fmt 0/1/2/3 --out test.inv
        gmat --process-grm --make-inv --grm test.agrm --out-fmt 0/1/2/3 --out test.inv --sparse

    Remove the elements of given iids from GRM
        gmat --process-grm --remove-iid --grm test.agrm --iid id.txt --out test2.agrm
        gmat --process-grm --remove-iid --grm test.agrm --iid id.txt --out test2.agrm --sparse

    Output the one of the paired individuals with grm larger than the cut value
        gmat --process-grm --out-close --grm test.agrm --out test
        gmat --process-grm --out-close --grm test.agrm --out test --sparse

    Group iids according to the relationship
        gmat --process-grm --group --grm test.agrm --cut-value 0.05 --out test
        gmat --process-grm --group --grm test.agrm --cut-value 0.05 --out test --sparse
    )";

    app.description(description);
    
    bool process_grm = false;
    app.add_flag("--process-grm", process_grm, "Calculate GRM")->required();

    int threads = 10;
    app.add_option("--threads", threads, "Number of threads (default: 10)");

    string grm_file, iid_file, out_file, epis_type;
    int sparse = 0, out_fmt = 2;
    long long npart = 0;
    double val = 0, cut_value = 0.05;

    bool merge = false, pca = false, make_epis = false, reformat = false, 
         check_positive = false, make_inv = false, remove_iid = false, 
         out_close = false, keep_iid = false, reformat_gcta = false, group = false;

    app.add_flag("--merge", merge, "Merge GRM files");
    app.add_flag("--pca", pca, "Perform PCA");
    app.add_flag("--make-epis", make_epis, "Calculate epistatic GRM");
    app.add_flag("--reformat", reformat, "Reformat GRM");
    app.add_flag("--check-positive", check_positive, "Check if the GRM is positive");
    app.add_flag("--make-inv", make_inv, "Calculate the inverse of GRM");
    app.add_flag("--remove-iid", remove_iid, "Remove given iids from GRM");
    app.add_flag("--out-close", out_close, "Output pairs with GRM larger than cut value");
    app.add_flag("--keep-iid", keep_iid, "Keep given iids in GRM");
    app.add_flag("--reformat-gcta", reformat_gcta, "Reformat GCTA GRM");
    app.add_flag("--group", group, "Group individuals by GRM relationships");

    app.add_option("--grm", grm_file, "GRM file")->required();
    app.add_option("--iid", iid_file, "IID file for filtering");
    app.add_option("--out", out_file, "Output file prefix");
    app.add_option("--epis-type", epis_type, "Epistatic GRM type");
    app.add_option("--npart", npart, "Partition count (for merging)");
    app.add_option("--out-fmt", out_fmt, "Output format (0-3)");
    app.add_option("--val", val, "Diagonal value adjustment");
    app.add_option("--cut-value", cut_value, "Cutoff value for filtering");
    app.add_flag("--sparse", sparse, "Indicate sparse GRM");

    CLI11_PARSE(app, argc, argv);

    spdlog::info("Using {} threads", threads);
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);

    if (merge) {
        spdlog::info("Merging GRM files");
        if (npart <= 1) {
            spdlog::error("--npart must be greater than 1");
            exit(1);
        }
        this->merge_grm(grm_file, sparse, npart, out_file);
    } else if (pca) {
        spdlog::info("Performing PCA");
        Eigen::MatrixXd mat;
        std::map<std::string, long long> grm_id_map;
        if (!sparse) {
            this->read_grm_bin(grm_file, grm_id_map, mat);
        } else {
            Eigen::SparseMatrix<double> mat2;
            this->read_grm_sp_bin(grm_file, grm_id_map, mat2);
            mat = Eigen::MatrixXd(mat2);
        }
        this->pca(mat, out_file);
    } else if (make_epis) {
        spdlog::info("Calculating epistatic GRM");
        this->epistatic_grm(grm_file, sparse, epis_type, out_file);
    } else if (reformat) {
        spdlog::info("Reformatting GRM");
        vector<string> grm_id_vec = this->read_grm_id(grm_file);
        std::map<std::string, long long> grm_id_map;
        if (!sparse) {
            Eigen::MatrixXd mat;
            this->read_grm_bin(grm_file, grm_id_map, mat);
            this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
        } else {
            Eigen::SparseMatrix<double> mat;
            this->read_grm_sp_bin(grm_file, grm_id_map, mat);
            this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
        }
    } else if (check_positive) {
        spdlog::info("Checking if GRM is positive");
        Eigen::MatrixXd mat;
        std::map<std::string, long long> grm_id_map;
        this->read_grm_bin(grm_file, grm_id_map, mat);
        for (long long i = 0; i < mat.rows(); i++) {
            mat(i, i) += val;
        }
        CustomLLT LLTA;
        LLTA.compute(mat);
        spdlog::info("GRM is positive");
    } else if (make_inv) {
        spdlog::info("Computing GRM inverse");
        Eigen::MatrixXd mat;
        std::map<std::string, long long> grm_id_map;
        this->read_grm_bin(grm_file, grm_id_map, mat);
        CustomLLT LLTA;
        LLTA.compute(mat);
        mat = LLTA.inverse();
        vector<string> grm_id_vec = this->read_grm_id(grm_file);
        this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
    } else if(group){
        spdlog::info("Grouping individuals by GRM relationships");
        this->group(grm_file, sparse, cut_value, out_file);
    } else if (remove_iid) {
        spdlog::info("Removing specified iids from GRM");
        if (iid_file.empty()) {
            spdlog::error("--iid is required for removing iids");
            exit(1);
        }
        this->remove_id(grm_file, sparse, iid_file, out_file);
    } else if (keep_iid) {
        spdlog::info("Keeping specified iids in GRM");
        if (iid_file.empty()) {
            spdlog::error("--iid is required for keeping iids");
            exit(1);
        }
        this->keep_id(grm_file, sparse, iid_file, out_file);
    } else if (out_close) {
        spdlog::info("Outputting one of each pair of close individuals");
        this->remove_close_id(grm_file, sparse, cut_value, out_file);
    } else {
        spdlog::error("No valid operation specified");
        exit(1);
    }

    return 0;
}
