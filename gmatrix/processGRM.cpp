/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-08-12 17:23:31
 * LastEditors: Chao Ning
 * LastEditTime: 2026-04-12 19:22:43
 */

#include <cstdint>

#define EIGEN_USE_MKL_ALL  // !must be before include Eigen

#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <queue>
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
using std::getline;

using std::string;
using std::map;
using std::vector;
using std::set;

namespace {

std::string grm_prefix(const std::string& prefix) {
    constexpr const char* suffix = ".grm";
    const std::size_t suffix_len = std::char_traits<char>::length(suffix);
    if(prefix.size() >= suffix_len &&
       prefix.compare(prefix.size() - suffix_len, suffix_len, suffix) == 0){
        return prefix;
    }
    return prefix + suffix;
}

const char* dense_output_suffix(const int out_fmt) {
    switch (out_fmt) {
        case 0: return ".matrix";
        case 1: return ".index_triplet";
        case 2: return ".iid_triplet";
        case 3: return ".bin";
        default: return nullptr;
    }
}

const char* sparse_output_suffix(const int out_fmt) {
    switch (out_fmt) {
        case 0: return ".matrix";
        case 1: return ".index_triplet";
        case 2: return ".iid_triplet";
        case 3: return ".sp.bin";
        default: return nullptr;
    }
}

bool read_sparse_grm_record(ifstream& fin, std::int64_t& index0, std::int64_t& index1, double& value, const string& input_file) {
    if(!fin.read(reinterpret_cast<char*>(&index0), sizeof(std::int64_t))){
        if(fin.eof()){
            return false;
        }
        spdlog::error("Failed while reading: {}", input_file);
        exit(1);
    }

    if(!fin.read(reinterpret_cast<char*>(&index1), sizeof(std::int64_t)) ||
       !fin.read(reinterpret_cast<char*>(&value), sizeof(double))){
        spdlog::error("Failed while reading: {}", input_file);
        exit(1);
    }
    return true;
}

} // namespace

ProcessGRM::ProcessGRM(){
}

ProcessGRM::~ProcessGRM(){
}


/**
 * @brief Read sample IDs from `<prefix>.grm.id` in file order.
 *
 * Blank lines are ignored after trimming. The returned vector is validated to
 * ensure there are no duplicated IDs, because downstream matrix readers assume
 * a one-to-one mapping between ID order and matrix rows.
 */
vector<string> ProcessGRM::read_grm_id(const std::string& grm_file){
    const string grm_id_file = grm_prefix(grm_file) + ".id";
    ifstream fin(grm_id_file);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}", grm_id_file);
        exit(EXIT_FAILURE);
    }

    vector<string> grm_id_vec;
    string line;
    while(getline(fin, line)){
        process_line(line);
        if(line.empty()){
            continue;
        }
        grm_id_vec.push_back(line);
    }

    if(is_duplicated(grm_id_vec)){
        spdlog::error("Duplicated iids exist in the grm file");
        exit(EXIT_FAILURE);
    }
    return grm_id_vec;
}


/**
 * @brief Write sample IDs to `<out>.grm.id` in GRM row order.
 *
 * The file is emitted as one ID per line so that downstream readers can reuse
 * the same row/column ordering without additional parsing logic.
 */
void ProcessGRM::out_grm_id(const std::string& out_file, const std::vector<std::string>& grm_id_vec){
    const string out_id_file = grm_prefix(out_file) + ".id";
    ofstream fout(out_id_file);
    if(!fout.is_open()){
        spdlog::error("Fail to open: {}", out_id_file);
        exit(EXIT_FAILURE);
    }

    for(const auto& grm_id : grm_id_vec){
        fout << grm_id << '\n';
    }

    if(!fout){
        spdlog::error("Failed while writing: {}", out_id_file);
        exit(EXIT_FAILURE);
    }
}


/**
 * @brief Read a dense lower-triangular GRM and materialize it as a full matrix.
 *
 * The on-disk `.grm.bin` file stores doubles in lower-triangular order. When
 * `grm_id_map` is empty, the matrix is read in its native order. Otherwise the
 * GRM is subset and reordered to match the caller-provided ID map.
 */
void ProcessGRM::read_grm_bin(const std::string& grm_file, const std::map<std::string, std::int64_t>& grm_id_map, Eigen::MatrixXd& mat, double val){
    const vector<string> grm_id_vec = this->read_grm_id(grm_file);
    vector<std::int64_t> grm_id_index_vec;
    grm_id_index_vec.reserve(grm_id_vec.size());
    map<string, std::int64_t> working_grm_id_map = grm_id_map;

    std::int64_t num_id_used = 0;
    if(working_grm_id_map.empty()){
        // Read the matrix in its original row/column order.
        num_id_used = static_cast<std::int64_t>(grm_id_vec.size());
        for(std::int64_t i = 0; i < num_id_used; ++i){
            working_grm_id_map[grm_id_vec[i]] = i;
            grm_id_index_vec.push_back(i);
        }
    }else{
        // Preserve the caller's requested order while skipping IDs not present
        // in the current GRM file.
        for(const auto& grm_id : grm_id_vec){
            const auto it = working_grm_id_map.find(grm_id);
            if(it != working_grm_id_map.end()){
                grm_id_index_vec.push_back(it->second);
                ++num_id_used;
            }else{
                grm_id_index_vec.push_back(-1);
            }
        }
    }

    if(working_grm_id_map.size() != static_cast<std::size_t>(num_id_used)){
        spdlog::error("iids in the grm do not match the data");
        exit(1);
    }

    const string grm_bin_file = grm_prefix(grm_file) + ".bin";
    ifstream fin(grm_bin_file, std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}", grm_bin_file);
        exit(1);
    }

    mat.setZero(num_id_used, num_id_used);
    const std::int64_t num_id = static_cast<std::int64_t>(grm_id_index_vec.size());
    for(std::int64_t i = 0; i < num_id; ++i){
        const std::int64_t index0 = grm_id_index_vec[i];
        for(std::int64_t j = 0; j <= i; ++j){
            double tmp_val = 0.0;
            fin.read(reinterpret_cast<char*>(&tmp_val), sizeof(double));
            if(!fin){
                spdlog::error("Failed while reading dense GRM file: {}", grm_bin_file);
                exit(1);
            }

            const std::int64_t index1 = grm_id_index_vec[j];
            if(index0 != -1 && index1 != -1){
                if(index0 == index1){
                    tmp_val += val;
                }
                mat(index0, index1) = mat(index1, index0) = tmp_val;
            }
        }
    }
}


/**
 * @brief Read a sparse triplet GRM and materialize it as a symmetric sparse matrix.
 *
 * The `.grm.sp.bin` payload stores lower-triangular triplets `(row, col, value)`
 * using zero-based indices. As with `read_grm_bin`, an empty ID map preserves
 * the native order, while a non-empty map subsets and reorders the matrix.
 */
void ProcessGRM::read_grm_sp_bin(const std::string& grm_sparse_file, const std::map<std::string, std::int64_t>& grm_id_map, Eigen::SparseMatrix<double>& mat, double val){
    const vector<string> grm_id_vec = this->read_grm_id(grm_sparse_file);
    vector<std::int64_t> grm_id_index_vec;
    grm_id_index_vec.reserve(grm_id_vec.size());
    map<string, std::int64_t> working_grm_id_map = grm_id_map;
    std::int64_t num_id_used = 0;

    if(working_grm_id_map.empty()){
        num_id_used = static_cast<std::int64_t>(grm_id_vec.size());
        for(std::int64_t i = 0; i < num_id_used; ++i){
            working_grm_id_map[grm_id_vec[i]] = i;
            grm_id_index_vec.push_back(i);
        }
    }else{
        for(const auto& grm_id : grm_id_vec){
            const auto it = working_grm_id_map.find(grm_id);
            if(it != working_grm_id_map.end()){
                grm_id_index_vec.push_back(it->second);
                ++num_id_used;
            }else{
                grm_id_index_vec.push_back(-1);
            }
        }
    }
    
    if(working_grm_id_map.size() != static_cast<std::size_t>(num_id_used)){
        spdlog::error("iids in the grm do not match the data");
        exit(1);
    }

    const string grm_sp_bin_file = grm_prefix(grm_sparse_file) + ".sp.bin";
    ifstream fin(grm_sp_bin_file, std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}", grm_sp_bin_file);
        exit(1);
    }


    std::vector<Eigen::Triplet<double>> triplet_list;
    double tmp_val = 0.0;
    std::int64_t index0 = 0, index1 = 0;
    while(read_sparse_grm_record(fin, index0, index1, tmp_val, grm_sp_bin_file)){
        const std::int64_t index_curr0 = grm_id_index_vec[index0];
        const std::int64_t index_curr1 = grm_id_index_vec[index1];
        if(index_curr0 != -1 && index_curr1 != -1){
            if(index_curr0 == index_curr1){
                triplet_list.emplace_back(index_curr0, index_curr1, tmp_val + val);
            }else{
                triplet_list.emplace_back(index_curr0, index_curr1, tmp_val);
                triplet_list.emplace_back(index_curr1, index_curr0, tmp_val);
            }
        }
    }
    
    mat.resize(num_id_used, num_id_used);
    mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
}


/**
 * @brief Merge partitioned GRM sidecar files back into a single GRM prefix.
 *
 * Partitioned GRM runs store `.id`, `.N`, and dense binary payloads with the
 * suffix `.<npart>_<ipart>`. This helper reconstructs the original outputs in
 * the same order they were produced.
 */
void ProcessGRM::merge_grm(const std::string& grm_file, int npart, const std::string& out_file){
    const auto build_part_files = [&](const std::string& suffix) {
        vector<string> file_vec;
        file_vec.reserve(static_cast<std::size_t>(npart));
        const string part_prefix = grm_file;
        for(int i = 1; i <= npart; ++i){
            file_vec.push_back(part_prefix + "." + std::to_string(npart) + "_" + std::to_string(i)+ ".grm" + suffix);
        }
        return file_vec;
    };

    // All `.id` fragments are concatenated because each partition contributes a
    // disjoint row block.
    spdlog::info("Merging *.id files");
    const vector<string> id_files = build_part_files(".id");
    ProcessGRM::merge_file(id_files, grm_prefix(out_file) + ".id", false);

    // `.N` is metadata shared across partitions, so one representative file is
    // enough for the merged output.
    spdlog::info("Merging *.N file");
    const vector<string> n_files = build_part_files(".N");
    ProcessGRM::merge_file({n_files.front()}, grm_prefix(out_file) + ".N", false);

    // The dense binary payload is appended in partition order to recover the
    // original lower-triangular serialization.
    const vector<string> dense_bin_files = build_part_files(".bin");
    if(!ifstream(dense_bin_files.front(), std::ios::binary).good()){
        spdlog::error("Failed to find partitioned binary GRM files for prefix {}", grm_file);
        exit(1);
    }

    spdlog::info("Merging *.bin files");
    ProcessGRM::merge_file(dense_bin_files, grm_prefix(out_file) + ".bin", true);
}

/**
 * @brief Concatenate multiple files into one output file.
 *
 * The merge order follows `file_vec` exactly. Binary files are copied byte by
 * byte to preserve the existing behavior, while text files are appended
 * line-by-line with trailing newlines.
 */
void ProcessGRM::merge_file(const vector<string>& file_vec, const std::string& out_file, bool bin){
    if(bin){
        ofstream fout(out_file, std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the: {}", out_file);
            exit(1);
        }

        for(const auto& input_file : file_vec){
            ifstream fin(input_file, std::ios::binary);
            if(!fin.is_open()){
                spdlog::error("Fail to open the: {}", input_file);
                exit(1);
            }
            spdlog::info("Merging file: {}", input_file);

            char tmp = '\0';
            while(fin.read(&tmp, sizeof(char))){
                fout.write(&tmp, sizeof(char));
            }

            if(!fin.eof()){
                spdlog::error("Failed while reading: {}", input_file);
                exit(1);
            }
        }

        if(!fout){
            spdlog::error("Failed while writing: {}", out_file);
            exit(1);
        }
        return;
    }

    ofstream fout(out_file);
    if(!fout.is_open()){
        spdlog::error("Fail to open the: {}", out_file);
        exit(1);
    }

    for(const auto& input_file : file_vec){
        ifstream fin(input_file);
        if(!fin.is_open()){
            spdlog::error("Fail to open the: {}", input_file);
            exit(1);
        }
        spdlog::info("Merging file: {}", input_file);

        string line;
        while(std::getline(fin, line)){
            fout << line << '\n';
        }

        if(!fin.eof()){
            spdlog::error("Failed while reading: {}", input_file);
            exit(1);
        }
    }

    if(!fout){
        spdlog::error("Failed while writing: {}", out_file);
        exit(1);
    }
}



void ProcessGRM::out_gmat(const Eigen::MatrixXd& gmat, const vector<string>& grm_id_vec, int out_fmt, const string& out_file){
    const std::int64_t num_id = static_cast<std::int64_t>(grm_id_vec.size());
    const char* suffix = dense_output_suffix(out_fmt);
    const string out_grm_prefix = grm_prefix(out_file);
    if(suffix == nullptr){
        spdlog::error("Unsupported output format: --out-fmt {}", out_fmt);
        exit(1);
    }

    if(out_fmt == 0){ // matrix
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        for(std::int64_t i = 0; i < num_id; ++i){
            for(std::int64_t j = 0; j < num_id; ++j){
                fout << gmat(i, j) << " ";
            }
            fout << '\n';
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 1){ // index_triplet
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        for(std::int64_t i = 0; i < num_id; ++i){
            for(std::int64_t j = 0; j <= i; ++j){
                fout << i + 1 << " " << j + 1 << " " << gmat(i, j) << '\n';
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 2){ // iid_triplet
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        for(std::int64_t i = 0; i < num_id; ++i){
            for(std::int64_t j = 0; j <= i; ++j){
                fout << grm_id_vec[i] << " " << grm_id_vec[j] << " " << gmat(i, j) << '\n';
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 3){ // .bin
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file, std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }

        for(std::int64_t i = 0; i < num_id; ++i){
            for(std::int64_t j = 0; j <= i; ++j){
                const double value = gmat(i, j);
                fout.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }
}


void ProcessGRM::out_gmat(const Eigen::SparseMatrix<double>& gmat, const vector<string>& grm_id_vec, int out_fmt, const std::string& out_file){
    const std::int64_t num_id = static_cast<std::int64_t>(grm_id_vec.size());
    const Eigen::SparseMatrix<double> upper_gmat = gmat.triangularView<Eigen::Upper>();
    const char* suffix = sparse_output_suffix(out_fmt);
    const string out_grm_prefix = grm_prefix(out_file);
    if(suffix == nullptr){
        spdlog::error("Unsupported output format: --out-fmt {}", out_fmt);
        exit(1);
    }

    if(out_fmt == 0){ // matrix
        spdlog::warn("Output sparse relationship matrix in dense format");
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        fout << Eigen::MatrixXd(gmat) << '\n';
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 1){ // index_triplet
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        for (std::int64_t k = 0; k < upper_gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(upper_gmat, k); it; ++it){
                fout << it.col() + 1 << " " << it.row() + 1 << " " << it.value() << '\n'; // output the transpose of upper triangular part
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 2){ // iid_triplet
        const string out_gmat_file = out_grm_prefix + suffix;
        ofstream fout(out_gmat_file);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }
        
        for (std::int64_t k = 0; k < upper_gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(upper_gmat, k); it; ++it){
                fout << grm_id_vec[it.col()] << " " << grm_id_vec[it.row()] << " " << it.value() << '\n'; // output the transpose of upper triangular part
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }else if(out_fmt == 3){
        const string out_gmat_file = out_grm_prefix + suffix;
        std::ofstream fout(out_gmat_file, std::ios::binary);
        if(!fout.is_open()){
            spdlog::error("Fail to open the output file: {}", out_gmat_file);
            exit(1);
        }

        for (std::int64_t k = 0; k < upper_gmat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(upper_gmat, k); it; ++it){
                //output the transpose of upper triangular part -> low triangular part
                const std::int64_t index0 = it.col();
                const std::int64_t index1 = it.row();
                const double value = it.value();
                fout.write(reinterpret_cast<const char*>(&index0), sizeof(std::int64_t));
                fout.write(reinterpret_cast<const char*>(&index1), sizeof(std::int64_t));
                fout.write(reinterpret_cast<const char*>(&value), sizeof(double));
            }
        }
        if(!fout){
            spdlog::error("Failed while writing the output file: {}", out_gmat_file);
            exit(1);
        }
    }
}



/**
 * @brief Compute PCA from a dense GRM and write eigenvalues/eigenvectors.
 *
 * Eigen returns eigenpairs in ascending order for self-adjoint matrices, so the
 * output is written in reverse order to keep the leading principal components
 * first in both files.
 */
void ProcessGRM::pca(const Eigen::MatrixXd& gmat, const std::string& out_file){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(gmat);
    if (eigensolver.info() != Eigen::Success) {
        spdlog::error("Fail to compute eigenvalues and eigenvectors");
        exit(1);
    }

    const string out_grm_prefix = grm_prefix(out_file);
    const string eigenvals_file = out_grm_prefix + ".eigenvals";
    const string eigenvecs_file = out_grm_prefix + ".eigenvecs";
    ofstream fout_eigenvals(eigenvals_file);
    ofstream fout_eigenvecs(eigenvecs_file);
    if((!fout_eigenvals.is_open()) || (!fout_eigenvecs.is_open())){
        spdlog::error("Fail to open the output file: {}, {}", eigenvals_file, eigenvecs_file);
        exit(1);
    }

    const auto& eigenvals = eigensolver.eigenvalues();
    const auto& eigenvecs = eigensolver.eigenvectors();
    const Eigen::Index num_element = eigenvals.size();

    for(Eigen::Index i = num_element - 1; i >= 0; --i){
        fout_eigenvals << eigenvals(i) << '\n';
    }

    for(Eigen::Index i = 0; i < eigenvecs.rows(); ++i){
        for(Eigen::Index j = num_element - 1; j > 0; --j){
            fout_eigenvecs << eigenvecs(i, j) << " ";
        }
        fout_eigenvecs << eigenvecs(i, 0) << '\n';
    }

    if((!fout_eigenvals) || (!fout_eigenvecs)){
        spdlog::error("Failed while writing PCA outputs: {}, {}", eigenvals_file, eigenvecs_file);
        exit(1);
    }
}

/**
 * @brief Multiply two dense GRMs elementwise to form an epistatic GRM.
 *
 * The two input prefixes must share the same `.grm.id` ordering. The output
 * keeps that ID order, concatenates the two `.grm.N` metadata files, and
 * writes the lower-triangular dense payload to `<out>.grm.bin`.
 */
void ProcessGRM::epistatic_grm(const std::string& prefix1, const std::string& prefix2, const std::string& out_file){
    const vector<string> grm_id_vec1 = ProcessGRM::read_grm_id(prefix1);
    const vector<string> grm_id_vec2 = ProcessGRM::read_grm_id(prefix2);
    if(grm_id_vec1 != grm_id_vec2){
        spdlog::error("The ID order in {}.id and {}.id does not match", prefix1, prefix2);
        exit(1);
    }

    // The output GRM reuses the first ID file after verifying that both inputs
    // describe the same samples in the same order.
    ProcessGRM::merge_file({grm_prefix(prefix1) + ".id"}, grm_prefix(out_file) + ".id", false);

    // Preserve both input metadata files in order for downstream inspection.
    ProcessGRM::merge_file({grm_prefix(prefix1) + ".N", grm_prefix(prefix2) + ".N"}, grm_prefix(out_file) + ".N", false);

    const string bin_file1 = grm_prefix(prefix1) + ".bin";
    const string bin_file2 = grm_prefix(prefix2) + ".bin";
    const string out_bin_file = grm_prefix(out_file) + ".bin";
    ifstream fin1(bin_file1, std::ios::binary);
    ifstream fin2(bin_file2, std::ios::binary);
    ofstream fout(out_bin_file, std::ios::binary);
    if(!fin1.is_open()){
        spdlog::error("Fail to open: {}", bin_file1);
        exit(1);
    }
    if(!fin2.is_open()){
        spdlog::error("Fail to open: {}", bin_file2);
        exit(1);
    }
    if(!fout.is_open()){
        spdlog::error("Fail to open: {}", out_bin_file);
        exit(1);
    }

    double val1 = 0.0;
    double val2 = 0.0;
    while(fin1.read(reinterpret_cast<char*>(&val1), sizeof(double))){
        if(!fin2.read(reinterpret_cast<char*>(&val2), sizeof(double))){
            spdlog::error("Dense GRM files have inconsistent lengths: {} and {}", bin_file1, bin_file2);
            exit(1);
        }
        const double product = val1 * val2;
        fout.write(reinterpret_cast<const char*>(&product), sizeof(double));
        if(!fout){
            spdlog::error("Failed while writing: {}", out_bin_file);
            exit(1);
        }
    }

    if(!fin1.eof()){
        spdlog::error("Failed while reading: {}", bin_file1);
        exit(1);
    }

    if(fin2.read(reinterpret_cast<char*>(&val2), sizeof(double))){
        spdlog::error("Dense GRM files have inconsistent lengths: {} and {}", bin_file1, bin_file2);
        exit(1);
    }
    if(!fin2.eof()){
        spdlog::error("Failed while reading: {}", bin_file2);
        exit(1);
    }
}

/**
 * @brief Greedily remove one IID from each set of overly related pairs.
 *
 * A graph edge is added whenever the pairwise GRM exceeds `cutoff`. The
 * routine repeatedly removes the IID with the largest current degree until no
 * close pairs remain, then writes the pruned IDs to `<out>.grm.pruned_iid`.
 */
void ProcessGRM::remove_close_id(const std::string& grm_file, int sparse, double cutoff, const std::string& out_file){
    const vector<string> id_in_grm_vec = this->read_grm_id(grm_file);
    const std::int64_t num_id = static_cast<std::int64_t>(id_in_grm_vec.size());
    map<string, set<string>> close_id_map;
    std::int64_t num_close = 0;

    for(const auto& id : id_in_grm_vec){
        close_id_map.emplace(id, set<string>{});
    }

    const auto prune_isolated_ids = [&close_id_map]() {
        vector<string> isolated_id_vec;
        isolated_id_vec.reserve(close_id_map.size());
        for(const auto& [id, close_ids] : close_id_map){
            if(close_ids.empty()){
                isolated_id_vec.push_back(id);
            }
        }
        for(const auto& id : isolated_id_vec){
            close_id_map.erase(id);
        }
    };

    if(sparse == 0){
        const string grm_bin_file = grm_prefix(grm_file) + ".bin";
        ifstream fin(grm_bin_file, std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: {}", grm_bin_file);
            exit(1);
        }

        double tmp_val = 0.0;
        for(std::int64_t i = 0; i < num_id; ++i){
            for(std::int64_t j = 0; j <= i; ++j){
                if(!fin.read(reinterpret_cast<char*>(&tmp_val), sizeof(double))){
                    spdlog::error("Failed while reading: {}", grm_bin_file);
                    exit(1);
                }
                if(tmp_val > cutoff && i != j){
                    ++num_close;
                    close_id_map[id_in_grm_vec[i]].insert(id_in_grm_vec[j]);
                    close_id_map[id_in_grm_vec[j]].insert(id_in_grm_vec[i]);
                }
            }
        }
    }else{
        const string grm_sp_bin_file = grm_prefix(grm_file) + ".sp.bin";
        ifstream fin(grm_sp_bin_file, std::ios::binary);
        if(!fin.is_open()){
            spdlog::error("Fail to open: {}", grm_sp_bin_file);
            exit(1);
        }

        double tmp_val = 0.0;
        std::int64_t index0 = 0, index1 = 0;
        while(read_sparse_grm_record(fin, index0, index1, tmp_val, grm_sp_bin_file)){
            if(tmp_val > cutoff && index0 != index1){
                ++num_close;
                close_id_map[id_in_grm_vec[index0]].insert(id_in_grm_vec[index1]);
                close_id_map[id_in_grm_vec[index1]].insert(id_in_grm_vec[index0]);
            }
        }
    }

    spdlog::info("There are {} pair of close individuals", num_close);

    prune_isolated_ids();

    vector<string> remove_id_vec;
    while(!close_id_map.empty()){
        const auto closest_it = std::max_element(
            close_id_map.begin(),
            close_id_map.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.second.size() < rhs.second.size();
            });

        const string closest_id = closest_it->first;
        remove_id_vec.push_back(closest_id);

        const set<string> removed_set = closest_it->second;
        close_id_map.erase(closest_it);

        for(const auto& removed_id : removed_set){
            close_id_map[removed_id].erase(closest_id);
        }

        prune_isolated_ids();
    }

    const string pruned_iid_file = grm_prefix(out_file) + ".pruned_iid";
    ofstream fout(pruned_iid_file);
    if(!fout.is_open()){
        spdlog::error("Fail to open: {}", pruned_iid_file);
        exit(1);
    }

    for(const auto& removed_id : remove_id_vec){
        fout << removed_id << '\n';
    }
    if(!fout){
        spdlog::error("Failed while writing: {}", pruned_iid_file);
        exit(1);
    }
}

/**
 * @brief Group individuals by connected components in the thresholded GRM graph.
 *
 * Every IID is a graph node and edges are created for dense GRM values above
 * `cutoff`. The output is written with the input GRM prefix and includes
 * per-group sizes, a per-IID group assignment, and a derived sparse GRM that
 * keeps only within-group pairs.
 */
void ProcessGRM::group_related_samples(const std::string& grm_file, double cutoff){
    const vector<string> id_in_grm_vec = this->read_grm_id(grm_file);
    const std::int64_t num_id_in_grm = static_cast<std::int64_t>(id_in_grm_vec.size());
    vector<vector<std::int64_t>> adjacency_list(num_id_in_grm);
    std::int64_t num_edges_above_cutoff = 0;

    spdlog::info(
        "Grouping {} samples from dense GRM input with cutoff {}",
        num_id_in_grm,
        cutoff);

    if(num_id_in_grm == 0){
        spdlog::warn("No sample IDs were found in {}.grm.id", grm_file);
        return;
    }

    const string grm_bin_file = grm_prefix(grm_file) + ".bin";
    ifstream fin(grm_bin_file, std::ios::binary);
    if(!fin.is_open()){
        spdlog::error("Fail to open: {}", grm_bin_file);
        exit(1);
    }

    double tmp_val = 0.0;
    for(std::int64_t i = 0; i < num_id_in_grm; ++i){
        for(std::int64_t j = 0; j <= i; ++j){
            if(!fin.read(reinterpret_cast<char*>(&tmp_val), sizeof(double))){
                spdlog::error("Failed while reading: {}", grm_bin_file);
                exit(1);
            }
            if(i != j && tmp_val > cutoff){
                adjacency_list[i].push_back(j);
                adjacency_list[j].push_back(i);
                ++num_edges_above_cutoff;
            }
        }
    }

    spdlog::info("Built relatedness graph with {} edges above cutoff", num_edges_above_cutoff);

    vector<std::int64_t> id_to_group(num_id_in_grm, -1);
    vector<std::int64_t> group_sizes;
    std::int64_t grouped_count = 0;
    std::int64_t num_groups = 0;

    for(std::int64_t seed_index = 0; seed_index < num_id_in_grm; ++seed_index){
        if(id_to_group[seed_index] != -1){
            continue;
        }

        std::queue<std::int64_t> frontier;
        frontier.push(seed_index);
        id_to_group[seed_index] = num_groups;
        std::int64_t current_group_size = 0;

        // Breadth-first expansion finds one connected component without the
        // repeated set copies used by the previous implementation.
        while(!frontier.empty()){
            const std::int64_t current_index = frontier.front();
            frontier.pop();
            ++current_group_size;

            for(const auto neighbor_index : adjacency_list[current_index]){
                if(id_to_group[neighbor_index] == -1){
                    id_to_group[neighbor_index] = num_groups;
                    frontier.push(neighbor_index);
                }
            }
        }

        group_sizes.push_back(current_group_size);
        grouped_count += current_group_size;
        ++num_groups;

        spdlog::info(
            "Identified group {} with {} samples ({} / {} assigned)",
            num_groups,
            current_group_size,
            grouped_count,
            num_id_in_grm);
    }

    spdlog::info("All iids are divided into {} groups", num_groups);

    const string group_size_file = grm_prefix(grm_file) + ".group.size";
    ofstream group_size_out(group_size_file);
    if(!group_size_out.is_open()){
        spdlog::error("Fail to open: {}", group_size_file);
        exit(1);
    }

    for(std::int64_t group_index = 0; group_index < num_groups; ++group_index){
        group_size_out << group_index + 1 << " " << group_sizes[group_index] << '\n';
    }
    if(!group_size_out){
        spdlog::error("Failed while writing: {}", group_size_file);
        exit(1);
    }

    const string group_file = grm_prefix(grm_file) + ".group";
    ofstream group_out(group_file);
    if(!group_out.is_open()){
        spdlog::error("Fail to open: {}", group_file);
        exit(1);
    }

    for(std::int64_t i = 0; i < num_id_in_grm; ++i){
        const std::int64_t group_index = id_to_group[i];
        group_out << i + 1 << " " << id_in_grm_vec[i] << " "
                  << group_index + 1 << " " << group_sizes[group_index] << '\n';
    }
    if(!group_out){
        spdlog::error("Failed while writing: {}", group_file);
        exit(1);
    }

    spdlog::info("Writing sparse GRM that keeps only within-group entries");

    ifstream dense_grm_in(grm_bin_file, std::ios::binary);
    if(!dense_grm_in.is_open()){
        spdlog::error("Fail to open: {}", grm_bin_file);
        exit(1);
    }

    const string grm_sp_bin_file = grm_prefix(grm_file) + ".sp.bin";
    ofstream fout(grm_sp_bin_file, std::ios::binary);
    if(!fout.is_open()){
        spdlog::error("Fail to open: {}", grm_sp_bin_file);
        exit(1);
    }

    tmp_val = 0.0;
    std::int64_t num_values_kept = 0;
    for(std::int64_t i = 0; i < num_id_in_grm; ++i){
        for(std::int64_t j = 0; j <= i; ++j){
            if(!dense_grm_in.read(reinterpret_cast<char*>(&tmp_val), sizeof(double))){
                spdlog::error("Failed while reading: {}", grm_bin_file);
                exit(1);
            }
            if(id_to_group[i] == id_to_group[j]){
                fout.write(reinterpret_cast<const char*>(&i), sizeof(std::int64_t));
                fout.write(reinterpret_cast<const char*>(&j), sizeof(std::int64_t));
                fout.write(reinterpret_cast<const char*>(&tmp_val), sizeof(double));
                ++num_values_kept;
            }
        }
    }
    if(!fout){
        spdlog::error("Failed while writing: {}", grm_sp_bin_file);
        exit(1);
    }

    spdlog::info("Wrote {} within-group lower-triangular values to {}", num_values_kept, grm_sp_bin_file);
}


int ProcessGRM::run(int argc, char* argv[]){
    CLI::App app{"GMAT: GRM Processing Tool"};

    // Keep the help text task-oriented: each block shows the primary flag plus
    // the minimal extra inputs typically needed to run that subcommand.
    app.description(R"(Quick start:
    Merge partitioned GRMs
        gmat --process-grm --merge --grm test --npart 5

    PCA
        gmat --process-grm --pca --grm test --out test

    Epistatic GRM from two dense GRMs
        gmat --process-grm --make-epis --prefix1 test1 --prefix2 test2 --out test

    Reformat GRM
        gmat --process-grm --reformat --grm test --out-fmt 0 --out test         # test.grm.matrix
        gmat --process-grm --reformat --grm test --out-fmt 1 --out test         # test.grm.index_triplet
        gmat --process-grm --reformat --grm test --out-fmt 2 --out test         # test.grm.iid_triplet
        gmat --process-grm --reformat --grm test --out-fmt 3 --out test         # test.grm.bin
        gmat --process-grm --reformat --sparse --grm test --out-fmt 0 --out test # test.grm.matrix
        gmat --process-grm --reformat --sparse --grm test --out-fmt 1 --out test # test.grm.index_triplet
        gmat --process-grm --reformat --sparse --grm test --out-fmt 2 --out test # test.grm.iid_triplet
        gmat --process-grm --reformat --sparse --grm test --out-fmt 3 --out test # test.grm.sp.bin

    GRM inverse
        gmat --process-grm --make-inv --grm test --out-fmt 0 --out test.inv     # test.inv.grm.matrix
        gmat --process-grm --make-inv --grm test --out-fmt 1 --out test.inv     # test.inv.grm.index_triplet
        gmat --process-grm --make-inv --grm test --out-fmt 2 --out test.inv     # test.inv.grm.iid_triplet
        gmat --process-grm --make-inv --grm test --out-fmt 3 --out test.inv     # test.inv.grm.bin

    Output one IID from each close pair
        gmat --process-grm --out-close --grm test --out test

    Group related IIDs
        gmat --process-grm --group --grm test --cut-value 0.05
    )");
    
    // `--process-grm` selects this top-level tool, while the remaining boolean
    // flags choose the concrete operation to run.
    bool process_grm = false;
    bool merge = false, pca = false, make_epis = false, reformat = false,
         make_inv = false, out_close = false, group = false;

    // Runtime controls and shared command inputs.
    int threads = 10;
    string grm_file, out_file, prefix1, prefix2;
    int sparse = 0, out_fmt = 2;
    std::int64_t npart = 0;
    double cut_value = 0.05;

    app.add_flag("--process-grm", process_grm, "Process an existing GRM")->required();
    app.add_option("--threads", threads, "Number of threads (default: 10)");

    // Task-selection flags. Exactly one of these should describe the requested
    // action; dispatch happens after argument parsing.
    app.add_flag("--merge", merge, "Merge GRM files");
    app.add_flag("--pca", pca, "Perform PCA");
    app.add_flag("--make-epis", make_epis, "Calculate epistatic GRM");
    app.add_flag("--reformat", reformat, "Reformat GRM");
    app.add_flag("--make-inv", make_inv, "Calculate the inverse of GRM");
    app.add_flag("--out-close", out_close, "Output pairs with GRM larger than cut value");
    app.add_flag("--group", group, "Group individuals by GRM relationships");

    // Shared file and numeric options. Individual subcommands only consume the
    // subset that is relevant to them.
    app.add_option("--grm", grm_file, "GRM prefix without the .grm suffix");
    app.add_option("--out", out_file, "Output prefix without the .grm suffix");
    app.add_option("--prefix1", prefix1, "First dense GRM prefix without the .grm suffix for --make-epis");
    app.add_option("--prefix2", prefix2, "Second dense GRM prefix without the .grm suffix for --make-epis");
    app.add_option("--npart", npart, "Partition count (for merging)");
    app.add_option("--out-fmt", out_fmt, "Output format: 0=matrix, 1=index_triplet, 2=iid_triplet, 3=binary");
    app.add_option("--cut-value", cut_value, "Cutoff value for filtering");
    app.add_flag("--sparse", sparse, "Indicate sparse GRM");

    CLI11_PARSE(app, argc, argv);

    spdlog::info("Using {} threads", threads);
    mkl_set_num_threads(threads);
    omp_set_num_threads(threads);

    const bool sparse_input = (sparse != 0);

    // Small local helpers keep the command branches focused on their own logic
    // instead of repeating the same input-validation and loading boilerplate.
    const auto require_grm_file = [&]() {
        if (grm_file.empty()) {
            spdlog::error("--grm is required for this operation");
            std::exit(1);
        }
    };

    const auto require_prefix_pair = [&]() {
        if (prefix1.empty() || prefix2.empty()) {
            spdlog::error("--prefix1 and --prefix2 are required for --make-epis");
            std::exit(1);
        }
    };

    const auto load_dense_grm = [&](bool allow_sparse_input) {
        Eigen::MatrixXd mat;
        std::map<std::string, std::int64_t> grm_id_map;
        if (!sparse_input || !allow_sparse_input) {
            this->read_grm_bin(grm_file, grm_id_map, mat);
        } else {
            Eigen::SparseMatrix<double> sparse_mat;
            this->read_grm_sp_bin(grm_file, grm_id_map, sparse_mat);
            mat = Eigen::MatrixXd(sparse_mat);
        }
        return mat;
    };

    // Merge is special because it operates directly on partitioned files instead
    // of materializing a matrix first.
    if (merge) {
        spdlog::info("Merging GRM files");
        require_grm_file();
        if (npart <= 1) {
            spdlog::error("--npart must be greater than 1");
            std::exit(1);
        }
        this->merge_grm(grm_file, npart, out_file);
        return 0;
    }

    // PCA always runs on a dense matrix. Sparse input is therefore expanded
    // before decomposition when `--sparse` is set.
    if (pca) {
        spdlog::info("Performing PCA");
        require_grm_file();
        Eigen::MatrixXd mat = load_dense_grm(true);
        this->pca(mat, out_file);
        return 0;
    }

    // Epistatic GRM construction multiplies two dense lower-triangular GRM
    // payloads element by element after verifying that their ID order matches.
    if (make_epis) {
        spdlog::info("Calculating epistatic GRM");
        require_prefix_pair();
        if (sparse_input) {
            spdlog::error("--make-epis no longer supports --sparse");
            std::exit(1);
        }
        this->epistatic_grm(prefix1, prefix2, out_file);
        return 0;
    }

    // Reformat preserves the original storage mode: dense input stays dense,
    // sparse input stays sparse, and both paths reuse the same ID ordering.
    if (reformat) {
        spdlog::info("Reformatting GRM");
        require_grm_file();
        const vector<string> grm_id_vec = this->read_grm_id(grm_file);
        std::map<std::string, std::int64_t> grm_id_map;
        if (!sparse_input) {
            Eigen::MatrixXd mat;
            this->read_grm_bin(grm_file, grm_id_map, mat);
            this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
        } else {
            Eigen::SparseMatrix<double> mat;
            this->read_grm_sp_bin(grm_file, grm_id_map, mat);
            this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
        }
        this->out_grm_id(out_file, grm_id_vec);
        return 0;
    }

    // Inversion also goes through the dense representation because the current
    // linear-algebra helper expects a dense positive-definite matrix.
    if (make_inv) {
        spdlog::info("Computing GRM inverse");
        require_grm_file();
        Eigen::MatrixXd mat = load_dense_grm(false);
        CustomLLT LLTA;
        LLTA.compute(mat);
        mat = LLTA.inverse();
        const vector<string> grm_id_vec = this->read_grm_id(grm_file);
        this->out_gmat(mat, grm_id_vec, out_fmt, out_file);
        return 0;
    }

    // Grouping is a graph-style post-processing step driven by the supplied
    // relatedness cutoff. It currently operates on dense GRM input only.
    if(group){
        spdlog::info("Grouping individuals by GRM relationships");
        require_grm_file();
        if (sparse_input) {
            spdlog::error("--group no longer supports --sparse");
            std::exit(1);
        }
        this->group_related_samples(grm_file, cut_value);
        return 0;
    }

    // `out_close` writes one representative from each overly related pair.
    if (out_close) {
        spdlog::info("Outputting one of each pair of close individuals");
        require_grm_file();
        this->remove_close_id(grm_file, sparse, cut_value, out_file);
        return 0;
    }

    spdlog::error("No valid operation specified");
    std::exit(1);
}
