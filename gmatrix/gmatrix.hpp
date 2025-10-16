/*
 * @Descripttion: Header file for Gmatrix class, used for genomic relationship matrix calculations
 * @version: 
 * @Author: Chao Ning
 * @Date: 2021-08-20 09:05:52
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-30 14:29:31
 */
#pragma once

#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using std::vector;
using std::string;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::RowMajor;


class Gmatrix{
public:
    /**
     * @brief Constructor for the Gmatrix class.
     */
    Gmatrix();

    /**
     * @brief Destructor for Gmatrix class.
     */
    ~Gmatrix();

    int run(int argc, char* argv[]);
    
    
    /**
     * @brief Partition the individuals (IID) for parallel computation.
     * @param npart Total number of partitions.
     * @param ipart Current partition index.
     */
    void iid_part(long long npart, long long ipart);

    /**
     * @brief Output individual identifiers.
     */
    void output_iid();

    /**
     * @brief Normalize weight values based on an input weight file.
     * @param weight_arr Vector to store normalized weights.
     * @param weight_file Path to the weight file.
     */
    void normalized_weights(VectorXd& weight_arr, const std::string& weight_file);

    /**
     * @brief Compute the additive genomic relationship matrix (AGRM).
     * @param npart_snp Number of SNPs per partition.
     * @param weight_arr Vector of SNP weights.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed genomic relationship matrix.
     */
    MatrixXd agmatrix(long long npart_snp, VectorXd& weight_arr, double maf, int code_type);

     /**
     * @brief Compute the dominance genomic relationship matrix (DGMRM) using allele sharing.
     * @param npart Number of partitions.
     * @param weight_arr Vector of SNP weights.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed dominance genomic relationship matrix.
     */
    MatrixXd dgmatrix_as(long long npart, VectorXd& weight_arr, double maf, int code_type);

    /**
     * @brief Compute the dominance genomic relationship matrix (DGMRM) using genotype similarity.
     * @param npart Number of partitions.
     * @param weight_arr Vector of SNP weights.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed dominance genomic relationship matrix.
     */
    MatrixXd dgmatrix_gs(long long npart, VectorXd& weight_arr, double maf, int code_type);

    /**
     * @brief Set the output file path.
     * @param file Path to the output file.
     */
    void set_out_file(const std::string& file);

     /**
     * @brief Output the computed genomic relationship matrix.
     * @param gmat Computed genomic relationship matrix.
     * @param sparse Whether to output in sparse format (1 = sparse, 0 = dense).
     * @param cut_val Threshold for filtering small values.
     * @param val Replacement value for filtered entries.
     */
    void out_gmat(const MatrixXd& gmat, const int sparse, const double cut_val, const double val);

    /**
     * @brief Output the number of SNPs used in computation.
     */
    void out_num_snp_used();

private:
    string m_geno_file;           ///< Path to the genotype file.
    string m_out_file;            ///< Path to the output file.
    vector<string> m_id_in_geno_vec;  ///< List of IDs in the genotype file.
    vector<string> m_missing_in_geno_vec; ///< List of missing genotype indicators.

    long long m_num_snp;              ///< Number of SNPs.
    long long m_num_id;               ///< Number of individuals (IDs).
    int m_input_geno_fmt;         ///< Format of the input genotype data.
    long long m_num_snp_used; ///< Number of SNPs used in computation.
    double m_scale;        ///< Scaling factor for genomic relationship matrix.
    long long m_start_id;  ///< Start index for processing individuals.
    long long m_num_id_read; ///< Number of individuals read.

};
