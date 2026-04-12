/*
 * @Descripttion: Header file for Gmatrix class, used for genomic relationship matrix calculations
 * @version: 
 * @Author: Chao Ning
 * @Date: 2021-08-20 09:05:52
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-01-30 14:29:31
 */
#include <cstdint>

#pragma once

#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

using std::vector;
using std::string;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowMajor;


class Gmatrix{
public:
    /**
     * @brief Constructor for the Gmatrix class.
     */
    Gmatrix() = default;

    /**
     * @brief Destructor for Gmatrix class.
     */
    ~Gmatrix() = default;

    int run(int argc, char* argv[]);
    
    
    /**
     * @brief Partition the individuals (IID) for parallel computation.
     * @param npart Total number of partitions.
     * @param ipart Current partition index.
     */
    void iid_part(std::int64_t npart, std::int64_t ipart);

    /**
     * @brief Output individual identifiers.
     */
    void output_iid();

    /**
     * @brief Compute the additive genomic relationship matrix (AGRM).
     * @param npart_snp Number of SNPs per partition.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed genomic relationship matrix.
     */
    MatrixXd agmatrix(std::int64_t npart_snp, double maf, double missing_rate_cut, int code_type);

     /**
     * @brief Compute the dominance genomic relationship matrix (DGMRM) using allele sharing.
     * @param npart Number of partitions.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed dominance genomic relationship matrix.
     */
    MatrixXd dgmatrix_as(std::int64_t npart, double maf, double missing_rate_cut, int code_type);

    /**
     * @brief Compute the dominance genomic relationship matrix (DGMRM) using genotype similarity.
     * @param npart Number of partitions.
     * @param maf Minor allele frequency threshold.
     * @param code_type Encoding type for genotypes.
     * @return Computed dominance genomic relationship matrix.
     */
    MatrixXd dgmatrix_gs(std::int64_t npart, double maf, double missing_rate_cut, int code_type);

    /**
     * @brief Set the output file path.
     * @param file Path to the output file.
     */
    void set_out_file(const std::string& file);

     /**
     * @brief Output the computed genomic relationship matrix.
     * @param gmat Computed genomic relationship matrix.
     * @param val Replacement value for filtered entries.
     */
    void out_gmat(const MatrixXd& gmat, const double val);

    /**
     * @brief Output the number of SNPs used in computation.
     */
    void out_num_snp_used();

private:
    string m_geno_file;           ///< Path to the genotype file.
    string m_out_file;            ///< Path to the output file.
    vector<string> m_id_in_geno_vec;  ///< List of IDs in the genotype file.
    vector<string> m_missing_in_geno_vec; ///< List of missing genotype indicators.

    std::int64_t m_num_snp = 0;              ///< Number of SNPs.
    std::int64_t m_num_id = 0;               ///< Number of individuals (IDs).
    int m_input_geno_fmt = 0;         ///< Format of the input genotype data.
    std::int64_t m_num_snp_used = 0; ///< Number of SNPs used in computation.
    double m_scale = 0.0;        ///< Scaling factor for genomic relationship matrix.
    std::int64_t m_start_id = 0;  ///< Start index for processing individuals.
    std::int64_t m_num_id_read = 0; ///< Number of individuals read.

};
