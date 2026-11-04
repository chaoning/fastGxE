/*
 * @Description: One or more PLINK filesets treated as a single ordered SNP list
 * @Author: Chao Ning
 *
 * Every SNP has a "global" index: its position in the concatenation of the
 * filesets, in the order they were given. A single whole-genome fileset and the
 * same genome split by chromosome therefore share the same global indices, so a
 * random draw (step 2) or an even split of the work (--split-task, step 3) does
 * not depend on how the genotypes are stored on disk.
 */
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "geno.hpp"

// A run of consecutive SNPs inside one fileset.
struct SnpSegment {
    std::size_t file_index = 0;
    std::int64_t local_start = 0;   // first SNP, counted within the fileset
    std::int64_t global_start = 0;  // first SNP, counted across all filesets
    std::int64_t num_snp = 0;
};

class GenoSet {
public:
    explicit GenoSet(const std::vector<std::string>& prefixes);

    // Reads a list file with one PLINK prefix per line (blank lines are skipped).
    static std::vector<std::string> read_prefix_list(const std::string& list_file);

    std::size_t num_files() const { return _genos.size(); }
    std::int64_t num_snp() const { return _offsets.back(); }
    const std::string& prefix(std::size_t file_index) const { return _prefixes[file_index]; }
    GENO& geno(std::size_t file_index) { return *_genos[file_index]; }

    // Maps the analysis samples to .fam rows and validates the .bed size for the
    // given filesets. Each fileset is matched by IID on its own, so .fam row order
    // may differ between filesets. Must run before genotypes are read from a file.
    void bind_samples(const std::vector<std::string>& sample_ids,
                      const std::vector<std::size_t>& file_indices);
    void bind_samples(const std::vector<std::string>& sample_ids);

    // .fam rows of the analysis samples for a fileset bound by bind_samples.
    const std::vector<std::int64_t>& sample_index(std::size_t file_index) const;

    // Splits the global range [start, end) at fileset boundaries.
    std::vector<SnpSegment> segments(std::int64_t start, std::int64_t end) const;

    // Reads SNPs by global index. Columns follow `global_indices`, whichever
    // fileset each SNP comes from. Missing genotypes are mean-imputed.
    void read_snps_by_index(const std::vector<std::int64_t>& global_indices, MatrixXd& snp_mat,
                            VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr);

private:
    std::size_t file_of(std::int64_t global_index) const;

    std::vector<std::string> _prefixes;
    std::vector<std::unique_ptr<GENO>> _genos;
    std::vector<std::int64_t> _offsets;  // cumulative SNP counts; size = num_files + 1
    std::vector<std::vector<std::int64_t>> _sample_index;
};
