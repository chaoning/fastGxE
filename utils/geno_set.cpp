/*
 * @Description: One or more PLINK filesets treated as a single ordered SNP list
 * @Author: Chao Ning
 */
#define EIGEN_USE_MKL_ALL  // must be before include Eigen

#include "geno_set.hpp"

#include <algorithm>
#include <fstream>
#include <numeric>

#include "fatal_error.hpp"
#include "string_utils.hpp"

GenoSet::GenoSet(const std::vector<std::string>& prefixes) : _prefixes(prefixes) {
    if (_prefixes.empty()) {
        fatal_error("No PLINK fileset was given.");
    }

    _genos.reserve(_prefixes.size());
    _offsets.reserve(_prefixes.size() + 1);
    _offsets.push_back(0);
    for (const std::string& prefix : _prefixes) {
        _genos.push_back(std::make_unique<GENO>(prefix));
        _offsets.push_back(_offsets.back() + _genos.back()->get_num_sid());
    }
    _sample_index.resize(_prefixes.size());

    if (num_snp() <= 0) {
        fatal_error("The PLINK input contains no SNPs.");
    }
    spdlog::info("Genotype input: {} fileset(s), {} SNPs", num_files(), num_snp());
}

std::vector<std::string> GenoSet::read_prefix_list(const std::string& list_file) {
    std::ifstream fin(list_file);
    if (!fin.is_open()) {
        fatal_error("Failed to open the fileset list: {}", list_file);
    }

    std::vector<std::string> prefixes;
    std::string one_line;
    while (std::getline(fin, one_line)) {
        process_line(one_line);
        if (one_line.empty()) continue;
        const std::vector<std::string> fields = split_string(one_line);
        if (fields.size() != 1) {
            fatal_error("Each line of {} must hold exactly one PLINK prefix: {}", list_file, one_line);
        }
        prefixes.push_back(fields[0]);
    }
    if (prefixes.empty()) {
        fatal_error("The fileset list is empty: {}", list_file);
    }
    return prefixes;
}

void GenoSet::bind_samples(const std::vector<std::string>& sample_ids,
                           const std::vector<std::size_t>& file_indices) {
    for (const std::size_t file_index : file_indices) {
        GENO& geno = *_genos[file_index];
        geno.validate_bed_size();
        spdlog::info("Matching samples to {}: {} individuals, {} SNPs",
                     _prefixes[file_index], geno.get_num_iid(), geno.get_num_sid());
        _sample_index[file_index] = geno.find_fam_index(sample_ids);
    }
}

void GenoSet::bind_samples(const std::vector<std::string>& sample_ids) {
    std::vector<std::size_t> all_files(num_files());
    std::iota(all_files.begin(), all_files.end(), std::size_t{0});
    bind_samples(sample_ids, all_files);
}

const std::vector<std::int64_t>& GenoSet::sample_index(std::size_t file_index) const {
    if (_sample_index[file_index].empty()) {
        fatal_error("Samples were not matched to fileset {}.", _prefixes[file_index]);
    }
    return _sample_index[file_index];
}

std::size_t GenoSet::file_of(std::int64_t global_index) const {
    const auto next_offset = std::upper_bound(_offsets.begin(), _offsets.end(), global_index);
    return static_cast<std::size_t>(next_offset - _offsets.begin()) - 1;
}

std::vector<SnpSegment> GenoSet::segments(std::int64_t start, std::int64_t end) const {
    if (start < 0 || end > num_snp() || start >= end) {
        fatal_error("Invalid SNP range [{}, {}) for {} SNPs.", start, end, num_snp());
    }

    std::vector<SnpSegment> result;
    std::int64_t position = start;
    while (position < end) {
        const std::size_t file_index = file_of(position);
        SnpSegment segment;
        segment.file_index = file_index;
        segment.local_start = position - _offsets[file_index];
        segment.global_start = position;
        segment.num_snp = std::min(end, _offsets[file_index + 1]) - position;
        result.push_back(segment);
        position += segment.num_snp;
    }
    return result;
}

void GenoSet::read_snps_by_index(const std::vector<std::int64_t>& global_indices, MatrixXd& snp_mat,
                                 VectorXd& freq_arr, VectorXd& missing_rate_arr, VectorXd& nobs_geno_arr) {
    // The filesets are read one at a time, so remember where each SNP belongs.
    std::vector<std::vector<std::int64_t>> local_indices(num_files());
    std::vector<std::vector<std::int64_t>> positions(num_files());
    for (std::size_t position = 0; position < global_indices.size(); ++position) {
        const std::int64_t global_index = global_indices[position];
        if (global_index < 0 || global_index >= num_snp()) {
            fatal_error("The requested SNP index {} exceeds the number of SNPs {}.", global_index, num_snp());
        }
        const std::size_t file_index = file_of(global_index);
        local_indices[file_index].push_back(global_index - _offsets[file_index]);
        positions[file_index].push_back(static_cast<std::int64_t>(position));
    }

    const std::int64_t num_requested = static_cast<std::int64_t>(global_indices.size());
    snp_mat.resize(0, 0);
    freq_arr.resize(num_requested);
    missing_rate_arr.resize(num_requested);
    nobs_geno_arr.resize(num_requested);

    for (std::size_t file_index = 0; file_index < num_files(); ++file_index) {
        if (local_indices[file_index].empty()) continue;

        MatrixXd part_mat;
        VectorXd part_freq, part_missing_rate, part_nobs;
        _genos[file_index]->read_bed_by_snp_indices(
            _prefixes[file_index] + ".bed", part_mat, part_freq, part_missing_rate, part_nobs,
            local_indices[file_index], sample_index(file_index));

        if (snp_mat.size() == 0) {
            snp_mat.resize(part_mat.rows(), num_requested);
        }
        for (std::size_t k = 0; k < positions[file_index].size(); ++k) {
            const std::int64_t position = positions[file_index][k];
            snp_mat.col(position) = part_mat.col(k);
            freq_arr(position) = part_freq(k);
            missing_rate_arr(position) = part_missing_rate(k);
            nobs_geno_arr(position) = part_nobs(k);
        }
    }
}
