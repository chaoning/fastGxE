/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 20:44:09
 * @LastEditTime: 2025-01-30 15:00:17
 * @LastEditors: Chao Ning
 */

#include <cstdint>

#include <stdexcept>

#include <spdlog/spdlog.h>
#include "partitionLowerTriangle.hpp"

namespace {

using int128 = __int128_t;

int128 triangularElements(std::int64_t rows) {
    return static_cast<int128>(rows) * (rows + 1) / 2;
}

int128 abs128(int128 value) {
    return value < 0 ? -value : value;
}

std::int64_t chooseBoundary(std::int64_t n, std::int64_t npart, std::int64_t boundary_index,
                         std::int64_t previous_boundary, int128 total_elements) {
    const std::int64_t min_boundary = previous_boundary + 1;
    const std::int64_t max_boundary = n - (npart - boundary_index);

    if (min_boundary > max_boundary) {
        throw std::runtime_error("Unable to assign at least one row to each partition.");
    }

    const int128 scaled_target = total_elements * boundary_index;

    std::int64_t left = min_boundary;
    std::int64_t right = max_boundary;
    while (left < right) {
        const std::int64_t mid = left + (right - left) / 2;
        if (triangularElements(mid) * npart < scaled_target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    std::int64_t upper = left;
    if (triangularElements(upper) * npart < scaled_target && upper < max_boundary) {
        ++upper;
    }

    std::int64_t lower = upper;
    if (lower > min_boundary) {
        --lower;
    }

    const int128 upper_distance = abs128(triangularElements(upper) * npart - scaled_target);
    const int128 lower_distance = abs128(triangularElements(lower) * npart - scaled_target);

    // Prefer the earlier boundary on ties so earlier partitions do not overrun.
    return lower_distance <= upper_distance ? lower : upper;
}

}  // namespace

/**
 * @brief Finds the partition bounds for a given part in a lower triangular matrix.
 *
 * @param n Total number of rows in the matrix.
 * @param npart Total number of partitions.
 * @param ipart The index of the partition (1-based).
 * @param start_row Output reference for the starting row index of the partition.
 * @param row_count Output reference for the number of rows in the partition.
 */
void partitionLowerTriangle(std::int64_t n, std::int64_t npart, std::int64_t ipart,
                            std::int64_t& start_row, std::int64_t& row_count) {
    if (n <= 0) {
        spdlog::error("n must be positive, got {}", n);
        throw std::invalid_argument("partitionLowerTriangle: n must be positive.");
    }

    if (npart <= 0) {
        spdlog::error("npart must be positive, got {}", npart);
        throw std::invalid_argument("partitionLowerTriangle: npart must be positive.");
    }

    if (ipart < 1 || ipart > npart) {
        spdlog::error("ipart must be within [1, {}], got {}", npart, ipart);
        throw std::invalid_argument("partitionLowerTriangle: ipart is out of range.");
    }

    if (npart > n) {
        spdlog::error("npart ({}) cannot exceed n ({}) when each partition must contain at least one row.", npart, n);
        throw std::invalid_argument("partitionLowerTriangle: npart cannot exceed n.");
    }

    const int128 total_elements = triangularElements(n);

    std::int64_t previous_boundary = 0;
    std::int64_t end_row = 0;
    start_row = 0;

    for (std::int64_t part = 1; part <= ipart; ++part) {
        std::int64_t current_boundary = n;
        if (part < npart) {
            current_boundary = chooseBoundary(n, npart, part, previous_boundary, total_elements);
        }

        if (part == ipart) {
            start_row = previous_boundary;
            end_row = current_boundary;
        }

        previous_boundary = current_boundary;
    }

    row_count = end_row - start_row;
    if (row_count <= 0) {
        spdlog::error("Partition {}/{} resolved to an empty row range [{}, {}).", ipart, npart, start_row, end_row);
        throw std::runtime_error("partitionLowerTriangle: partition has zero rows.");
    }

}
