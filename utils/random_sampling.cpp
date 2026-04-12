/*
 * @Description:
 * @Author: Chao Ning
 * @Date: 2025-01-31 21:15:06
 * LastEditTime: 2026-03-31 23:18:26
 * LastEditors: Chao Ning
 */

#include <cstdint>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>

#include "random_sampling.hpp"

using std::cout;
using std::endl;
using std::vector;

/**
 * @brief Sample distinct integers from the half-open interval [`min`, `max`).
 *
 * The returned values are unique and stored in randomized order.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `test/random_sampling_example.cpp`
 */
vector<std::int64_t> sample_unique_integers(std::int64_t min, std::int64_t max, std::int64_t num)
{
    if (num > max - min) {
        cout << num << " exceeds the upper limit" << endl;
        exit(1);
    }

    vector<std::int64_t> sampled_values(static_cast<size_t>(max - min));
    std::iota(sampled_values.begin(), sampled_values.end(), min);

    // Seed once per thread to avoid repeated same-second reseeding.
    static thread_local std::mt19937_64 rng(std::random_device{}());
    std::shuffle(sampled_values.begin(), sampled_values.end(), rng);

    if (num < static_cast<std::int64_t>(sampled_values.size())) {
        sampled_values.resize(static_cast<size_t>(num));
    }

    return sampled_values;
}
