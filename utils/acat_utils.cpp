/*
 * Descripttion: 
 * version: 
 * Author: Chao Ning
 * Date: 2026-03-31 23:10:59
 * LastEditors: Chao Ning
 * LastEditTime: 2026-03-31 23:11:00
 */

#include "acat_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

/**
 * @brief Combine p-values with the ACAT/Cauchy aggregation statistic.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `test/someCommonFun_test.cpp`
 */
Eigen::VectorXd acat_combine_pvalues(const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights) {
    Eigen::VectorXd clipped_p_values = p_values;

    // Clip p-values to avoid extreme values.
    double min_p = clipped_p_values.minCoeff();
    double max_p = 0.99999999;
    for (int i = 0; i < clipped_p_values.size(); ++i) {
        clipped_p_values[i] = std::clamp(clipped_p_values[i], min_p, max_p);
    }

    // If weights are not provided, use equal weights.
    Eigen::VectorXd used_weights;
    if (weights.size() == 0) {
        used_weights = Eigen::VectorXd::Ones(clipped_p_values.size());
    } else {
        used_weights = weights;
        assert(used_weights.size() == clipped_p_values.size() && "Weights and p-values must have the same length.");
    }

    // Transform each p-value using the tangent function.
    Eigen::VectorXd transformed_values(clipped_p_values.size());
    for (int i = 0; i < clipped_p_values.size(); ++i) {
        transformed_values[i] = std::tan((0.5 - clipped_p_values[i]) * 3.14159265358979323846);
        if (clipped_p_values[i] < 1.0e-15) {
            transformed_values[i] = 1 / (clipped_p_values[i] * 3.14159265358979323846);
        }
    }

    double T_ACAT = (used_weights.array() * transformed_values.array()).sum();
    double w = used_weights.sum();

    double combined_p_value = 0.5 - (std::atan(T_ACAT / w) / 3.14159265358979323846);
    if (combined_p_value < 1e-50) {
        combined_p_value = w / (T_ACAT * 3.14159265358979323846);
    }

    Eigen::VectorXd result(2);
    result << T_ACAT, combined_p_value;
    return result;
}
