/*
 * Descripttion: 
 * version: 
 * Author: Chao Ning
 * Date: 2026-03-31 23:11:16
 * LastEditors: Chao Ning
 * LastEditTime: 2026-03-31 23:11:16
 */


#pragma once

#include <Eigen/Dense>

/**
 * @brief Combine p-values with the ACAT/Cauchy aggregation statistic.
 *
 * The returned length-2 vector stores the ACAT statistic and its combined
 * p-value.
 */
Eigen::VectorXd acat_combine_pvalues(const Eigen::VectorXd& p_values,
                                     const Eigen::VectorXd& weights = Eigen::VectorXd());
