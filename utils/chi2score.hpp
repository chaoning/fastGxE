/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2023-05-26 12:23:40
 * @LastEditTime: 2025-02-01 20:30:59
 * @LastEditors: Chao Ning
 */
#pragma once

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Dense>

double LiuSF(double t, Eigen::VectorXd &lambs, Eigen::VectorXd &dofs, Eigen::VectorXd &deltas, 
                        bool kurtosis = false);

double saddle(double x, Eigen::VectorXd& lambda);
