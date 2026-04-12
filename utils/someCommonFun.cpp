/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-02 18:24:20
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-20 16:48:57
 */

#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <cmath>


using Eigen::VectorXd;

/**
 * @brief Check whether iterative optimization has met any stopping criterion.
 *
 * This helper updates the current parameter-change and gradient norms, and
 * optionally logs the current log-likelihood and convergence diagnostics.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 */
bool convergence_criteria(double logL_curr, double logL, double cc_logL, 
        const VectorXd& delta_Vec, double& cc_par_val, double cc_par, 
        const VectorXd& fd_Vec, double& cc_gra_val, double cc_gra, bool isPrint){
    cc_par_val = delta_Vec.norm();
    cc_gra_val = fd_Vec.norm();

    if(isPrint){
        spdlog::info("logL: {:.5f}", logL_curr);
        spdlog::info("Change in parameters: {}", cc_par_val);
        spdlog::info("Norm of gradient vector: {}", cc_gra_val);
    }

    return std::fabs(logL_curr - logL) < cc_logL ||
           cc_par_val < cc_par ||
           cc_gra_val < cc_gra;
}
