#pragma once

#include <Eigen/Dense>


bool convergence_criteria(double logL_curr, double logL, double cc_logL,
        const Eigen::VectorXd& delta_Vec, double& cc_par_val, double cc_par,
        const Eigen::VectorXd& fd_Vec, double& cc_gra_val, double cc_gra, bool isPrint = true);
