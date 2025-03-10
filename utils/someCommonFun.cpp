/*
 * @Descripttion: 
 * @version: 
 * @Author: Chao Ning
 * @Date: 2022-12-02 18:24:20
 * @LastEditors: Chao Ning
 * @LastEditTime: 2025-02-20 16:48:57
 */

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <string_utils.hpp>
#include <cmath>


using Eigen::VectorXd;

bool convergence_criteria(double logL_curr, double logL, double cc_logL, 
        VectorXd delta_Vec, double& cc_par_val, double cc_par, 
        VectorXd fd_Vec, double& cc_gra_val, double cc_gra, bool isPrint){
    
    cc_par_val = sqrt(delta_Vec.array().square().sum());
    cc_gra_val = sqrt(fd_Vec.array().square().sum());
    
    if(isPrint){
        std::ostringstream buf;
        buf.precision(5);
        buf.setf(std::ios::fixed);
        buf << logL_curr;
        spdlog::info("logL: {}", buf.str());
        spdlog::info("Change in parameters: {}", cc_par_val);
        spdlog::info("Norm of gradient vector: ", cc_gra_val);
    }

    if(std::fabs(logL_curr - logL) < cc_logL || cc_par_val < cc_par || cc_gra_val < cc_gra){
        return true;
    }else{
        return false;
    }
}


// Custom clamp function
double clamp(double val, double min_val, double max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

// Function to combine p-values using the cauchy method
Eigen::VectorXd combine_pvalues_cauchy(const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights) {
    Eigen::VectorXd clipped_p_values = p_values;

    // Clip p-values to avoid extreme values
    double min_p = clipped_p_values.minCoeff();
    double max_p = 0.99999999;
    for (int i = 0; i < clipped_p_values.size(); ++i) {
        clipped_p_values[i] = clamp(clipped_p_values[i], min_p, max_p);
    }

    // If weights are not provided, use equal weights
    Eigen::VectorXd used_weights;
    if (weights.size() == 0) {
        used_weights = Eigen::VectorXd::Ones(clipped_p_values.size());
    } else {
        used_weights = weights;
        assert(used_weights.size() == clipped_p_values.size() && "Weights and p-values must have the same length.");
    }

    // Transform each p-value using the tangent function
    Eigen::VectorXd transformed_values(clipped_p_values.size());
    for (int i = 0; i < clipped_p_values.size(); ++i) {
        transformed_values[i] = std::tan((0.5 - clipped_p_values[i]) * 3.14159265358979323846);
        if(clipped_p_values[i] < 1.0e-15){
            transformed_values[i] = 1 / (clipped_p_values[i] * 3.14159265358979323846);
        }
    }

    // Calculate the weighted sum of the transformed values
    double T_ACAT = (used_weights.array() * transformed_values.array()).sum();

    // Scale parameter w
    double w = used_weights.sum();

    // Calculate the combined p-value from the ACAT statistic
    double combined_p_value = 0.5 - (std::atan(T_ACAT / w) / 3.14159265358979323846);
    // combined_p_value = clamp(combined_p_value, min_p, max_p);

    if(combined_p_value < 1e-50){
        combined_p_value = w / (T_ACAT * 3.14159265358979323846);
    }

    // Create a vector to store the results
    Eigen::VectorXd result(2);
    result << T_ACAT, combined_p_value;

    return result;
}


Eigen::VectorXd calWeightZscore(Eigen::VectorXd& beta_Vec, Eigen::VectorXd& se_Vec, Eigen::MatrixXd& EMat){
    Eigen::VectorXd weight_Vec = beta_Vec.array() / se_Vec.array();
    weight_Vec = weight_Vec.cwiseProduct(weight_Vec);
    Eigen::VectorXd direct_Vec = beta_Vec.array() / beta_Vec.array().abs();
    weight_Vec = (weight_Vec).cwiseProduct(direct_Vec);
    weight_Vec = weight_Vec.array() + 30;
    double sum_weight = weight_Vec.array().abs().sum();
    weight_Vec = weight_Vec.array() / sum_weight;
    Eigen::VectorXd envi_score_Vec = EMat * weight_Vec;
    return envi_score_Vec;
}


std::vector<std::string> multiple_options_vec(char* _optarg, int _argc, int _optind, char* _argv[]){
    std::vector <std::string> seq;
    seq.push_back(_optarg);
    while (_optind < _argc && _argv[_optind][0] != '-') {
        seq.push_back(_argv[_optind]);
        _optind++;
    }
    return seq;
}

