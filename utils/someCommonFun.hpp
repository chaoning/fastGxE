#pragma once


#include<Eigen/Dense>


bool convergence_criteria(double logL_curr, double logL, double cc_logL, 
        VectorXd delta_Vec, double& cc_par_val, double cc_par, 
        VectorXd fd_Vec, double& cc_gra_val, double cc_gra, bool isPrint=true);

double clamp(double val, double min_val, double max_val);

Eigen::VectorXd combine_pvalues_cauchy(const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights = Eigen::VectorXd());

Eigen::VectorXd calWeightZscore(Eigen::VectorXd& beta_Vec, Eigen::VectorXd& se_Vec, Eigen::MatrixXd& EMat);

vector<string> multiple_options_vec(char* _optarg, int _argc, int _optind, char* _argv[]);
