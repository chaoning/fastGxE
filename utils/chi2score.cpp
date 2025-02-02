/*
 * @Description: Calculate the P values from linear combination of chi2 distributions
 * @Author: Chao Ning
 * @Date: 2023-05-25 21:05:03
 * @LastEditTime: 2023-06-13 16:05:34
 * @LastEditors: Chao Ning
 */

#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_cdf.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/normal.hpp>



double LiuSF(double t, Eigen::VectorXd &lambs, Eigen::VectorXd &dofs, Eigen::VectorXd &deltas, 
                        bool kurtosis){
    // lambs**1, lambs**2, lambs**3, lambs**4
    std::vector <Eigen::VectorXd> lambs_vec(4);
    for(long long i = 0; i < 4; i++){
        lambs_vec[i] = lambs.array().pow(i+1);
    }

    // c[0..3]
    std::vector <double> c_vec(4);
    for(long long i = 0; i < 4; i++){
        c_vec[i] = (lambs_vec[i].cwiseProduct(dofs)).sum() + (i + 1) * (lambs_vec[i].cwiseProduct(deltas)).sum();
    }

    // s1, s2, s1**2
    double s1 = c_vec[2] / std::pow(c_vec[1], 3.0/2.0);
    double s2 = c_vec[3] / std::pow(c_vec[1], 2.0);
    double s12 = std::pow(s1, 2);

    double delta_x, dof_x, a;
    if(s12 > s2){
        a = 1 / (s1 - std::sqrt(s12 - s2));
        delta_x = s1 * std::pow(a, 3) - std::pow(a, 2);
        dof_x = std::pow(a, 2) - 2 * delta_x;
    }else{
        delta_x = 0;
        if(kurtosis){
            a = 1 / std::sqrt(s2);
            dof_x = 1 / s2;
        }
        else{
            a = 1 / s1;
            dof_x = 1 / s12;
        }
    }

    double mu_q = c_vec[0];
    double sigma_q = std::sqrt(2 * c_vec[1]);

    double mu_x = dof_x + delta_x;
    double sigma_x = std::sqrt(2 * (dof_x + 2 * delta_x));

    double t_star = (t - mu_q) / sigma_q;
    double tfinal = t_star * sigma_x + mu_x;
    
    if(tfinal < 0) return 1.0;

    if(delta_x < 1.e-9) delta_x = 1.0e-9;
    boost::math::non_central_chi_squared_distribution<> dist(dof_x, delta_x);
    double p = boost::math::cdf(complement(dist, tfinal));

    return p;
}



double saddle(double x, Eigen::VectorXd& lambda) {
    double d = lambda.maxCoeff();
    Eigen::VectorXd lambda_norm = lambda / d;
    double x_norm = x / d;

    double lmin, lmax, hatzeta, w, v;
    if (lambda_norm.minCoeff() < 0) {
        std::vector<double> lambda_norm_vec;
        for(auto tmp:lambda_norm){
            if(tmp < 0) lambda_norm_vec.push_back(tmp);
        }
        auto max_element_iter = std::max_element(lambda_norm_vec.begin(), lambda_norm_vec.end());
        lmin = 1 / (2 * *max_element_iter) * 0.99999;
    } else if (x_norm > lambda_norm.sum()) {
        lmin = -0.01;
    } else {
        lmin = -lambda_norm.size() / (2 * x_norm);
    }

    lmax = 1 / (2 * lambda_norm.maxCoeff()) * 0.99999;

    // root-finding function
    int digits = std::numeric_limits<double>::digits;
    boost::uintmax_t max_iter = 500;
    auto f = [&](double zeta) {
        return (lambda_norm.array() / (1 - 2*zeta*lambda_norm.array())).sum() - x_norm;
    };
    std::pair<double, double> result = boost::math::tools::toms748_solve(f, lmin, lmax, boost::math::tools::eps_tolerance<double>(digits), max_iter);

    hatzeta = result.first + (result.second - result.first) / 2;

    double k0 = -(1 - 2 * hatzeta * lambda_norm.array()).log().sum() / 2;
    w = std::copysign(std::sqrt(2*(hatzeta*x_norm - k0)), hatzeta);
    double kpprime0 = 2 * (lambda_norm.array().square() / (1 - 2*hatzeta*lambda_norm.array()).square()).sum();
    v = hatzeta*std::sqrt(kpprime0);

    if (std::abs(hatzeta) < 1e-4) {
        double tr = lambda.mean();
        double tr2 = lambda.array().square().mean() / (tr * tr);
        double scale = tr*tr2;
        double df = lambda.size() / tr2;
        double guess = gsl_cdf_chisq_Q(x / scale, df);
        return guess;
    } else {
        boost::math::normal_distribution<> dist;
        return boost::math::cdf(complement(dist, w + std::log(v / w) / w));
    }
}
