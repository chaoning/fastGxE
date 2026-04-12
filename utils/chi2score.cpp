/*
 * @Description: Calculate the P values from linear combination of chi2 distributions
 * @Author: Chao Ning
 * @Date: 2023-05-25 21:05:03
 * @LastEditTime: 2023-06-13 16:05:34
 * @LastEditors: Chao Ning
 */

#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <gsl/gsl_cdf.h>
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/normal.hpp>

namespace {

double clamp_probability(double p) {
    if (!std::isfinite(p)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::min(1.0, std::max(0.0, p));
}

double weighted_chisq_sf_fallback(double x, const Eigen::VectorXd& lambda) {
    if (!std::isfinite(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x <= 0.0) {
        return 1.0;
    }

    Eigen::ArrayXd lambda_pos = lambda.array().max(0.0);
    double sum1 = lambda_pos.sum();
    double sum2 = lambda_pos.square().sum();
    if (sum1 <= 0.0 || sum2 <= 0.0) {
        return 0.0;
    }

    double scale = sum2 / sum1;
    double df = sum1 * sum1 / sum2;
    if (!(scale > 0.0) || !(df > 0.0) || !std::isfinite(scale) || !std::isfinite(df)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return clamp_probability(gsl_cdf_chisq_Q(x / scale, df));
}

Eigen::VectorXd sanitize_lambda(const Eigen::VectorXd& lambda) {
    if (lambda.size() == 0) {
        return lambda;
    }

    double max_abs = lambda.cwiseAbs().maxCoeff();
    double zero_tol = std::max(1.0, max_abs) * 1e-12;
    Eigen::VectorXd sanitized(lambda.size());
    for (Eigen::Index i = 0; i < lambda.size(); ++i) {
        double value = lambda(i);
        if (!std::isfinite(value) || value <= zero_tol) {
            sanitized(i) = 0.0;
        } else {
            sanitized(i) = value;
        }
    }
    return sanitized;
}

}  // namespace

double saddle(double x, const Eigen::VectorXd& lambda) {
    if (!std::isfinite(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x <= 0.0) {
        return 1.0;
    }
    if (lambda.size() == 0) {
        return 0.0;
    }

    Eigen::VectorXd lambda_clean = sanitize_lambda(lambda);
    if (lambda_clean.maxCoeff() <= 0.0) {
        return 0.0;
    }

    double d = lambda_clean.maxCoeff();
    Eigen::VectorXd lambda_norm = lambda_clean / d;
    double x_norm = x / d;

    double lmin;
    if (x_norm > lambda_norm.sum()) {
        lmin = -0.01;
    } else {
        lmin = -static_cast<double>(lambda_norm.size()) / (2.0 * x_norm);
    }

    double lmax = 0.5 / lambda_norm.maxCoeff() * 0.99999;
    if (!std::isfinite(lmin) || !std::isfinite(lmax) || !(lmin < lmax)) {
        return weighted_chisq_sf_fallback(x, lambda_clean);
    }

    int digits = std::numeric_limits<double>::digits;
    boost::uintmax_t max_iter = 500;
    auto f = [&](double zeta) {
        return (lambda_norm.array() / (1 - 2 * zeta * lambda_norm.array())).sum() - x_norm;
    };

    double f_lmin = f(lmin);
    double f_lmax = f(lmax);
    if (!std::isfinite(f_lmin) || !std::isfinite(f_lmax) || f_lmin * f_lmax > 0.0) {
        return weighted_chisq_sf_fallback(x, lambda_clean);
    }

    double hatzeta;
    try {
        std::pair<double, double> result = boost::math::tools::toms748_solve(
            f, lmin, lmax, boost::math::tools::eps_tolerance<double>(digits), max_iter);
        hatzeta = result.first + (result.second - result.first) / 2.0;
    } catch (...) {
        return weighted_chisq_sf_fallback(x, lambda_clean);
    }

    double k0 = -(1 - 2 * hatzeta * lambda_norm.array()).log().sum() / 2.0;
    double w_sq = 2.0 * (hatzeta * x_norm - k0);
    double kpprime0 = 2.0 * (lambda_norm.array().square() / (1 - 2 * hatzeta * lambda_norm.array()).square()).sum();
    if (!(w_sq > 0.0) || !(kpprime0 > 0.0) || !std::isfinite(w_sq) || !std::isfinite(kpprime0)) {
        return weighted_chisq_sf_fallback(x, lambda_clean);
    }

    double w = std::copysign(std::sqrt(w_sq), hatzeta);
    double v = hatzeta * std::sqrt(kpprime0);

    if (std::abs(hatzeta) < 1e-4 || std::abs(w) < 1e-8 || !(v / w > 0.0) || !std::isfinite(v)) {
        return weighted_chisq_sf_fallback(x, lambda_clean);
    }

    boost::math::normal_distribution<> dist;
    return clamp_probability(boost::math::cdf(complement(dist, w + std::log(v / w) / w)));
}
