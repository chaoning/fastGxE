/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:18:58
 * LastEditTime: 2026-03-31 22:44:33
 * LastEditors: Chao Ning
 */

#pragma once

#include <vector>
#include <string>
#include <Eigen/Dense>

/**
 * @brief Normalize a line by replacing common delimiters, trimming spaces,
 * and removing the trailing comment segment.
 */
void process_line(std::string &line, const std::string &comment_str = "#");


/**
 * @brief Split a string on runs of whitespace and discard empty fields.
 */
std::vector<std::string> split_string(const std::string &str);


/**
 * @brief Join a string vector with a caller-provided separator.
 */
std::string join_string(const std::vector<std::string> &str_vec, const std::string &split_str = " ");

/**
 * @brief Join an Eigen vector into a single string.
 */
std::string join_string(const Eigen::VectorXd &Vec, const std::string &split_str = " ");

/**
 * @brief Parse a string as double and fail on malformed or non-finite input.
 */
double string_to_double(const std::string &str);


/**
 * @brief Return whether a string matches one of the configured missing-value tags.
 */
bool is_nan(const std::string &str, const std::vector<std::string> &nan_vec);
