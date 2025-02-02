/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:18:58
 * @LastEditTime: 2025-02-01 21:57:43
 * @LastEditors: Chao Ning
 */
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <Eigen/Dense>


void process_line(std::string &line, const std::string comment_str = "#");


std::vector<std::string> split_string(const std::string &str);



std::string join_string(const std::vector<std::string> &str_vec, const std::string &split_str = " ");

std::string join_string(const Eigen::VectorXd &Vec, const std::string &split_str = " ");


double string_to_double(const std::string &str, char *endptr = nullptr);



long long string_to_longlong(const std::string &str, char *endptr = nullptr);


bool is_nan(const std::string &str, const std::vector<std::string> &nan_vec);


std::string double_to_string(double value, int decimal = 4);
