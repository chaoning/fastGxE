/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 14:18:58
 * LastEditTime: 2026-03-31 22:44:55
 * LastEditors: Chao Ning
 */
 

#include "string_utils.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <charconv>
#include <algorithm>
#include <cmath>
#include "fatal_error.hpp"

using std::string;
using std::vector;
using std::ostringstream;
using std::cout;
using std::endl;

/**
 * @brief Normalize a text line for downstream token parsing.
 *
 * This helper removes the trailing comment segment once `comment_str` is
 * matched, converts common separators to spaces, and trims leading/trailing
 * blanks in place.
 *
 * Used in:
 * - `gmatrix/processGRM.cpp`
 * - `gmatrix/gmatrix.cpp`
 * - `utils/geno.cpp`
 * - `utils/phen.cpp`
 */
void process_line(string &line, const string &comment_str) {
    if (!comment_str.empty()) {
        // Match the full comment token rather than any single character in it.
        size_t comment_pos = line.find(comment_str);
        if (comment_pos != string::npos) {
            line.resize(comment_pos);
        }
    }

    size_t write = 0;
    size_t first_non_space = string::npos;
    size_t last_non_space = string::npos;

    // Normalize delimiters and track the non-space window in a single pass.
    for (size_t read = 0; read < line.size(); ++read) {
        char c = line[read];
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n') {
            c = ' ';
        }
        line[write] = c;
        if (c != ' ') {
            if (first_non_space == string::npos) {
                first_non_space = write;
            }
            last_non_space = write;
        }
        ++write;
    }

    if (first_non_space == string::npos) {
        line.clear();
        return;
    }

    line.erase(last_non_space + 1);
    line.erase(0, first_non_space);
}


/**
 * @brief Split a string on runs of whitespace.
 *
 * Empty fields are skipped, so the returned vector only contains non-empty
 * tokens.
 *
 * Used in:
 * - `gmatrix/processGRM.cpp`
 * - `gmatrix/gmatrix.cpp`
 * - `utils/geno.cpp`
 * - `utils/phen.cpp`
 * - `utils/random_sampling.cpp`
 */
vector<string> split_string(const std::string& s) {
    vector<string> result;
    if (s.empty()) {
        return result;
    }

    size_t start = 0;
    const size_t n = s.size();

    while (start < n) {
        while (start < n && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
        if (start >= n) break;

        size_t end = start;
        while (end < n && !std::isspace(static_cast<unsigned char>(s[end]))) ++end;

        result.emplace_back(s, start, end - start);
        start = end;
    }
    return result;
}


/**
 * @brief Join a string vector with a caller-provided separator.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
string join_string(const std::vector<std::string> &str_vec, const std::string &split_str) {
    if (str_vec.empty()) return "";
    if (str_vec.size() == 1) return str_vec.front();

    size_t total_size = split_str.size() * (str_vec.size() - 1);
    for (const auto &item : str_vec) {
        total_size += item.size();
    }

    string result;
    result.reserve(total_size);
    result += str_vec.front();
    for (size_t i = 1; i < str_vec.size(); ++i) {
        result += split_str;
        result += str_vec[i];
    }
    return result;
}

/**
 * @brief Join an Eigen vector into a single string.
 *
 * Used in:
 * - `fastgxe/mmsusie.cpp`
 */
std::string join_string(const Eigen::VectorXd &Vec, const std::string &split_str){
    if (Vec.size() == 0) return "";

    std::ostringstream oss;
    oss << Vec[0];
    for (Eigen::Index i = 1; i < Vec.size(); ++i) {
        oss << split_str;
        oss << Vec[i];
    }
    return oss.str();
}

/**
 * @brief Parse a string as double and fail fast on malformed input.
 *
 * The conversion is strict: empty input, trailing non-numeric characters, and
 * non-finite results are all treated as errors.
 *
 * Used in:
 * - `gmatrix/gmatrix.cpp`
 * - `utils/geno.cpp`
 * - `utils/phen.cpp`
 */
double string_to_double(const std::string &str) {
    double result = 0.0;
    const char *begin = str.data();
    const char *end = begin + str.size();
    auto [parse_end, ec] = std::from_chars(begin, end, result);

    if (ec != std::errc() || parse_end == begin || parse_end != end) {
        string invalid_suffix;
        if (parse_end < end) {
            invalid_suffix.assign(parse_end, end);
        }
        cout << "ERROR: '" << str << "' contains invalid characters: '" << invalid_suffix
             << "', cannot be converted to double." << endl;
        exit(EXIT_FAILURE);
    }
    if (!std::isfinite(result)) {
        cout << "ERROR: The number is not finite: " << result << endl;
        exit(EXIT_FAILURE);
    }
    return result;
}


/**
 * @brief Return whether `str` matches one of the configured missing-value tags.
 *
 * Used in:
 * - `utils/geno.cpp`
 * - `utils/phen.cpp`
 */
bool is_nan(const std::string &str, const std::vector<std::string> &nan_vec) {
    return std::find(nan_vec.begin(), nan_vec.end(), str) != nan_vec.end();
}
