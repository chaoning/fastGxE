/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-30 14:13:33
 * @LastEditTime: 2025-02-01 15:10:42
 * @LastEditors: Chao Ning
 */
#include <cstdint>

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <spdlog/spdlog.h>

#include "iterator_utils.hpp"

using std::vector;
using std::string;
using std::set;
using std::map;

/**
 * @brief Return whether a string vector contains duplicate elements.
 *
 * Used in:
 * - `gmatrix/processGRM.cpp`
 */
bool is_duplicated(std::vector<std::string> vec){
    std::set<std::string> st(vec.begin(), vec.end());
    if(vec.size() != st.size()){
        return true;
    }else{
        return false;
    }
}


/**
 * @brief Remove duplicate strings while preserving the first occurrence order.
 *
 * This helper keeps the original traversal order, unlike the `vector_unique`
 * overloads below that return sorted unique values through `std::set`.
 */
vector<string> vector_unique_keep_order(const vector<string>& vec){
    vector<string> vec2;
    vec2.reserve(vec.size());
    std::unordered_set<string> seen;
    seen.reserve(vec.size());
    for (const auto& value : vec) {
        if (seen.insert(value).second) {
            vec2.push_back(value);
        }
    }
    return vec2;
}

/**
 * @brief Return sorted unique strings.
 *
 * Used in:
 * - `utils/phen.cpp`
 */
vector<string> vector_unique(const vector<string>& vec){
    vector<string> vec2 = vec;
    std::sort(vec2.begin(), vec2.end());
    vec2.erase(std::unique(vec2.begin(), vec2.end()), vec2.end());
    return vec2;
}

/**
 * @brief Return sorted unique integer values.
 */
vector<std::int64_t> vector_unique(const vector<std::int64_t>& vec){
    vector<std::int64_t> vec2 = vec;
    std::sort(vec2.begin(), vec2.end());
    vec2.erase(std::unique(vec2.begin(), vec2.end()), vec2.end());
    return vec2;
}


/**
 * @brief Merge nested integer vectors into one sorted unique vector.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 */
vector<std::int64_t> merge_unique_vectors(const vector<vector<std::int64_t>>& vec2D){
    set<std::int64_t> merged_set;
    for (const auto& vec : vec2D) {
        merged_set.insert(vec.begin(), vec.end());
    }
    vector<std::int64_t> merged_vec;
    merged_vec.assign(merged_set.begin(), merged_set.end());
    return merged_vec;
}

/**
 * @brief Merge nested string vectors into one sorted unique vector.
 */
vector<string> merge_unique_vectors(const vector<vector<string>>& vec2D){
    set<string> merged_set;
    for (const auto& vec : vec2D) {
        merged_set.insert(vec.begin(), vec.end());
    }
    vector<string> merged_vec;
    merged_vec.assign(merged_set.begin(), merged_set.end());
    return merged_vec;
}


/**
 * @brief Find the indices of source strings in a destination string vector.
 *
 * Missing elements are appended to `strNoFound`.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 * - `fastgxe/mmsusie.cpp`
 */
vector<std::int64_t> find_index(const vector<string>& strDestination, const vector<string>& strSource, vector<string>& strNoFound){
    std::unordered_map<string, std::int64_t> strDestination_map;
    strDestination_map.reserve(strDestination.size());
    std::int64_t i = 0;
    for(auto tmp:strDestination){
        strDestination_map[tmp] = i++;
    }

    vector<std::int64_t> index_vec;
    index_vec.reserve(strSource.size());
    for(auto tmp:strSource){
        auto it = strDestination_map.find(tmp);
        if(it == strDestination_map.end()){
            strNoFound.push_back(tmp);
        }else{
            index_vec.push_back(it->second);
        }
    }
    return index_vec;
}


/**
 * @brief Expands ranges specified using ':' notation within a list of variables.
 *
 * This function takes a list of user-specified variables, some of which may be in range notation
 * (e.g., "var1:var2"), and expands them into individual variables based on a predefined ordered list.
 * It also ensures that no duplicate variables appear in the final result.
 *
 * Used in:
 * - `fastgxe/fastgxe.cpp`
 * - `fastgxe/mom.cpp`
 * - `fastgxe/mmsusie.cpp`
 *
 * @param input_vars A vector of strings representing the user-specified variables.
 *                   Can include single variables (e.g., "A") and ranges (e.g., "A:C").
 * @param all_vars A vector of strings representing the full list of valid variables in a predefined order.
 *                 This order determines how ranges are expanded.
 * @return std::vector<std::string> A vector containing the expanded list of variables.
 *
 * @note If an invalid range is encountered (e.g., start variable appears after end variable in all_vars),
 *       or if an unknown variable is provided, or if duplicates are found in the final output, the function logs an error and exits.
 */
std::vector<std::string> expand_variable_ranges(const std::vector<std::string>& input_vars, 
                                                const std::vector<std::string>& all_vars) {
    std::vector<std::string> expanded_vars;
    expanded_vars.reserve(input_vars.size());
    std::unordered_set<std::string> seen_vars;
    seen_vars.reserve(all_vars.size());

    std::unordered_map<std::string, std::size_t> all_var_index;
    all_var_index.reserve(all_vars.size());
    for (std::size_t i = 0; i < all_vars.size(); ++i) {
        all_var_index[all_vars[i]] = i;
    }

    for (const auto& var : input_vars) {
        std::size_t colon_pos = var.find(':');
        if (colon_pos != std::string::npos) {
            std::string start_var = var.substr(0, colon_pos);
            std::string end_var = var.substr(colon_pos + 1);

            auto start_it = all_var_index.find(start_var);
            auto end_it = all_var_index.find(end_var);

            if (start_it != all_var_index.end() &&
                end_it != all_var_index.end() &&
                start_it->second <= end_it->second) {
                for (std::size_t i = start_it->second; i <= end_it->second; ++i) {
                    const auto& expanded = all_vars[i];
                    if (!seen_vars.insert(expanded).second) {
                        spdlog::error("Duplicate variable detected: {}", expanded);
                        exit(1);
                    }
                    expanded_vars.push_back(expanded);
                }
            } else {
                spdlog::error("Invalid range: {}", var);
                exit(1);
            }
        } else {
            if (all_var_index.find(var) != all_var_index.end()) {
                if (!seen_vars.insert(var).second) {
                    spdlog::error("Duplicate variable detected: {}", var);
                    exit(1);
                }
                expanded_vars.push_back(var);
            } else {
                spdlog::error("Invalid variable: {}", var);
                exit(1);
            }
        }
    }
    return expanded_vars;
}
