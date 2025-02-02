/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-30 14:14:47
 * @LastEditTime: 2025-02-01 14:11:07
 * @LastEditors: Chao Ning
 */
#pragma once

#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <map>


template <typename T>
T set_intersection_(T v1, T v2){
    T v_intersection;
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::insert_iterator<T>(v_intersection, v_intersection.begin()));
    return v_intersection;
}


template <typename T>
T set_union_(T v1, T v2){
    T v_union;
    std::set_union(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::insert_iterator<T>(v_union, v_union.begin()));
    return v_union;
}


template <typename T>
T set_difference_(T v1, T v2){
    T v_diff;
    std::set_difference(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::insert_iterator<T>(v_diff, v_diff.begin()));
    return v_diff;
}



bool is_duplicated(std::vector<std::string> vec);

std::vector<std::string> vector_unique_keep_order(std::vector<std::string> vec);

std::vector<std::string> vector_unique(std::vector<std::string> vec);

std::vector<long long> vector_unique(std::vector<long long> vec);

std::vector<long long> vector_merged(std::vector<std::vector<long long>> vec2D);

std::vector<std::string> vector_merged(std::vector<std::vector<std::string>> vec2D);

/**
 * @brief Finds the indices of elements in strSource within strDestination. 
 *        Elements that are not found are stored in strNoFound.
 * 
 * @param strDestination Target vector containing reference elements.
 * @param strSource Source vector containing elements to locate.
 * @param strNoFound Reference vector to store elements that could not be found.
 * @return vector<long long> Vector of indices corresponding to strSource elements in strDestination.
 */
std::vector<long long> find_index(const std::vector<std::string>& strDestination, const std::vector<std::string>& strSource, std::vector<std::string>& strNoFound);


/**
 * @brief Expands ranges specified using ':' notation within a list of variables.
 *
 * This function takes a list of user-specified variables, some of which may be in range notation
 * (e.g., "var1:var2"), and expands them into individual variables based on a predefined ordered list.
 * It also ensures that no duplicate variables appear in the final result.
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
                                                const std::vector<std::string>& all_vars);
