/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-30 14:13:33
 * @LastEditTime: 2025-02-01 15:10:42
 * @LastEditors: Chao Ning
 */
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <spdlog/spdlog.h>

#include "iterator_utils.hpp"

using std::vector;
using std::string;
using std::set;
using std::map;

bool is_duplicated(std::vector<std::string> vec){
    std::set<std::string> st(vec.begin(), vec.end());
    if(vec.size() != st.size()){
        return true;
    }else{
        return false;
    }
}



vector<string> vector_unique_keep_order(vector<string> vec){
    vector<string> vec2;
    for (vector<string>::iterator it = vec.begin(); it != vec.end(); it++){
        if(find(vec2.begin(), vec2.end(), *it) == vec2.end())
            vec2.push_back(*it);
    }
    return vec2;
}

vector<string> vector_unique(vector<string> vec){
    set<string> set2(vec.begin(), vec.end());
    vector<string> vec2;
    vec2.assign(set2.begin(), set2.end());
    return vec2;
}

vector<long long> vector_unique(vector<long long> vec){
    set<long long> set2(vec.begin(), vec.end());
    vector<long long> vec2;
    vec2.assign(set2.begin(), set2.end());
    return vec2;
}


vector<long long> vector_merged(vector<vector<long long>> vec2D){
    set<long long> merged_set;
    for(vector<vector<long long>>::iterator it = vec2D.begin(); it != vec2D.end(); it++){
        for(vector<long long>::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++){
            merged_set.insert(*it2);
        }
    }
    vector<long long> merged_vec;
    merged_vec.assign(merged_set.begin(), merged_set.end());
    return merged_vec;
}

vector<string> vector_merged(vector<vector<string>> vec2D){
    set<string> merged_set;
    for(vector<vector<string>>::iterator it = vec2D.begin(); it != vec2D.end(); it++){
        for(vector<string>::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++){
            merged_set.insert(*it2);
        }
    }
    vector<string> merged_vec;
    merged_vec.assign(merged_set.begin(), merged_set.end());
    return merged_vec;
}


vector<long long> find_index(const vector<string>& strDestination, const vector<string>& strSource, vector<string>& strNoFound){
    std::map<string, long long> strDestination_map;
    long long i = 0;
    for(auto tmp:strDestination){
        strDestination_map[tmp] = i++;
    }

    vector<long long> index_vec;
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
    std::unordered_set<std::string> seen_vars; // Set to track duplicates

    for (const auto& var : input_vars) {
        size_t colon_pos = var.find(':');
        if (colon_pos != std::string::npos) {
            // Handle range format, e.g., "A:C"
            std::string start_var = var.substr(0, colon_pos);
            std::string end_var = var.substr(colon_pos + 1);

            // Find positions of start and end variables in the predefined list
            auto start_it = std::find(all_vars.begin(), all_vars.end(), start_var);
            auto end_it = std::find(all_vars.begin(), all_vars.end(), end_var);

            if (start_it != all_vars.end() && end_it != all_vars.end() && start_it <= end_it) {
                // Expand the range by adding all variables in between (inclusive)
                for (auto it = start_it; it <= end_it; ++it) {
                    if (seen_vars.count(*it)) {
                        spdlog::error("Duplicate variable detected: {}", *it);
                        exit(1);
                    }
                    seen_vars.insert(*it);
                    expanded_vars.push_back(*it);
                }
            } else {
                spdlog::error("Invalid range: {}", var);
                exit(1); // Exit with an error
            }
        } else {
            // Directly add individual variables if valid
            if (std::find(all_vars.begin(), all_vars.end(), var) != all_vars.end()) {
                if (seen_vars.count(var)) {
                    spdlog::error("Duplicate variable detected: {}", var);
                    exit(1);
                }
                seen_vars.insert(var);
                expanded_vars.push_back(var);
            } else {
                spdlog::error("Invalid variable: {}", var);
                exit(1);
            }
        }
    }
    return expanded_vars;
}
