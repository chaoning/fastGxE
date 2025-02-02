/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-31 21:15:06
 * @LastEditTime: 2025-01-31 21:15:27
 * @LastEditors: Chao Ning
 */
#pragma once

#include <vector>

std::vector<long long> GenerateDiffNumber(long long min, long long max, long long num);
void GenerateDiffPairNumber(long long min, long long max, long long num, std::vector<long long>& res_vec0, std::vector<long long>& res_vec1);
void GenerateDiffPairNumber2(long long min, long long max, long long num, std::vector<long long>& res_vec0, std::vector<long long>& res_vec1);
