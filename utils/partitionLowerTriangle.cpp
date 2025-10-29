/*
 * @Description: 
 * @Author: Chao Ning
 * @Date: 2025-01-29 20:44:09
 * @LastEditTime: 2025-01-30 15:00:17
 * @LastEditors: Chao Ning
 */

#include <spdlog/spdlog.h>
#include "partitionLowerTriangle.hpp"

/**
 * @brief Finds the partition bounds for a given part in a lower triangular matrix.
 *
 * @param n Total number of rows in the matrix.
 * @param npart Total number of partitions.
 * @param ipart The index of the partition (1-based).
 * @param start_row Output reference for the starting row index of the partition.
 * @param row_count Output reference for the number of rows in the partition.
 */
void partitionLowerTriangle(long long n, long long npart, long long ipart, 
                            long long& start_row, long long& row_count) {

     // Calculate elements per partition
    long long num_elements_each_part = (n * (n + 1)) / (2 * npart);
    long long start_elements = num_elements_each_part * (ipart - 1);
    long long end_elements = num_elements_each_part * ipart;

    // Compute start_row
    long long num_elements_tmp = 0;
    start_row = 0;
    for(long long i = 0; i < n; i++){
        if(num_elements_tmp >= start_elements){
            start_row = i;
            break;
        }
        num_elements_tmp += i + 1;
    }

    // Compute end_row
    num_elements_tmp = 0;
    long long end_row = 0;
    for(long long i = 0; i < n; i++){
        if(num_elements_tmp >= end_elements){
            end_row = i;
            break;
        }
        num_elements_tmp += i + 1;
    }

    // Ensure last partition captures all remaining rows
    if(ipart == npart) end_row = n;

    // Compute row_count
    row_count = end_row - start_row;
    if(row_count <= 0){
        spdlog::error("Partition {} has zero rows. Reduce --npart argument.", ipart);
        exit(1);
    }

    // Log partition details
    spdlog::info("Partition {}/{} -> Start row: {}, Row count: {}", ipart, npart, start_row, row_count);
    spdlog::info("GRM elements in partition: {}", (start_row + start_row + row_count + 1) * row_count / 2);
}
