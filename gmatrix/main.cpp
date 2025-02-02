/*
 * @Description: Genomic Relationship Matrix (GRM) Computation CLI
 * @Author: Chao Ning
 * @Date: 2025-01-30
 * @LastEditTime: 2025-01-30 23:56:28
 * @LastEditors: Chao Ning
 */

#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <spdlog/spdlog.h>

#include "gmatrix.hpp"
#include "processGRM.hpp"

// Display the general help message
void show_help_general() {
    std::cout << "Usage: gmatrix [options]\n"
              << "Options:\n"
              << "  -h 1   Compute GRM\n"
              << "  -h 2   Process GRM\n";
}

// Compute the Genomic Relationship Matrix (GRM)
void compute_grm(int argc, char* argv[]) {
    Gmatrix GRM;
    GRM.run(argc, argv);
}

// Process the Genomic Relationship Matrix (GRM)
void process_grm(int argc, char* argv[]) {
    ProcessGRM ProcessGRMA;
    ProcessGRMA.run(argc, argv);
}

int main(int argc, char* argv[]) {
    // Define command-line options
    struct option long_options[] = {
        {"make-grm", no_argument, nullptr, 'm'},        // Compute GRM
        {"process-grm", no_argument, nullptr, 'p'},     // Process GRM
        {"help", optional_argument, nullptr, 'h'},      // Help (-h [1|2])
        {nullptr, 0, nullptr, 0}
    };

    if (argc == 1) {
        show_help_general();
        return 1;
    }

    int opt;
    while ((opt = getopt_long(argc, argv, "h::mp", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h': {
                // Manually check argv if optarg is nullptr
                if (optarg == nullptr) {
                    if (optind < argc && argv[optind][0] != '-') {
                        optarg = argv[optind];  // Assign the next argument manually
                        optind++;
                    }
                }

                if (optarg == nullptr) {
                    show_help_general();  // No argument, show general help
                } else {
                    int help_type = std::atoi(optarg);
                    if (help_type == 1) {
                        compute_grm(argc, argv);
                    } else if (help_type == 2) {
                        process_grm(argc, argv);
                    } else {
                        spdlog::error("Invalid argument for -h");
                        show_help_general();
                    }
                }
                return 0;
            }
            case 'm':  // Compute GRM
                compute_grm(argc, argv);
                return 0;
            case 'p':  // Process GRM
                process_grm(argc, argv);
                return 0;
            default:
                spdlog::error("Unknown option provided.");
                show_help_general();
                return 1;
        }
    }

    return 1;
}
