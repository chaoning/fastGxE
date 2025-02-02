/*
 * @Description: Genomic Relationship Matrix (GRM) Computation CLI
 * @Author: Chao Ning
 * @Date: 2025-01-30
 * @LastEditTime: 2025-02-01 22:08:13
 * @LastEditors: Chao Ning
 */

#define EIGEN_USE_MKL_ALL  // Must be defined before including Eigen

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <spdlog/spdlog.h>

#include "../gmatrix/gmatrix.hpp"
#include "../gmatrix/processGRM.hpp"
#include "fastgxe.hpp"
#include "mom.hpp"
#include "mmsusie.hpp"

// Display the general help message
void show_help_general() {
    std::cout << "Usage: gmatrix [options]\n"
              << "Options:\n"
              << "  -h 1   Compute GRM\n"
              << "  -h 2   Process GRM\n"
              << "  -h 3   test the SNP main effects\n"
              << "  -h 4   test GxE\n"
              << "  -h 5   GxE heritability estimation with Method-of-moment\n"
              << "  -h 6   mmSuSiE\n";
              
}

// Compute the Genomic Relationship Matrix (GRM)
void compute_grm_(int argc, char* argv[]) {
    Gmatrix GRM;
    GRM.run(argc, argv);
}

// Process the Genomic Relationship Matrix (GRM)
void process_grm_(int argc, char* argv[]) {
    ProcessGRM ProcessGRMA;
    ProcessGRMA.run(argc, argv);
}

void test_main_(int argc, char* argv[]) {
    fastGxE fastA;
    fastA.run(argc, argv);
}

void test_gxe_(int argc, char* argv[]) {
    fastGxE fastA;
    fastA.run(argc, argv);
}

void mom_(int argc, char* argv[]) {
    MoM momA;
    momA.run(argc, argv);
}

void mmsusie_(int argc, char* argv[]) {
    MMSUSIE mmsusieA;
    mmsusieA.run(argc, argv);
}



int main(int argc, char* argv[]) {
    // Define command-line options
    struct option long_options[] = {
        {"make-grm", no_argument, nullptr, 'm'},        // Compute GRM
        {"process-grm", no_argument, nullptr, 'd'},     // Process GRM
        {"test-main", no_argument, nullptr, 'g'},     // test the SNP main effects
        {"test-gxe", no_argument, nullptr, 'x'},     // test GxE
        {"mom", no_argument, nullptr, 'o'},     // GxE heritability estimation with Method-of-moment
        {"mmsusie", no_argument, nullptr, 's'},     // GxE heritability estimation with Method-of-moment
        {"help", optional_argument, nullptr, 'h'}, 
        {nullptr, 0, nullptr, 0}
    };

    if (argc == 1) {
        show_help_general();
        return 1;
    }

    int opt;
    while ((opt = getopt_long(argc, argv, "h::mdgxos", long_options, nullptr)) != -1) {
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
                        compute_grm_(argc, argv);
                    } else if (help_type == 2) {
                        process_grm_(argc, argv);
                    } else if (help_type == 3)
                    {
                        test_main_(argc, argv);
                    } else if (help_type == 4)
                    {
                        test_gxe_(argc, argv);
                    }else if (help_type == 5)
                    {
                        mom_(argc, argv);
                    }else if (help_type == 6)
                    {
                        mmsusie_(argc, argv);
                    }else {
                        spdlog::error("Invalid argument for -h");
                        show_help_general();
                    }
                }
                return 0;
            }
            case 'm':  // Compute GRM
                compute_grm_(argc, argv);
                return 0;
            case 'd':  // Process GRM
                process_grm_(argc, argv);
                return 0;
            case 'g':  // test the SNP main effects
                test_main_(argc, argv);
                return 0;
            case 'x':  // test GxE
                test_gxe_(argc, argv);
                return 0;
            case 'o':
                mom_(argc, argv);
                return 0;
            case 's':
                mmsusie_(argc, argv);
                return 0;
            default:
                spdlog::error("Unknown option provided.");
                show_help_general();
                return 1;
        }
    }

    return 0;
}
