/*
 * @Description: Top-level CLI dispatcher for the fastGxE toolkit
 * @Author: Chao Ning
 * @Date: 2025-01-30
 * LastEditors: Chao Ning
 * LastEditTime: 2026-04-12 22:40:00
 */

#include <cstdlib>
#include <iostream>
#include <string_view>

#include <spdlog/spdlog.h>

#include "fastgxe.hpp"
#include "gmatrix.hpp"
#include "mom.hpp"
#include "processGRM.hpp"

namespace {

enum class DispatchTarget {
    none,
    grm,
    process_grm,
    association,
    mom
};

struct DispatchSelection {
    DispatchTarget target = DispatchTarget::none;
    bool conflict = false;
};

bool matches_option(std::string_view arg, std::string_view short_opt, std::string_view long_opt) {
    return arg == short_opt || arg == long_opt;
}

bool is_help_flag(std::string_view arg) {
    return arg == "-h" || arg == "--help";
}

DispatchTarget target_for_argument(std::string_view arg) {
    if (matches_option(arg, "-m", "--make-grm")) {
        return DispatchTarget::grm;
    }
    if (matches_option(arg, "-d", "--process-grm")) {
        return DispatchTarget::process_grm;
    }
    if (arg == "-g" || arg == "-x" ||
            arg == "--test-main" || arg == "--test-main-binary" ||
            arg == "--test-main-binary-continuous" ||
            arg == "--test-main-multitrait-continuous" ||
            arg == "--test-gxe") {
        return DispatchTarget::association;
    }
    if (matches_option(arg, "-o", "--mom")) {
        return DispatchTarget::mom;
    }
    return DispatchTarget::none;
}

DispatchSelection resolve_dispatch_target(int argc, char* argv[]) {
    DispatchSelection selection;
    for (int i = 1; i < argc; ++i) {
        const DispatchTarget arg_target = target_for_argument(argv[i]);
        if (arg_target == DispatchTarget::none) {
            continue;
        }

        if (selection.target == DispatchTarget::none) {
            selection.target = arg_target;
            continue;
        }

        if (selection.target != arg_target) {
            selection.conflict = true;
            return selection;
        }
    }
    return selection;
}

int dispatch_to_target(DispatchTarget target, int argc, char* argv[]) {
    switch (target) {
        case DispatchTarget::grm: {
            Gmatrix grm_runner;
            return grm_runner.run(argc, argv);
        }
        case DispatchTarget::process_grm: {
            ProcessGRM process_grm_runner;
            return process_grm_runner.run(argc, argv);
        }
        case DispatchTarget::association: {
            fastGxE association_runner;
            return association_runner.run(argc, argv);
        }
        case DispatchTarget::mom: {
            MoM mom_runner;
            return mom_runner.run(argc, argv);
        }
        case DispatchTarget::none:
            return 1;
    }

    return 1;
}

int dispatch_help_to_target(DispatchTarget target, const char* program_name) {
    char help_flag[] = "--help";
    char* help_argv[] = {const_cast<char*>(program_name), help_flag, nullptr};
    return dispatch_to_target(target, 2, help_argv);
}

DispatchTarget help_topic_to_target(int help_topic) {
    switch (help_topic) {
        case 1:
            return DispatchTarget::grm;
        case 2:
            return DispatchTarget::process_grm;
        case 3:
        case 4:
        case 6:
            return DispatchTarget::association;
        case 5:
            return DispatchTarget::mom;
        default:
            return DispatchTarget::none;
    }
}

bool parse_help_topic(std::string_view arg, int& help_topic) {
    if (arg.empty()) {
        return false;
    }

    char* parse_end = nullptr;
    const long parsed_value = std::strtol(std::string(arg).c_str(), &parse_end, 10);
    if (parse_end == nullptr || *parse_end != '\0') {
        return false;
    }

    help_topic = static_cast<int>(parsed_value);
    return true;
}

int resolve_help_topic(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (!is_help_flag(argv[i])) {
            continue;
        }

        if (i + 1 < argc) {
            int help_topic = 0;
            if (parse_help_topic(argv[i + 1], help_topic)) {
                return help_topic;
            }
        }
        return 0;
    }

    return -1;
}

void show_help_general() {
    std::cout
        << "Usage:\n"
        << "  fastgxe --make-grm [options]\n"
        << "  fastgxe --process-grm [options]\n"
        << "  fastgxe --test-main [options]\n"
        << "  fastgxe --test-main-binary [options]\n"
        << "  fastgxe --test-main-binary-continuous [options]\n"
        << "  fastgxe --test-main-multitrait-continuous [options]\n"
        << "  fastgxe --test-gxe [options]\n"
        << "  fastgxe --mom [options]\n\n"
        << "Top-Level Modes:\n"
        << "  -m, --make-grm           Compute the genomic relationship matrix (GRM)\n"
        << "  -d, --process-grm        Post-process an existing GRM\n"
        << "  -g, --test-main          Run the association module for SNP main effects\n"
        << "      --test-main-binary   Run the association module for binary traits\n"
        << "      --test-main-binary-continuous\n"
        << "                           Run the joint binary + continuous association module\n"
        << "      --test-main-multitrait-continuous\n"
        << "                           Run the multitrait continuous association module\n"
        << "  -x, --test-gxe           Run the association module for GxE testing\n"
        << "  -o, --mom                Estimate GxE heritability with Method of Moments\n\n"
        << "Help:\n"
        << "  -h, --help               Show this general help\n"
        << "  -h 1                     Show GRM-construction help\n"
        << "  -h 2                     Show GRM post-processing help\n"
        << "  -h 3                     Show association-module help\n"
        << "  -h 5                     Show MoM help\n\n"
        << "Examples:\n"
        << "  fastgxe --make-grm --bfile data --out data_grm\n"
        << "  fastgxe --process-grm --grm data_grm --out data_grm_processed\n"
        << "  fastgxe --test-main --grm data_grm --bfile data --data pheno.txt --trait BMI --out test_main\n"
        << "  fastgxe --test-main-binary --grm data_grm --bfile data --data pheno.txt --trait disease --out test_main_binary\n"
        << "  fastgxe --test-main-binary-continuous --grm data_grm --bfile data --data pheno.txt --trait disease BMI --out test_main_binary_continuous\n"
        << "  fastgxe --test-main-multitrait-continuous --grm data_grm --bfile data --data pheno.txt --trait trait1 trait2 trait3 --out test_main_multitrait_continuous\n"
        << "  fastgxe --test-gxe --grm data_grm --bfile data --data pheno.txt --trait BMI --env-int age smoking --out test_gxe\n"
        << "  fastgxe --mom --grm data_grm --data pheno.txt --trait BMI --env-int age smoking --out test_mom\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc == 1) {
        show_help_general();
        return 1;
    }

    const DispatchSelection selection = resolve_dispatch_target(argc, argv);
    if (selection.conflict) {
        spdlog::error(
            "Choose only one top-level mode among --make-grm, --process-grm, "
            "--test-main/--test-main-binary/--test-main-binary-continuous/"
            "--test-main-multitrait-continuous/--test-gxe, and --mom.");
        show_help_general();
        return 1;
    }

    const int help_topic = resolve_help_topic(argc, argv);
    if (help_topic >= 0) {
        if (help_topic == 0) {
            if (selection.target != DispatchTarget::none) {
                return dispatch_to_target(selection.target, argc, argv);
            }
            show_help_general();
            return 0;
        }

        const DispatchTarget help_target = help_topic_to_target(help_topic);
        if (help_target == DispatchTarget::none) {
            spdlog::error("Invalid argument for -h/--help: {}", help_topic);
            show_help_general();
            return 1;
        }
        return dispatch_help_to_target(help_target, argv[0]);
    }

    if (selection.target == DispatchTarget::none) {
        spdlog::error("No top-level analysis mode was specified.");
        show_help_general();
        return 1;
    }

    return dispatch_to_target(selection.target, argc, argv);
}
