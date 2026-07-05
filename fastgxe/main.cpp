/*
 * @Description: Top-level CLI dispatcher for the fastGxE toolkit
 * @Author: Chao Ning
 * @Date: 2025-01-30
 * LastEditors: Chao Ning
 * LastEditTime: 2026-04-12 22:40:00
 *
 * What this file does
 * -------------------
 * fastGxE is one executable that contains several independent sub-tools. This
 * file is just the "front desk": it looks at the command line, decides which
 * sub-tool the user wants, and forwards the original arguments to that tool.
 * It does NOT parse the detailed options (each sub-tool does that itself).
 *
 * The sub-tools and the flag that selects each one:
 *
 *     --make-grm      ->  Gmatrix      build a genomic relationship matrix (GRM)
 *     --process-grm   ->  ProcessGRM   post-process an existing GRM
 *     --test-main     ->  fastGxE      SNP main-effect association scan
 *     --test-gxe      ->  fastGxE      SNP-by-environment (GxE) association scan
 *     --mom           ->  MoM          Method-of-Moments GxE heritability
 *
 * --test-main and --test-gxe both run the fastGxE class, so here they count as
 * the SAME choice (Tool::Association); fastGxE::run() tells them apart later.
 *
 * The rule enforced here: the user must pick exactly one tool.
 */

#include <cstdlib>      // std::strtol
#include <iostream>     // std::cout
#include <string_view>  // std::string_view (cheap string comparisons, no copies)

#include <spdlog/spdlog.h>  // logging for error messages

#include "fastgxe.hpp"     // fastGxE     (--test-main / --test-gxe)
#include "gmatrix.hpp"     // Gmatrix     (--make-grm)
#include "mom.hpp"         // MoM         (--mom)
#include "processGRM.hpp"  // ProcessGRM  (--process-grm)

namespace {  // everything below is private to this file

// ---------------------------------------------------------------------------
// Which sub-tool to run
// ---------------------------------------------------------------------------

// The five possible choices. `None` means "no tool was requested".
enum class Tool {
    None,
    MakeGrm,
    ProcessGrm,
    Association,  // shared by --test-main and --test-gxe
    Mom
};

// Looks at ONE command-line word and returns the tool it selects (or None if
// the word is not a tool flag, e.g. a file name or a value).
Tool tool_for_flag(std::string_view arg) {
    if (arg == "-m" || arg == "--make-grm")     return Tool::MakeGrm;
    if (arg == "-d" || arg == "--process-grm")  return Tool::ProcessGrm;
    if (arg == "-g" || arg == "--test-main" ||
        arg == "-x" || arg == "--test-gxe")     return Tool::Association;
    if (arg == "-o" || arg == "--mom")          return Tool::Mom;
    return Tool::None;
}

// Outcome of scanning the whole command line for a tool flag.
struct ToolChoice {
    Tool tool = Tool::None;
    bool conflict = false;  // true if two DIFFERENT tools were both requested
};

// Walks every argument and decides which single tool was asked for.
//   - words that are not tool flags are ignored
//   - the first tool flag wins
//   - a second, DIFFERENT tool flag is an error (conflict = true)
//   - repeating the same tool flag is fine
ToolChoice choose_tool(int argc, char* argv[]) {
    ToolChoice choice;
    for (int i = 1; i < argc; ++i) {  // start at 1: argv[0] is the program name
        const Tool t = tool_for_flag(argv[i]);
        if (t == Tool::None) continue;          // not a tool flag

        if (choice.tool == Tool::None) {
            choice.tool = t;                    // remember the first tool
        } else if (choice.tool != t) {
            choice.conflict = true;             // a different tool -> conflict
            return choice;
        }
    }
    return choice;
}

// Creates the requested tool object, runs it with the ORIGINAL command line,
// and returns its exit code. Each tool parses its own options inside run().
int run_tool(Tool tool, int argc, char* argv[]) {
    switch (tool) {
        case Tool::MakeGrm:     return Gmatrix{}.run(argc, argv);
        case Tool::ProcessGrm:  return ProcessGRM{}.run(argc, argv);
        case Tool::Association: return fastGxE{}.run(argc, argv);
        case Tool::Mom:         return MoM{}.run(argc, argv);
        case Tool::None:        return 1;  // should never happen; fail safely
    }
    return 1;
}

// ---------------------------------------------------------------------------
// Help handling
// ---------------------------------------------------------------------------

// What kind of help (if any) the user asked for.
//   asked == false           -> no -h/--help on the command line
//   asked == true,  topic 0  -> general help (bare "-h", or "-h" + non-number)
//   asked == true,  topic >0 -> help for a specific numbered topic ("-h 3")
struct HelpRequest {
    bool asked = false;
    int topic = 0;
};

// True for "-h" or "--help".
bool is_help_flag(std::string_view arg) {
    return arg == "-h" || arg == "--help";
}

// Tries to read `arg` as a whole base-10 integer (e.g. the "3" in "-h 3").
// Succeeds only if the ENTIRE word is a number ("3" yes, "3x" no, "abc" no).
bool try_parse_int(std::string_view arg, int& out) {
    if (arg.empty()) return false;
    char* end = nullptr;
    const long value = std::strtol(std::string(arg).c_str(), &end, 10);
    if (end == nullptr || *end != '\0') return false;  // leftover characters -> not a pure int
    out = static_cast<int>(value);
    return true;
}

// Scans the command line for the first help flag and, if present, the numeric
// topic that follows it. See HelpRequest above for the meaning of the fields.
//
// Note (kept identical to the original tool): a help flag followed by a number
// <= 0 (e.g. "-h 0" or "-h -3") is treated as general help.
HelpRequest scan_help(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (!is_help_flag(argv[i])) continue;

        HelpRequest help;
        help.asked = true;

        // If the next word is a positive number, treat it as a topic.
        int topic = 0;
        if (i + 1 < argc && try_parse_int(argv[i + 1], topic) && topic > 0) {
            help.topic = topic;
        }
        return help;  // only the first help flag matters
    }
    return HelpRequest{};  // asked == false
}

// Maps a numeric help topic to the tool whose help should be shown.
// Topics 3, 4 and 6 all refer to the association module. Unknown topics
// return Tool::None so the caller can report an error.
Tool tool_for_help_topic(int topic) {
    switch (topic) {
        case 1:            return Tool::MakeGrm;
        case 2:            return Tool::ProcessGrm;
        case 3: case 4: case 6: return Tool::Association;
        case 5:            return Tool::Mom;
        default:           return Tool::None;
    }
}

// Prints a tool's OWN detailed help by running it with a synthetic "--help"
// command line, so each tool formats its own help text.
int print_tool_help(Tool tool, const char* program_name) {
    char help_flag[] = "--help";
    char* help_argv[] = {const_cast<char*>(program_name), help_flag, nullptr};
    return run_tool(tool, 2, help_argv);
}

// Prints the program-wide usage banner (the list of tools and a few examples).
void print_general_help() {
    std::cout
        << "Usage:\n"
        << "  fastgxe --make-grm [options]\n"
        << "  fastgxe --process-grm [options]\n"
        << "  fastgxe --test-main [options]\n"
        << "  fastgxe --test-gxe [options]\n"
        << "  fastgxe --mom [options]\n\n"
        << "Top-Level Modes:\n"
        << "  -m, --make-grm           Compute the genomic relationship matrix (GRM)\n"
        << "  -d, --process-grm        Post-process an existing GRM\n"
        << "  -g, --test-main          Run the association module for SNP main effects\n"
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
        << "  fastgxe --test-gxe --grm data_grm --bfile data --data pheno.txt --trait BMI --env-int age smoking --out test_gxe\n"
        << "  fastgxe --mom --grm data_grm --data pheno.txt --trait BMI --env-int age smoking --out test_mom\n";
}

}  // namespace

// ---------------------------------------------------------------------------
// Entry point: a straight sequence of guard clauses, top to bottom.
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // 1. No arguments at all -> show help and report misuse.
    if (argc == 1) {
        print_general_help();
        return 1;
    }

    // 2. Decide which tool was requested; reject mixing two tools.
    const ToolChoice choice = choose_tool(argc, argv);
    if (choice.conflict) {
        spdlog::error(
            "Choose only one top-level mode among --make-grm, --process-grm, "
            "--test-main/--test-gxe, and --mom.");
        print_general_help();
        return 1;
    }

    // 3. Help requested? Handle it and stop here.
    const HelpRequest help = scan_help(argc, argv);
    if (help.asked) {
        if (help.topic > 0) {
            // "-h N": show that topic's tool-specific help.
            const Tool help_tool = tool_for_help_topic(help.topic);
            if (help_tool == Tool::None) {
                spdlog::error("Invalid argument for -h/--help: {}", help.topic);
                print_general_help();
                return 1;
            }
            return print_tool_help(help_tool, argv[0]);
        }

        // Bare "--help": if a tool was also named (e.g. "--test-gxe --help"),
        // let that tool print its own help; otherwise show the general help.
        if (choice.tool != Tool::None) {
            return run_tool(choice.tool, argc, argv);
        }
        print_general_help();
        return 0;
    }

    // 4. No help and no tool -> nothing to do.
    if (choice.tool == Tool::None) {
        spdlog::error("No top-level analysis mode was specified.");
        print_general_help();
        return 1;
    }

    // 5. Normal case: run the chosen tool with the original arguments.
    return run_tool(choice.tool, argc, argv);
}
