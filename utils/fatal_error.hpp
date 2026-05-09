#pragma once

#include <cstdlib>
#include <spdlog/spdlog.h>

template<typename... Args>
[[noreturn]] inline void fatal_error(spdlog::format_string_t<Args...> fmt_str, Args&&... args) {
    spdlog::error(fmt_str, std::forward<Args>(args)...);
    std::exit(1);
}
