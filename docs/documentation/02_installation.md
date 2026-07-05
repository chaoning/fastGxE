---
layout: page
title: Installation
description: ~
order: 2
---

`fastGxE` is a C++ command-line tool for genome-wide genotype-by-environment analysis. This page describes two installation paths:

- use the prebuilt Linux executable
- build `fastGxE` from source with CMake

## 1. Prebuilt Linux executable

A prebuilt 64-bit Linux executable is available here:
[**fastGxE Linux Executable**](https://github.com/chaoning/fastGxE/raw/refs/heads/main/app/linux/fastgxe)

After downloading it, make it executable and check that it runs:

```bash
chmod +x fastgxe
./fastgxe -h
```

This is the quickest option if you only need to run the program and your system is compatible with the published binary.

## 2. Build from source

Building from source is recommended when you want to modify the code, change compiler settings, or adjust library paths for your local environment.

### Prerequisites

The current build system expects:

- a C++17 compiler
- CMake 3.16 or newer
- Intel oneAPI compiler and Intel MKL
- OpenMP support

The repository already vendors several dependencies under `external/`, including:

- `GSL`
- `Eigen`
- `LBFGSpp`
- `CLI`
- `spdlog`
- `Boost`

### Clone the repository

```bash
git clone https://github.com/chaoning/fastGxE.git
cd fastGxE
```

### Check and update build paths

Before compiling, review `CMakeLists.txt` and update the local paths if needed. In particular, these variables often need to match your machine:

- `CMAKE_CXX_COMPILER`
- `COMPILERROOT`
- `MKLROOT`
- `GSLROOT`
- `EIGENROOT`

The current project configuration is written for an Intel oneAPI-based environment. If you build with a different compiler or install location, adjust the compiler and library paths first.

### Configure and compile

First make the Intel oneAPI toolchain available in your shell (this puts `icpx` and the MKL / OpenMP runtime on the environment). Then configure and build:

```bash
source /path/to/intel/oneapi/setvars.sh
cmake -S . -B build
cmake --build build -j
```

`fastGxE` is C++ only, and the CMake project selects the Intel `icpx` compiler by default (it is set before `project()`), so a plain `cmake -S . -B build` is enough once the oneAPI environment is sourced. To use a different compiler, override it explicitly:

```bash
cmake -S . -B build -DCMAKE_CXX_COMPILER=/path/to/your/compiler
```

The main executable is generated at `build/fastgxe`:

```bash
./build/fastgxe -h
```

### Static linking option

The project currently enables full static linking by default through:

```cmake
FASTGXE_FULL_STATIC_LINK=ON
```

To keep the default static-link behavior:

```bash
cmake .. -DFASTGXE_FULL_STATIC_LINK=ON
cmake --build . -j
```

If static linking is not available on your system, switch to dynamic linking:

```bash
cmake .. -DFASTGXE_FULL_STATIC_LINK=OFF
cmake --build . -j
```

You can inspect the resulting binary with:

```bash
file ./fastgxe
ldd ./fastgxe
```

## 3. Quick validation

After compilation, confirm that the executable starts correctly:

```bash
./fastgxe -h
```

If you want a minimal end-to-end smoke test, the example files in `example/` can be used together with the commands shown in the tutorial page.

## 4. Troubleshooting

- If CMake cannot find MKL or the Intel compiler, recheck the path variables in `CMakeLists.txt`.
- If static linking fails, rebuild with `-DFASTGXE_FULL_STATIC_LINK=OFF`.
- If you changed compiler or dependency paths, it is usually safest to remove the old `build/` directory and configure again from scratch.
- If the executable starts but fails at runtime, verify that the compiler, MKL, and OpenMP settings are consistent with your local environment.
