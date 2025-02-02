---
layout: page
title: Installation
description: ~
---

`fastGxE` is implemented as an C++ package, which can be installed from GitHub by:

### Linux Executable

A statically compiled executable for 64-bit Linux systems is available: [**fastGxE Linux Executable**](https://github.com/chaoning/fastGxE/raw/refs/heads/main/app/linux/fastgxe). This can be used directly on compatible systems.

```
chmod +x fastgxe
./fastgxe -h
```



### Prerequisites

Ensure the following dependencies are installed on your system:

- **C++ Compiler** (GCC 9+ or Intel C++ Compiler)
- **CMake** (Version 3.16 or higher)
- **Intel MKL** (2024.1)
- **GSL** (2.7)
- **Eigen** (3.4.0)
- **LBFGSpp**
- **OpenMP**
- **CLI**
- **spdlog**

### Installation Steps

#### 1. Clone the Repository

```bash
git clone https://github.com/chaoning/fastGxE.git
cd fastGxE
```

#### 2. Set Up Dependencies

Modify `CMakeLists.txt` to update the paths of external libraries (MKL, GSL, Eigen, LBFGSpp, etc.) according to your system. 

#### 3. Build the Project

Run the following commands to compile fastGxE:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

#### 4. Run fastGxE

After compilation, you can execute FastGxE:

```bash
./fastgxe -h
```

### Notes

- If you encounter missing library errors, check that all paths in `CMakeLists.txt` are correctly configured.
- Use `make clean && make` to rebuild after modifications.
