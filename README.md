# fastGxE: **Powering genome-wide detection of genotype-environment interactions in biobank studies**



## Overview
![scheme](./docs/images/overview.v1.0.png)

fastGxE is a scalable and effective method designed for genome-wide GxE association analysis. fastGxE handles multiple environmental factors and examines one SNP at a time, decomposing the phenotype into SNP main effect, environmental main effects, GxE interaction effects, while controlling for both polygenic effects and polygenic interaction effects (a). fastGxE evaluates various GxE effect size configurations and combine the resulting p-values into a single p-value to test whether the SNP interacts with at least one environmental factor (b). With controlled polygenic and polygenic interaction effects, fastGxE generates calibrated p-values for identifying candidate GxE loci (c). Additionally, it utilizes mmSuSiE, an extension of the SuSiE algorithm, to identify the environmental factors driving the detected GxE interactions and employs the stratified Wald test to validate and visualize these interactions (d).

## Installation

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

## How to use `fastGxE`



Check our [vignettes](https://chaoning.github.io/fastGxE).

## Citing the work



If you find the `fastGxE` package or any of the source code in this repository useful for your work, please cite:

> Chao Ning and Xiang Zhou (2025). Powering genome-wide detection of genotype-environment interactions in biobank studies.



Visit our [group website](https://xiangzhou.github.io/) for more statistical tools on analyzing genetics, genomics and transcriptomics data.

