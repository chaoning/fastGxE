---
layout: page
title: Tutorial
description:
---



## 1. Data input

Example files are provided in the **example** directory:
 `test.bed`, `test.bim`, `test.fam`, and `test_simu_pheno.txt`.

**Note:** The example dataset is intentionally small and should only be used for testing program functionality.
 To obtain **stable and reliable results**, analyses should be performed on **large-scale datasets**.

### PLINK Binary File

PLINK binary files consist of three components:

- **`*.bed` (binary genotype file)**: Stores genotype data in a compact binary format.
- **`*.bim` (variant information file)**: Contains SNP identifiers, chromosome positions, and allele information.
- **`*.fam` (sample information file)**: Includes individual IDs, family structure, sex, and phenotype data.

In fastGxE, missing genotype values for the focus SNP are imputed using the mean genotype. However, it is recommended to preprocess missing genotypes with Beagle or other imputation software for improved accuracy.



### Data File Format

The data file should contain individual IDs, confounding covariates, interacting environmental covariates, and phenotypic values. It must follow the structure below:

1. **Header Line**: A header row is required.
2. **Delimiter**: Data should be separated by blanks or tabs.
3. **First Column**: The first column must be **individual IDs**.
4. **Matching Individuals**:
   - Only individuals present in **both** the PLINK binary file and the data file will be used in the analysis.
   - Individuals missing in either the PLINK binary file or the data file will be automatically removed by fastGxE during the analysis.

An example data file

| iid     | pheno    | E1       | E2       | E3       | E4       | E5       | E6       | E7       | E8       | E9       | E10      | E11      | E12      | E13      | E14      | E15      | E16      | E17      | E18      | E19      | E20      | E21      | E22      | E23      | E24      | E25      | E26      | E27      | E28      | E29      | E30      | E31      | E32      | E33      | E34      | E35      | E36      | E37      | E38      | E39      | E40      |
| ------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| 1110005 | -0.75286 | 0.21092  | 0.648552 | 0.468537 | 0.285131 | -0.15092 | -0.1555  | -0.56958 | 0.07954  | -0.1891  | 0.338504 | -0.07503 | 0.744465 | -1.2262  | -0.77347 | -0.9495  | 0.204434 | 0.287995 | -0.36906 | 0.023615 | 0.100999 | 0.913586 | 0.281957 | -0.26894 | 1.076054 | 0.631956 | -0.25934 | -0.12063 | 0.293751 | -0.14226 | 0.269535 | -0.69097 | 0.12688  | 1.172146 | -0.77402 | -1.75724 | -2.34637 | -0.57211 | 1.166077 | -0.4581  | -0.55331 |
| 1110006 | -0.62861 | 0.000148 | -0.12916 | -1.77284 | -0.31648 | 0.461651 | 0.281532 | -0.74519 | -1.43904 | -0.94778 | -0.1699  | 0.276132 | 0.292023 | -0.52048 | 0.963399 | 0.311522 | -1.27339 | 0.302758 | -0.91646 | 1.055273 | 0.500981 | -0.87986 | -0.48948 | 0.78681  | 1.010843 | -1.76933 | 1.697091 | -0.32637 | -0.83483 | -1.82932 | -1.17107 | 0.733051 | -0.82921 | -0.22705 | 0.509752 | -0.06756 | 1.728553 | -0.03009 | -1.43085 | 0.59895  | 0.836599 |
| 1110008 | 1.132841 | -1.14961 | -0.29059 | 0.446468 | -0.15161 | 0.222543 | -0.33517 | 0.81958  | 0.502829 | -1.87353 | -0.82221 | 0.332078 | 1.77179  | -0.12034 | -0.81934 | 0.500914 | 0.265973 | -0.81532 | 0.794444 | 1.517311 | -1.1751  | 0.228761 | -0.07199 | 0.446747 | 0.091397 | -1.82258 | 0.542994 | -0.16897 | 0.352514 | -1.63777 | -0.10023 | 1.466236 | 0.746707 | 0.59923  | 0.515903 | 0.473775 | 1.697436 | 1.165627 | 0.500921 | 1.237719 | -0.10859 |
| 1110009 | 0.773396 | -0.90416 | -0.00513 | -0.045   | -0.90065 | -0.46537 | 0.752511 | 0.355544 | 0.359747 | 1.585845 | -0.38118 | -0.50373 | 0.253581 | -0.44434 | -0.5287  | -0.15618 | 1.483534 | 1.320471 | 0.593645 | 0.506708 | -0.65951 | 0.033098 | -1.52481 | 1.100987 | 0.33127  | -0.45424 | -0.18036 | 1.700584 | 1.22     | -0.86434 | -0.15836 | 0.332595 | 1.73738  | 1.11772  | -0.92361 | -0.39209 | -3.61622 | 0.128307 | -0.20065 | 0.291977 | -1.28622 |
| 1110011 | 0.209078 | -0.87784 | -2.24717 | -0.33057 | -0.87823 | 1.68817  | 0.193299 | 0.753634 | -0.22589 | -2.1622  | 1.294862 | 0.898332 | -0.33214 | -1.23525 | -1.06617 | -1.29512 | 1.44334  | -1.03666 | -0.44052 | -1.22453 | -0.03294 | 1.713709 | -0.53588 | -0.11923 | -0.68673 | 0.378322 | -1.08003 | 0.897898 | -2.00473 | 1.621039 | -1.26384 | 0.352975 | -0.45193 | -0.36424 | -0.07237 | -0.4594  | 0.971471 | -0.84527 | 0.207248 | 2.920991 | -0.81485 |
| 1110023 | -0.13773 | 0.99387  | 0.577044 | 0.074119 | -0.34685 | -1.34186 | 0.63063  | 0.828125 | 1.127962 | 0.679971 | 0.242218 | -0.06751 | -0.13105 | -1.01052 | -0.86894 | -1.46575 | 0.20377  | -0.93781 | -1.36779 | 0.476924 | -1.15502 | -0.0056  | 0.40834  | -1.19064 | -0.12594 | -0.13059 | -0.50763 | 0.904232 | -1.18535 | -0.87793 | -1.05899 | -0.01426 | 1.457778 | 1.412801 | 0.00192  | 0.074879 | -0.10291 | 1.164542 | 0.203062 | 0.478138 | 1.223259 |
| 1110033 | -1.16732 | 0.14276  | -0.4078  | 0.215731 | 1.581048 | 0.077753 | -0.12957 | -0.91409 | -1.17939 | -0.74846 | 1.230236 | -0.41687 | 0.228574 | 0.667259 | 1.042392 | 0.16277  | -0.22133 | 1.092244 | -0.20791 | 0.527073 | 0.03231  | -0.46105 | 0.084871 | -1.07208 | 0.936474 | 0.117742 | 0.472216 | 0.086622 | 0.072974 | 0.02184  | -0.47137 | 0.089335 | -1.54968 | -0.10841 | 1.040897 | -0.19092 | -0.95165 | 0.342614 | 0.174287 | 0.104902 | 0.20719  |
| 1110035 | 1.52892  | 0.380602 | -0.23611 | -0.25621 | 1.251567 | -1.07136 | 1.600939 | -1.5843  | 0.54906  | 1.703418 | -1.6009  | 0.126096 | -0.13079 | -0.30247 | -1.79334 | 0.608061 | -1.37453 | -0.59196 | -1.2016  | 0.451544 | 0.170108 | -1.18104 | 0.517898 | -0.84502 | 1.49383  | 0.509155 | -0.28951 | 0.18437  | -0.27631 | 0.032461 | 2.121288 | 1.155691 | -0.84489 | 0.037005 | -0.19401 | -0.7891  | -1.29008 | -0.49291 | -0.82919 | -0.94171 | 0.200087 |
| 1110043 | -1.63005 | 1.625271 | 1.558561 | -0.36283 | -1.31468 | -0.38357 | 0.820452 | 1.067546 | 1.625511 | 1.289687 | -1.14444 | 1.275622 | 1.990491 | -0.69936 | -0.08353 | -0.63966 | 2.483601 | -1.80873 | 0.815655 | 0.134742 | 1.330531 | -0.32461 | -0.33861 | 0.008995 | -0.55535 | 1.039129 | -0.01007 | -0.38419 | 0.107745 | -2.12606 | -0.17434 | 1.990626 | 4.504256 | 1.592914 | 0.327584 | -0.29862 | -1.24944 | 1.32345  | -0.89109 | 1.247789 | 2.851869 |
| 1110044 | 1.27608  | -1.35986 | -2.07641 | -1.65412 | 0.267007 | 0.046312 | -1.33002 | -0.13738 | -0.05922 | -0.651   | -0.0998  | -0.3337  | 0.188556 | -0.72331 | -0.84339 | 0.091237 | -0.23129 | 0.493803 | -1.79401 | -0.75501 | -0.30332 | -0.32378 | 0.658686 | 1.003795 | -0.26546 | 1.366845 | 1.35157  | -0.17636 | 0.769366 | 0.844615 | -1.59522 | 0.466167 | -0.83507 | -0.00921 | -0.15057 | 0.169063 | 0.677308 | 0.109467 | 0.440098 | 0.611541 | -1.35669 |
| 1110063 | -0.40119 | 1.028213 | 0.746369 | 0.000715 | 0.037617 | 0.65775  | 0.234105 | 0.660881 | -0.09425 | -0.22022 | 0.299793 | -0.31184 | -1.48914 | -0.82789 | -0.88923 | -1.36938 | 0.712816 | 0.124924 | -0.44722 | -1.54906 | 2.592713 | -0.25346 | -1.10651 | 0.757067 | -0.87045 | -1.89933 | 0.686628 | -0.25726 | -1.73575 | -0.51228 | 1.107235 | 1.016678 | -0.60243 | 0.017886 | -2.17325 | -0.19513 | -0.41595 | -1.96924 | -1.52661 | -0.19715 | -0.65824 |

**Associated options in fastGxE**

<font color=red>**--missing-data**</font> NA Na na NAN NaN nan -NAN -NaN -nan N/A n/a

Any symbol listed above will be treated as a missing value. Symbols are separated by spaces, and users can add custom missing values to the list. If a variable is included in the model and contains a missing value in this field, the entire row will be discarded.

## 2. Calculate the Genomic Relationship Matrix (GRM)

### Single-Step Calculation  

For standard datasets, compute the GRM in a single step:  

```
fastgxe --make-grm --bfile test --code-type 2 --out test
```

### Parallel Computation for Large Datasets (e.g., Biobank Data)  

For large-scale datasets, partition the GRM into 100 parts, compute them in parallel, and then merge the results:  

```
fastgxe --make-grm --bfile test --code-type 2 --npart 100 1 --out test
fastgxe --make-grm --bfile test --code-type 2 --npart 100 2 --out test
...
fastgxe --make-grm --bfile test --code-type 2 --npart 100 100 --out test
```

Once all parts are computed, merge them into a single GRM file:  

```
fastgxe --process-grm --merge --grm test.agrm --npart 100
```

### Sample Clustering Based on Relatedness  

To group samples such that any pair with an estimated relatedness greater than 0.05 belongs to the same subgroup:  

```
fastgxe --process-grm --group --grm test.agrm --cut-value 0.05 --out test.agrm
```

---

## 3. Genome-Wide GxE Analysis

Using the previously generated sparse GRM, along with a PLINK binary file and a data file, perform a GxE analysis as follows:  

```
fastgxe --test-gxe --grm test --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_gxe
```

- The column `pheno` represents the trait values.  
- Columns from `E1` to `E40` are considered interacting environmental factors.  

### Parallel Computation for Large Datasets  

For large datasets, partition the SNPs into 10 parts and compute them in parallel. The results can then be merged using a custom Python or R script.  

```
fastgxe --test-gxe --split-task 10 1 --grm test --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_gxe 

fastgxe --test-gxe --split-task 10 2 --grm test --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_gxe 

...

fastgxe --test-gxe --split-task 10 10 --grm test --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_gxe
```

Once all tasks are completed, the results can be merged using a custom Python or R script.

---

## 4. GxE Heritability Estimation Using the Method of Moments  

To obtain an unbiased estimate of GxE heritability, use the method of moments:  

Note: we need large sample size to get stable results. The example data have too small size.

```
fastgxe --mom --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_mom
```

---

## 5. Identifying Environmental Factors Driving GxE  

Once GxE genomic loci are identified, the environmental factors driving the GxE effect of the lead SNP can be determined using **mmSuSiE**.  mmSuSiE is written in python can be found at https://github.com/chaoning/mmsusie

First reformat the genomic relationship matrix with fastgxe

```
# This step generates the reformatted GRM used by mmSuSiE.
fastgxe --process-grm --reformat --grm test.agrm --out-fmt 1 --out test.agrm
```



```python
from mmsusie import MMSuSiE
model = MMSuSiE()

# Define input files
bed_file = "test"              # PLINK binary genotype file prefix (.bed/.bim/.fam)
grm_file = "test"              # GRM (genetic relationship matrix) file prefix
pheno_file = "test_simu_pheno.txt"  # Phenotype file with environmental covariates and trait

# Define a list of interacting environmental variables: E1, E2, ..., E40
env_int = [f"E{i}" for i in range(1, 41)]

# Perform LD pruning to identify leading SNPs for GxE signals
assoc_file = "test_gxe.res"    # Association results from fastGxE or similar method
df_leading = model.ld_pure(assoc_file, bed_file, ld_r2=0.1, snp="SNP", p="p_gxe", p_cutoff=5e-8)
print(df_leading)              # Print the pruned set of lead SNPs

# Specify parameters for mmSuSiE fine-mapping
snp_id = "rs550011"            # Lead SNP to analyze
trait = "pheno"                # Target phenotype column name
varcom_file = "test_gxe.var"   # Variance components file estimated from fastGxE --test-gxe
out_file = "test_mmsusie"      # Output file prefix for mmSuSiE results

# Run the mmSuSiE model to identify environmental factors driving the GxE effect
res_dct = model.run(pheno_file, trait, env_int, grm_file, bed_file, snp_id, varcom_file, out_file)

# Print the credible sets (CS) identified by mmSuSiE (index start from 0)
print(res_dct['cs'])

# Extract posterior inclusion probabilities (PIPs)
df_pip = res_dct["pip"]
print(df_pip)

# Display environmental factors with PIP > 0.5
print(df_pip[df_pip["pip"] > 0.5])

```
