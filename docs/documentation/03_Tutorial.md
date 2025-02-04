---
layout: page
title: Tutorial
description:
---



## 1. Data input

### PLINK Binary File

PLINK binary files consist of three components:

- **`*.bed` (binary genotype file)**: Stores genotype data in a compact binary format.
- **`*.bim` (variant information file)**: Contains SNP identifiers, chromosome positions, and allele information.
- **`*.fam` (sample information file)**: Includes individual IDs, family structure, sex, and phenotype data.

In fastGxE, missing genotype values for the focus SNP are imputed using the mean genotype. However, it is recommended to preprocess missing genotypes with Beagle or other mputation software for improved accuracy.



### Data File Format

The data file should contain individual IDs, confounding covariates, interacting environmental covariates, and phenotypic values. It must follow the structure below:

1. **Header Line**: A header row is required.
2. **Delimiter**: Data should be separated by blanks or tabs.
3. **First Column**: The first column must be **individual IDs**.
4. **Matching Individuals**:
   - Only individuals present in **both** the PLINK binary file and the data file will be used in the analysis.
   - Individuals missing in either the PLINK binary file or the data file will be automatically removed by fastGxE during the analysis.

An example data file

| eid  | trait | PCA1     | PCA2    | PCA3     | PCA4     | PCA5     | Sex  | Age  | TDI      | Stair_climbing | Moderate_PA | Vigorous_PA | Walked   | Driving  | Sleep_duration2 | Using_computer | Watching_TV | Walking_pace | Phone_use | Computer_games | Sleep_duration | Getting_up | Nap      | Sleeplessness | Daytime_dozing | Smoking  | Oily_fish | Non-oily_fish | Processed_meat | Poultry  | Beef     | Lamb     | Pork     | Cheese   | Added_salt | Hot_drink | Variation_diet | Cooked_vegetable | Raw_vegetable | Fresh_fruit | Dried_fruit | Bread    | Cereal   | Tea      | Coffee   | Water    | Alcohol  | Friend_visits | Confide  |
| ---- | ----- | -------- | ------- | -------- | -------- | -------- | ---- | ---- | -------- | -------------- | ----------- | ----------- | -------- | -------- | --------------- | -------------- | ----------- | ------------ | --------- | -------------- | -------------- | ---------- | -------- | ------------- | -------------- | -------- | --------- | ------------- | -------------- | -------- | -------- | -------- | -------- | -------- | ---------- | --------- | -------------- | ---------------- | ------------- | ----------- | ----------- | -------- | -------- | -------- | -------- | -------- | -------- | ------------- | -------- |
| 1    | 110   | -12.2135 | 3.72002 | -0.65186 | 0.273789 | -10.0724 | 0    | 35   | -0.66726 | -0.98706       | 1.462283    | -0.9705     | 0.828779 | -0.43317 | -0.18776        | -0.55128       | -0.46106    | 1.068678     | -2.26352  | -0.46436       | 0.842472       | 1.158258   | -0.79583 | -1.40547      | -0.51618       | -0.88689 | -0.74539  | 0.26308       | 0.093726       | -0.37088 | -0.51403 | -0.19141 | 1.219972 | 0.372192 | -0.71954   | 0.001141  | 0.479207       | -0.45328         | -0.65821      | 0.702135    | -0.2159     | -0.38355 | -1.73875 | 1.04156  | -1.11171 | -1.26658 | -1.55927 | -0.70583      | -1.47168 |
| 2    | 109   | -12.1099 | 3.53791 | -0.5843  | 0.604633 | -2.46493 | 1    | 40   | -0.0521  | 0.588444       | -0.28367    | 0.084353    | -2.36427 | 0.133726 | -0.62812        | 0.733788       | -1.798      | 1.068678     | 0.033876  | -0.46436       | -0.17057       | -0.15634   | 0.924186 | 1.369392      | 1.568155       | -0.88689 | 1.504599  | 0.26308       | -1.83878       | -2.70955 | -1.79651 | -1.67671 | -1.70418 | 0.372192 | -0.71954   | -1.68802  | 0.479207       | -0.45328         | 1.169738      | 2.25021     | -0.64274    | -0.25194 | -0.98    | 0.234354 | -1.11171 | 0.252315 | -2.25141 | 2.063359      | 0.771595 |
| 3    | 102   | -12.3776 | 2.45864 | -1.43371 | 3.91774  | 2.92263  | 0    | 20   | 0.817235 | -0.19931       | -0.28367    | -0.44307    | -0.23557 | 0.133726 | -0.18776        | -0.12292       | 1.544351    | 1.068678     | 0.799674  | -0.46436       | 0.842472       | 1.158258   | 0.924186 | -1.40547      | -0.51618       | 1.127516 | 0.379604  | 0.26308       | -0.87253       | 0.798463 | 2.050936 | -0.19141 | -0.24211 | 0.372192 | 1.748492   | 0.001141  | -1.21437       | -0.45328         | -0.04889      | -0.84594    | 0.210932    | -1.30488 | -0.60062 | -1.38006 | 1.675894 | 0.252315 | -0.17499 | -1.62889      | 0.771595 |
| 4    | 101   | -12.8004 | 3.787   | -4.67861 | 4.58364  | 3.23267  | 1    | 36   | -0.01894 | -0.19931       | 1.462283    | 1.666632    | 0.828779 | 0.133726 | 4.676591        | -0.97963       | -1.12953    | -0.63705     | 0.033876  | -0.46436       | 2.86856        | 1.158258   | 0.924186 | -1.40547      | -0.51618       | -0.88689 | 0.379604  | 0.26308       | 0.093726       | -0.37088 | 0.768455 | 1.293892 | -0.24211 | -0.58183 | -0.71954   | 0.001141  | 0.479207       | -1.58983         | 0.560423      | -1.61998    | -0.64274    | 0.274537 | -0.98    | 0.637957 | -1.11171 | -0.76028 | -2.25141 | 0.217236      | -1.47168 |
| 5    | 96    | -11.218  | 4.20388 | -1.86853 | 0.805394 | 5.23192  | 1    | 70   | -0.32468 | -0.19931       | -0.28367    | 1.139206    | 0.828779 | -0.43317 | 1.580476        | -0.12292       | 0.20741     | -0.63705     | 0.033876  | -0.46436       | 1.855516       | -0.15634   | -0.79583 | -0.01804      | -0.51618       | 1.127516 | -1.87038  | -2.37842      | -1.83878       | -2.70955 | -1.79651 | -1.67671 | -1.70418 | 0.372192 | -0.71954   | 0.001141  | 0.479207       | -0.45328         | 1.779053      | 1.476172    | -0.2159     | -1.50231 | -1.54906 | 0.234354 | 1.118373 | -0.25398 | 0.517155 | 1.140297      | -0.35004 |
| 6    | 103   | -11.4399 | 6.39207 | -2.76728 | 1.92316  | -5.91012 | 0    | 30   | -0.52728 | 1.376196       | -0.72016    | 0.084353    | -0.76775 | -0.43317 | -0.62812        | -0.12292       | -1.12953    | -0.63705     | 0.799674  | -0.46436       | -0.17057       | -0.15634   | -0.79583 | -1.40547      | -0.51618       | 1.127516 | -0.74539  | -1.05767      | 0.093726       | 0.798463 | 2.050936 | 1.293892 | 1.219972 | 0.372192 | -0.71954   | 0.001141  | 0.479207       | -0.45328         | 0.560423      | -1.61998    | -0.64274    | -0.25194 | -1.73875 | -0.16925 | -0.55419 | 1.264911 | -0.17499 | 0.217236      | -0.35004 |
| 7    | 114   | -13.5546 | 2.63848 | -1.49553 | 4.12438  | 11.761   | 1    | 45   | 0.412038 | -0.19931       | 1.462283    | 1.139206    | 0.828779 | -0.43317 | -0.18776        | -0.55128       | 0.20741     | 1.068678     | 0.033876  | -0.46436       | 0.842472       | 1.158258   | -0.79583 | -0.01804      | -0.51618       | -0.88689 | 0.379604  | 0.26308       | -1.83878       | -2.70955 | -1.79651 | -1.67671 | -1.70418 | 0.372192 | 0.514474   | 1.690297  | -1.21437       | -1.21098         | -0.65821      | 1.476172    | 0.210932    | 0.537773 | 0.916875 | -0.16925 | -0.55419 | -1.01343 | -0.17499 | -0.70583      | 0.771595 |
| 8    | 110   | -9.17969 | 2.43813 | -0.15923 | 1.56791  | -6.39978 | 0    | 27   | -1.05404 | -0.19931       | -1.15664    | -0.44307    | -0.23557 | -0.43317 | -0.62812        | -0.55128       | 1.544351    | -0.63705     | 0.033876  | -0.46436       | -0.17057       | -0.15634   | -0.79583 | 1.369392      | 1.568155       | 1.127516 | -0.74539  | -1.05767      | 1.059981       | 0.798463 | -0.51403 | -0.19141 | 1.219972 | 0.372192 | 0.514474   | 0.001141  | 0.479207       | -0.45328         | -0.65821      | -0.0719     | 0.210932    | -1.30488 | 0.916875 | 0.637957 | -0.55419 | -0.25398 | -0.17499 | 0.217236      | -2.03249 |
| 9    | 121   | -5.70495 | 2.64119 | 4.12728  | -22.8834 | 2.91856  | 0    | 55   | 1.472918 | 0.588444       | -0.28367    | 0.084353    | -0.23557 | 1.267524 | 0.259408        | 0.733788       | 0.20741     | -0.63705     | 0.799674  | -0.46436       | -1.18362       | -1.47094   | 0.924186 | 1.369392      | 1.568155       | -0.88689 | 1.504599  | -1.05767      | 0.093726       | 0.798463 | -0.51403 | -0.19141 | -0.24211 | 0.372192 | 0.514474   | 0.001141  | 0.479207       | -0.45328         | -0.04889      | -0.0719     | -0.64274    | -1.17326 | 0.5375   | -0.16925 | -1.11171 | 0.252315 | 0.517155 | -0.70583      | -0.35004 |

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
fastgxe --test-gxe --grm test --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_gxe
```

- The column `trait` represents the trait values.  
- Columns from `Age` to `Confide` are considered interacting environmental factors.  
- Columns `PCA1` to `PCA5` and `Sex` are treated as confounding covariates.  

### Parallel Computation for Large Datasets  

For large datasets, partition the SNPs into 10 parts and compute them in parallel. The results can then be merged using a custom Python or R script.  

```
fastgxe --test-gxe --split-task 10 1 --grm test --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_gxe 

fastgxe --test-gxe --split-task 10 2 --grm test --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_gxe 

...

fastgxe --test-gxe --split-task 10 10 --grm test --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_gxe 
```

Once all tasks are completed, the results can be merged using a custom Python or R script.

---

## 4. GxE Heritability Estimation Using the Method of Moments  

To obtain an unbiased estimate of GxE heritability, use the method of moments:  

```
fastgxe --mom --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_mom
```

- The column `trait` represents the trait values.  
- Columns from `Age` to `Confide` are considered interacting environmental factors.  
- Columns `PCA1` to `PCA5` and `Sex` are treated as confounding covariates.  

This approach ensures an unbiased estimation of GxE heritability using the method of moments.

---

## 5. Identifying Environmental Factors Driving GxE  

Once GxE genomic loci are identified, the environmental factors driving the GxE effect of the lead SNP can be determined using **mmSuSiE**.  

```
fastgxe --mmsusie --bfile bed_file --snp-focus snp_name --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_mmsusie
```

- `snp_name` refers to the lead SNP identified in the GxE analysis.  
- The column `trait` represents the trait values.  
- Columns from `Age` to `Confide` are treated as interacting environmental factors.  
- Columns `PCA1` to `PCA5` and `Sex` are considered confounding covariates.  

This approach helps pinpoint the environmental factors driving the GxE effect of the lead SNP.
