---
layout: page
title: Tutorial
description:
---



# 1. Calculate the Genomic Relationship Matrix (GRM)

## Single-Step Calculation  

For standard datasets, compute the GRM in a single step:  

```
fastgxe --make-grm --bfile test --code-type 2 --out test
```

## Parallel Computation for Large Datasets (e.g., Biobank Data)  

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

## Sample Clustering Based on Relatedness  

To group samples such that any pair with an estimated relatedness greater than 0.05 belongs to the same subgroup:  

```
fastgxe --process-grm --group --grm test.agrm --cut-value 0.05 --out test.agrm
```

---

# 2. Genome-Wide GxE Analysis

Using the previously generated sparse GRM, along with a PLINK binary file and a data file, perform a GxE analysis as follows:  

```
fastgxe --test-gxe --grm test --bfile bed_file --data data_file --trait trait \
--env-int Age:Confide --covar PCA1:PCA5 Sex --out test_gxe
```

- The column `trait` represents the trait values.  
- Columns from `Age` to `Confide` are considered interacting environmental factors.  
- Columns `PCA1` to `PCA5` and `Sex` are treated as confounding covariates.  

## Parallel Computation for Large Datasets  

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

# 3. GxE Heritability Estimation Using the Method of Moments  

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

# 4. Identifying Environmental Factors Driving GxE  

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
