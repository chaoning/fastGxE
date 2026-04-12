---
layout: page
title: Tutorial
description:
---

This tutorial walks through the example data bundled with `fastGxE` and follows the current command-line workflow implemented in the codebase.

Unless noted otherwise, the commands below assume:

- you are running from the repository root
- the executable is available as `fastgxe`
- example inputs are referenced with the `example/` prefix

If you built from source and the binary is not on your `PATH`, replace `fastgxe` with `./build/fastgxe`.

## 1. Example inputs

The `example/` directory contains a minimal dataset for testing the software:

- `example/test.bed`
- `example/test.bim`
- `example/test.fam`
- `example/test_simu_pheno.txt`

This dataset is intentionally small. It is useful for checking that the workflow runs end to end, but it is not large enough to provide stable real-data inference.

### PLINK genotype input

`fastGxE` reads PLINK binary genotype data through a prefix passed to `--bfile`.

For example:

```bash
--bfile example/test
```

means the program will use:

- `example/test.bed`
- `example/test.bim`
- `example/test.fam`

The same prefix convention is used for `--grm`. For example, `--grm example/test` refers to the GRM prefix, not to a single file such as `example/test.grm.bin`.

### Phenotype and covariate table

The phenotype table passed to `--data` should satisfy the following rules:

1. A header row is required.
2. Fields should be separated by blanks or tabs.
3. The first column must contain individual IDs.
4. Only individuals present in both the phenotype table and the GRM will be retained.
5. Duplicate sample IDs are not allowed.

A minimal layout looks like this:

| iid | pheno | age | sex | E1 | E2 | E3 | ... | E40 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1110005 | -0.75286 | 54 | F | 0.21092 | 0.648552 | 0.468537 | ... | -0.55331 |
| 1110006 | -0.62861 | 49 | M | 0.000148 | -0.12916 | -1.77284 | ... | 0.836599 |
| 1110008 | 1.132841 | 61 | F | -1.14961 | -0.29059 | 0.446468 | ... | -0.10859 |

Typical column usage is:

- `--trait pheno` for the phenotype
- `--covar age` for continuous adjustment covariates
- `--class sex` for categorical covariates
- `--env-int E1:E40` for interacting environments

### Important option behavior

- `--env-int` accepts either explicit names or a header range such as `E1:E40`.
- Variables listed in `--env-int` should not also be listed in `--covar` or `--class`.
- Interacting environments are standardized by default before GxE testing.
- Use `--no-standardize-env` if you want to keep the raw scale.
- Missing phenotype or covariate tokens default to:
  `NA Na na NAN NaN nan -NAN -NaN -nan <NA> <na> N/A n/a`
- You can override or extend the missing-value list with `--missing-data`.

## 2. Build and prepare the GRM

### Step 1: build a dense additive GRM

Generate the additive GRM from the example PLINK files:

```bash
fastgxe --make-grm --bfile example/test --out example/test
```

This produces the standard GRM sidecar files under the chosen prefix, including:

- `example/test.grm.bin`
- `example/test.grm.id`
- `example/test.grm.N`

### Step 2: split and merge for large datasets

For large datasets, the GRM can be computed in partitions with `--npart <total> <part>`:

```bash
fastgxe --make-grm --bfile example/test --npart 4 1 --out example/test
fastgxe --make-grm --bfile example/test --npart 4 2 --out example/test
fastgxe --make-grm --bfile example/test --npart 4 3 --out example/test
fastgxe --make-grm --bfile example/test --npart 4 4 --out example/test
```

Each job writes files such as:

- `example/test.4_1.grm.bin`
- `example/test.4_1.grm.id`
- `example/test.4_1.grm.N`

After all parts are finished, merge them into a new GRM prefix:

```bash
fastgxe --process-grm --merge --grm example/test --npart 4 --out example/test_merged
```

Using an explicit `--out` prefix is recommended for merge jobs.

### Step 3: group related samples for downstream fastGxE analysis

The downstream `--test-main` and `--test-gxe` workflows rely on grouped GRM sidecar files. Create them with:

```bash
fastgxe --process-grm --group --grm example/test --cut-value 0.05
```

This step identifies connected components in the relatedness graph and writes grouped sidecar files, including:

- `example/test.grm.group`
- `example/test.grm.group.size`
- a sparse within-group GRM representation used by downstream analysis

In practice, the GRM preparation workflow is:

1. `--make-grm`
2. `--process-grm --group`
3. `--test-main` or `--test-gxe`

### Optional: reformat a GRM

For inspection or downstream tools, you can reformat a GRM:

```bash
fastgxe --process-grm --reformat --grm example/test --out-fmt 1 --out example/test_triplet
```

This writes an index-triplet representation:

- `example/test_triplet.grm.index_triplet`
- `example/test_triplet.grm.id`

## 3. Test SNP main effects

Run a genome-wide SNP main-effect scan with:

```bash
fastgxe --test-main \
  --grm example/test \
  --bfile example/test \
  --data example/test_simu_pheno.txt \
  --trait pheno \
  --out example/test_main
```

Main outputs are:

- `example/test_main.var`
- `example/test_main.random.res`
- `example/test_main.res`

`test_main.res` contains one row per SNP with the columns:

- `order chrom SNP cm base allele1 allele2 af missing beta se p`

### Useful options

- `--covar` adds continuous covariates.
- `--class` adds categorical covariates.
- `--p-cut` filters the final output by p-value. The default `2` keeps all SNPs.
- `--num-random-snp` controls the number of random SNPs used in calibration.

### Variance components only

If you omit `--bfile`, the command still fits the variance components and writes `*.var`, then stops before the genome-wide SNP scan.

## 4. Test GxE interactions

Run a genome-wide GxE scan with:

```bash
fastgxe --test-gxe \
  --grm example/test \
  --bfile example/test \
  --data example/test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --out example/test_gxe
```

This command treats `E1` through `E40` as the interacting environments.

Main outputs are:

- `example/test_gxe.var`
- `example/test_gxe.main.random.res`
- `example/test_gxe.GxEnoMain.random.res`
- `example/test_gxe.GxE.random.res`
- `example/test_gxe.res`

The final `test_gxe.res` file contains:

- SNP main-effect estimates: `beta_main`, `se_main`, `p_main`
- one `beta`, `se`, `p` triplet for each interacting environment
- combined GxE summary columns: `p_single`, `p_multi`, `p_gxe`

### Important defaults

- Interacting environments are standardized by default.
- Noise-by-environment interaction terms are enabled by default.

To disable environment standardization:

```bash
fastgxe --test-gxe \
  --grm example/test \
  --bfile example/test \
  --data example/test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --no-standardize-env \
  --out example/test_gxe_rawE
```

To disable noise-by-environment terms:

```bash
fastgxe --test-gxe \
  --grm example/test \
  --bfile example/test \
  --data example/test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --no-noisebye \
  --out example/test_gxe_no_noisebye
```

## 5. Parallelize genome-wide scans

Both `--test-main` and `--test-gxe` support task splitting with:

```bash
--split-task <total_parts> <current_part>
```

For example, split a GxE scan into 10 jobs:

```bash
fastgxe --test-gxe \
  --grm example/test \
  --bfile example/test \
  --data example/test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --split-task 10 1 \
  --out example/test_gxe_part1
```

Repeat the same command with part indices `2` through `10`.

Important notes:

- `--split-task` and `--snp-range` cannot be used together.
- Give each parallel job a distinct `--out` prefix.
- The result files are written under task-specific output names after the split is resolved.
- Merging association result files is currently left to an external Python or R script.

## 6. Use fastGxE output with mmSuSiE

Once candidate GxE loci are identified, you can use **mmSuSiE** to prioritize the environmental variables driving the signal:

https://github.com/chaoning/mmsusie

First, export the GRM in index-triplet form:

```bash
fastgxe --process-grm --reformat --grm example/test --out-fmt 1 --out example/test_triplet
```

This creates a GRM representation suitable for downstream Python workflows.

An example `mmSuSiE` workflow is:

```python
from mmsusie import MMSuSiE

model = MMSuSiE()

bed_file = "example/test"
grm_file = "example/test_triplet"
pheno_file = "example/test_simu_pheno.txt"
env_int = [f"E{i}" for i in range(1, 41)]

assoc_file = "example/test_gxe.res"
df_leading = model.ld_pure(
    assoc_file,
    bed_file,
    ld_r2=0.1,
    snp="SNP",
    p="p_gxe",
    p_cutoff=5e-8,
)
print(df_leading)

snp_id = "rs550011"
trait = "pheno"
varcom_file = "example/test_gxe.var"
out_file = "example/test_mmsusie"

res_dct = model.run(
    pheno_file,
    trait,
    env_int,
    grm_file,
    bed_file,
    snp_id,
    varcom_file,
    out_file,
)

print(res_dct["cs"])
df_pip = res_dct["pip"]
print(df_pip)
print(df_pip[df_pip["pip"] > 0.5])
```

## 7. Summary workflow

For a full example run, the main steps are:

1. Build the GRM with `--make-grm`.
2. Group related samples with `--process-grm --group`.
3. Run `--test-main` or `--test-gxe`.
4. Reformat the GRM and follow up with `mmSuSiE` if you want environment-level fine-mapping.
