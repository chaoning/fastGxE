---
layout: page
title: Tutorial
description:
order: 3
---

This tutorial walks through the example data bundled with `fastGxE` and follows the current command-line workflow implemented in the codebase.

All commands below assume you have changed into the `example/` directory and created an `output/` subdirectory for results:

```bash
cd example
mkdir -p output
```

The executable is assumed to be on your `PATH` as `fastgxe`. If you built from source and the binary is not on your `PATH`, replace `fastgxe` with the full path, e.g. `../build/fastgxe`.

## 1. Example inputs

The `example/` directory contains a minimal dataset for testing the software:

- `test.bed`
- `test.bim`
- `test.fam`
- `test_simu_pheno.txt`

This dataset is intentionally small. It is useful for checking that the workflow runs end to end, but it is not large enough to provide stable real-data inference.

### PLINK genotype input

`fastGxE` reads PLINK binary genotype data through a prefix passed to `--bfile`.

For example:

```bash
--bfile test
```

means the program will use:

- `test.bed`
- `test.bim`
- `test.fam`

The same prefix convention is used for `--grm`. For example, `--grm test` refers to the GRM prefix, not to a single file such as `test.grm.bin`.

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
fastgxe --make-grm --bfile test --out output/test
```

This produces the standard GRM sidecar files under the chosen prefix, including:

- `output/test.grm.bin`
- `output/test.grm.id`
- `output/test.grm.N`

For large datasets, the GRM can be computed in partitions with `--npart <total> <part>` and then merged:

```bash
fastgxe --make-grm --bfile test --npart 4 1 --out output/test
fastgxe --make-grm --bfile test --npart 4 2 --out output/test
fastgxe --make-grm --bfile test --npart 4 3 --out output/test
fastgxe --make-grm --bfile test --npart 4 4 --out output/test
fastgxe --process-grm --merge --grm output/test --npart 4
```

Each partition writes files such as `output/test.4_1.grm.bin`, `output/test.4_1.grm.id`, and `output/test.4_1.grm.N`. The merge step combines them back into `output/test.grm.bin`.

### Step 2: group related samples for downstream fastGxE analysis

The downstream `--test-main` and `--test-gxe` workflows rely on grouped GRM sidecar files. Create them with:

```bash
fastgxe --process-grm --group --grm output/test --cut-value 0.05
```

This step identifies connected components in the relatedness graph and writes grouped sidecar files, including:

- `output/test.grm.group`
- `output/test.grm.group.size`
- a sparse within-group GRM representation used by downstream analysis

## 3. Test SNP main effects

Run a genome-wide SNP main-effect scan with:

```bash
fastgxe --test-main \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --out output/test_main
```

Main outputs are:

- `output/test_main.var`
- `output/test_main.random.res`
- `output/test_main.res`

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
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --out output/test_gxe
```

This command treats `E1` through `E40` as the interacting environments.

Main outputs are:

- `output/test_gxe.var`
- `output/test_gxe.main.random.res`
- `output/test_gxe.GxEnoMain.random.res`
- `output/test_gxe.GxE.random.res`
- `output/test_gxe.res`

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
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --no-standardize-env \
  --out output/test_gxe_rawE
```

To disable noise-by-environment terms:

```bash
fastgxe --test-gxe \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --no-noisebye \
  --out output/test_gxe_no_noisebye
```

## 5. Parallelize genome-wide scans

Both `--test-main` and `--test-gxe` support task splitting with:

```bash
--split-task <total_parts> <current_part>
```

For example, split a GxE scan into 10 jobs:

```bash
fastgxe --test-gxe \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --split-task 10 1 \
  --out output/test_gxe_part1
```

Repeat the same command with part indices `2` through `10`.

Important notes:

- `--split-task` and `--snp-range` cannot be used together.
- Give each parallel job a distinct `--out` prefix.
- The result files are written under task-specific output names after the split is resolved.
- Merging association result files is currently left to an external Python or R script.

## 6. Identify environments driving lead SNP GxE signals with mmSuSiE

Once candidate GxE loci are identified, you can use **mmSuSiE** to prioritize the environmental variables driving the signal:

https://github.com/chaoning/mmsusie

### Installation

```bash
git clone https://github.com/chaoning/mmsusie.git
cd mmsusie
pip install .
```

### Prepare the sparse GRM

Export the GRM in index-triplet form required by mmSuSiE:

```bash
fastgxe --process-grm --reformat --grm output/test --out-fmt 1 --out output/test
```

This writes the files used by the sparse block-diagonal workflow:

- `output/test.grm.id`
- `output/test.grm.group`
- `output/test.grm.index_triplet`

### Run mmSuSiE

The sparse-GRM workflow uses the `MMSuSiESp` class:

```python
from mmsusie import MMSuSiESp

model = MMSuSiESp()

# Input files (run from the example/ directory)
bed_file = "test"                           # PLINK genotype prefix
grm_file = "output/test"                    # sparse GRM prefix from --process-grm --reformat
pheno_file = "test_simu_pheno.txt"          # phenotype/covariate table
env_int = [f"E{i}" for i in range(1, 41)]  # environment columns to fine-map

# Step 1: LD-prune significant GxE hits to obtain one lead SNP per locus.
# ld_pure keeps the most significant SNP in each LD block (r² < ld_r2).
assoc_file = "output/test_gxe.res"
df_leading = model.ld_pure(
    assoc_file,
    bed_file,
    ld_r2=0.1,       # prune SNPs with pairwise r² >= 0.1
    snp="SNP",       # column name for SNP ID in assoc_file
    p="p_gxe",       # column name for GxE p-value
    p_cutoff=5e-8,   # genome-wide significance threshold
)
print(df_leading)

# Step 2: Fine-map environments for a chosen lead SNP.
snp_id = "rs550011"                    # lead SNP from df_leading
trait = "pheno"                        # phenotype column in pheno_file
varcom_file = "output/test_gxe.var"   # variance components estimated by fastGxE
out_file = "output/test_mmsusie"      # output file prefix

res_dct = model.mmsusie_lead_gxe(
    pheno_file,
    trait,
    env_int,
    grm_file,
    bed_file,
    snp_id,
    varcom_file,
    out_file,
    L=10,                # maximum number of causal environments
    maxiter=100,         # maximum ELBO iterations
    tol=1e-3,            # convergence tolerance on ELBO
    coverage=0.95,       # credible set coverage level
    min_abs_corr=0.5,    # minimum purity for a credible set to be reported
    estimate_sigma=False, # use fixed prior variance during fitting
)

# res_dct["cs"]             — credible sets: environments with cumulative PIP >= coverage
# res_dct["pip"]            — per-environment posterior inclusion probabilities (DataFrame)
# res_dct["alpha"]          — L x p posterior assignment probabilities per component
# res_dct["mu"]             — L x p posterior mean effects per component
# res_dct["elbo"]           — ELBO values across iterations (convergence diagnostic)
print(res_dct["cs"])
df_pip = res_dct["pip"]
print(df_pip)
print(df_pip[df_pip["pip"] > 0.5])  # environments with strong evidence of interaction
```

### Output files

`model.mmsusie_lead_gxe()` writes four result files under the given `out_file` prefix:

| File | Contents |
| --- | --- |
| `output/test_mmsusie.pip.txt` | Posterior inclusion probability for each environment |
| `output/test_mmsusie.alpha.txt` | Per-component posterior assignment probabilities |
| `output/test_mmsusie.mu.txt` | Posterior mean effects |
| `output/test_mmsusie.cs.txt` | Credible sets |

## 7. Summary workflow

For a full example run, the main steps are:

1. `cd example && mkdir -p output`
2. Build the GRM with `--make-grm`.
3. Group related samples with `--process-grm --group`.
4. Run `--test-main` or `--test-gxe`.
5. Reformat the GRM and follow up with `mmSuSiE` if you want environment-level fine-mapping.
