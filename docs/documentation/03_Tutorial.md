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

The commands use `fastgxe` for the executable. They work with both version 1.0 (`fastgxe`) and version 2.0 (`fastgxe2.0`), except for the three-step workflow in section 5, which needs `fastgxe2.0`. If the binary is not on your `PATH`, use its full path, e.g. `../build/fastgxe2.0`.

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

All result files are tab-separated. The main outputs are:

- `output/test_main.assoc` — genome-wide association results
- `output/test_main.var` — estimated variance components with standard errors
- `output/test_main.h2` — heritability decomposition (`G`, `residual`, `total`)
- `output/test_main.gamma` — random-SNP calibration summary

`test_main.assoc` contains one row per SNP with the columns:

- `order chrom SNP base allele1 allele2 af N beta se p`

Here `af` is the allele frequency and `N` is the number of non-missing genotypes for that SNP.

### Useful options

- `--covar` adds continuous covariates.
- `--class` adds categorical covariates.
- `--num-random-snp` controls the number of random SNPs used for calibration (default `1000`).
- `--log-transform` applies a natural-log transform to the raw phenotype (all values must be positive).
- `--inv-normal` applies a rank-based inverse-normal transform to the phenotype.
- `--keep-random` additionally writes the per-SNP calibration file `output/test_main.random.assoc` (off by default; the scalar calibration factor is always written to `*.gamma`).

### Variance components only

If you omit `--bfile`, the command still fits the variance components and writes `*.var` and `*.h2`, then stops before the genome-wide SNP scan.

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

All result files are tab-separated. The main outputs are:

- `output/test_gxe.assoc` — main results: SNP main effect and combined GxE p-values
- `output/test_gxe.perE.assoc` — single-environment interaction tests (one triplet per environment)
- `output/test_gxe.var` — estimated variance components with standard errors
- `output/test_gxe.h2` — heritability decomposition (`G`, `GxE`, `NxE`, `residual`, `total`)
- `output/test_gxe.gamma` — random-SNP calibration summary

`test_gxe.assoc` contains, per SNP:

- SNP annotation: `order chrom SNP base allele1 allele2 af N`
- SNP main-effect estimates: `beta_main`, `se_main`, `p_main`
- combined GxE summary columns: `p_single`, `p_multi`, `p_gxe`

`test_gxe.perE.assoc` contains the same SNP annotation, followed by a `beta_<env>`, `se_<env>`, `p_<env>` triplet for **each** interacting environment (columns are named after the environment, e.g. `beta_E1 se_E1 p_E1`).

### Useful options

- `--exact-cut` sets the adaptive exact-refit threshold (default `1e-5`): a promising SNP whose fast approximate p-value falls below this is recomputed exactly, so top hits get exact p-values while the genome-wide scan stays fast. Set to `0` for pure approximation or `1` to always refit exactly.
- `--log-transform` / `--inv-normal` transform the phenotype. For `--test-gxe` the inverse-normal transform is applied to the residual **after** covariate/environment correction.
- `--keep-random` also writes the per-SNP calibration files (`*.main.random.assoc`, `*.GxEnoMain.random.assoc`, `*.GxE.random.assoc`); off by default.

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

## 5. Run the GxE scan in three steps

Version 2.0 (`fastgxe2.0`) can run the GxE analysis as three separate steps, similar to [regenie](https://rgcgithub.github.io/regenie/options/). Each step reads the files written by the step before. The null model is fitted once, and every scan job uses the same gamma factors, however the genome is divided.

| Step | What it does | Inputs | Outputs |
| --- | --- | --- | --- |
| `--step 1` | Fits the null model and estimates the variance components | GRM, phenotype table | `.var`, `.h2`, `.null.bin` |
| `--step 2` | Draws random SNPs from the whole genome and calibrates gamma | step 1, genotypes | `.gamma` |
| `--step 3` | Scans the SNPs for GxE | step 1, step 2, genotypes | `.assoc`, `.perE.assoc` |

The output files have the same format as those of `--test-gxe`. Use `fastgxe2.0 --step N -h` to list the options of a step. An option that does not belong to a step is rejected.

### Step 1: fit the null model

```bash
fastgxe2.0 --step 1 \
  --grm output/test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --env-int E1:E40 \
  --out output/s1
```

This step needs no genotypes. It takes the same phenotype, covariate and environment options as `--test-gxe`, such as `--covar`, `--class`, `--inv-normal` and `--no-noisebye`.

`output/s1.null.bin` holds what the later steps need, so they take no phenotype options. It contains sample IDs and environment values, so treat it as individual-level data.

### Step 2: calibrate gamma

With one fileset for the whole genome:

```bash
fastgxe2.0 --step 2 --null output/s1 --bfile test --seed 2026 --out output/s2
```

With one fileset per chromosome, list the PLINK prefixes in a text file, one per line, and pass it to `--mbfile`:

```bash
for c in {1..22}; do echo "geno/chr$c"; done > chr_list.txt
fastgxe2.0 --step 2 --null output/s1 --mbfile chr_list.txt --seed 2026 --out output/s2
```

The SNPs are drawn evenly across all filesets, so each chromosome contributes in proportion to its number of SNPs. For the same seed, the same SNPs are drawn whether the genome is stored in one fileset or several, provided the SNPs are listed in the same order.

- `--num-random-snp` (default `1000`) counts SNPs that pass `--maf` and `--missing-rate`. SNPs that fail are replaced by further random SNPs.
- `--seed` makes the draw reproducible. Without it a random seed is chosen. The seed is always recorded in `output/s2.gamma`.
- `--keep-random` also writes the per-SNP calibration files.

### Step 3: scan the genome

Whole genome in one job:

```bash
fastgxe2.0 --step 3 --null output/s1 --gamma output/s2 --bfile test --out output/s3
```

One job per chromosome:

```bash
for c in {1..22}; do
  fastgxe2.0 --step 3 --null output/s1 --gamma output/s2 --bfile geno/chr$c --out output/s3.chr$c
done
```

One chromosome split into 20 jobs, which write `output/s3.chr1.20_1.assoc` to `output/s3.chr1.20_20.assoc`:

```bash
for i in {1..20}; do
  fastgxe2.0 --step 3 --null output/s1 --gamma output/s2 --bfile geno/chr1 --split-task 20 $i --out output/s3.chr1
done
```

The whole genome split into 200 jobs of equal size. A job may span two chromosomes:

```bash
for i in {1..200}; do
  fastgxe2.0 --step 3 --null output/s1 --gamma output/s2 --mbfile chr_list.txt --split-task 200 $i --out output/s3
done
```

Each loop body can be submitted as a job array. Merge the results, keeping one header line:

```bash
awk 'FNR>1 || NR==1' output/s3.chr{1..22}.assoc > output/s3.assoc
awk 'FNR>1 || NR==1' output/s3.chr{1..22}.perE.assoc > output/s3.perE.assoc
```

Useful options:

- `--bsize` (default `1000`) sets how many SNPs are read at a time. Genotypes take about `samples * bsize * 8` bytes of memory.
- `--exact-cut`, `--p-approx-cut`, `--maf`, `--missing-rate` and `--snp-range` work as for `--test-gxe`. `--snp-range` needs `--bfile`.
- `--exact-cut 0` skips the exact refit and does not load `V^-1`, which saves memory, but the top hits then only have approximate p-values.

Notes:

- The `order` column counts SNPs within the input of that run. With `--bfile` it restarts in each file; with `--mbfile` it counts across all filesets.
- `.gamma` records which null model it belongs to. Step 3 stops if `--null` and `--gamma` do not match.
- Splitting the work changes results only through rounding. `p_multi` and `p_gxe` can differ by about `1e-6` (relative) for SNPs near `p_multi = 0.46`, where the saddlepoint approximation is sensitive to the last digits of the score. Significant SNPs are not affected.

## 6. Parallelize single-command scans

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
- Each job fits the null model again and draws its own random SNPs, so the gamma factors differ slightly between jobs. Use the three-step workflow in section 5 if all jobs should share them.

## 7. Identify environments driving lead SNP GxE signals with mmSuSiE

Once candidate GxE loci are identified, you can use **mmSuSiE** to prioritize the environmental variables driving the signal:

https://github.com/chaoning/mmsusie

### Installation

```bash
git clone https://github.com/chaoning/mmsusie.git
cd mmsusie
pip install .
```

### Prepare the sparse GRM

mmSuSiE's sparse workflow reads two files under one prefix: `output/test.grm.group`
(sample IDs + relatedness groups) and `output/test.grm.index_triplet` (the GRM
values). The `.grm.group` file is written by the grouping step; if you have not run
it yet, do so first:

```bash
fastgxe --process-grm --group --grm output/test --cut-value 0.05
```

Then reformat the GRM to index-triplet form:

```bash
fastgxe --process-grm --reformat --sparse --grm output/test --out-fmt 1 --out output/test
```

`--reformat` writes `output/test.grm.index_triplet` (and `output/test.grm.id`); keep
`--out` on the same prefix as `--grm` so both files sit under one prefix.
`MMSuSiESp.read_sp_grm` then reads `.grm.group` and `.grm.index_triplet` — the
`.grm.id` file is produced but not required.

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
assoc_file = "output/test_gxe.assoc"
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
varcom_file = "output/test_gxe.var"   # variance components from --test-gxe, or the .var of --step 1
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
    estimate_sigma=False, # keep the fastGxE variance components fixed (recommended)
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

> **`estimate_sigma=False` is recommended here.** It reuses the variance components
> that `--test-gxe` estimated **genome-wide** (over all SNPs) — the best-calibrated
> values for this locus. `estimate_sigma=True` instead re-estimates the components in
> the loop from this single locus; that refit is now REML (it corrects for the many
> fixed effects mmSuSiE projects out — intercept, covariates, all environment main
> effects, and the lead-SNP main effect, ≈ `2 + n_env` columns), so it is stable and no
> longer biases the components downward. Prefer `False` for the genome-wide calibration;
> use `True` only if you specifically want the components re-estimated jointly with the
> fine-mapping.

### Output files

`model.mmsusie_lead_gxe()` writes four result files under the given `out_file` prefix:

| File | Contents |
| --- | --- |
| `output/test_mmsusie.pip.txt` | Posterior inclusion probability for each environment |
| `output/test_mmsusie.alpha.txt` | Per-component posterior assignment probabilities |
| `output/test_mmsusie.mu.txt` | Posterior mean effects |
| `output/test_mmsusie.cs.txt` | Credible sets |

## 8. Summary workflow

For a full example run, the main steps are:

1. `cd example && mkdir -p output`
2. Build the GRM with `--make-grm`.
3. Group related samples with `--process-grm --group`.
4. Run `--test-main`, or run the GxE scan with `--test-gxe` or with `--step 1`, `--step 2` and `--step 3`.
5. Reformat the GRM and follow up with `mmSuSiE` if you want environment-level fine-mapping.
