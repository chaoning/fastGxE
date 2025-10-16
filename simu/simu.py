import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
import os

# --------------------------
# utilities
# --------------------------
def is_pos_def(A):
    """Return True if matrix A is positive definite (Cholesky succeeds)."""
    try:
        np.linalg.cholesky(A)
        return True
    except np.linalg.LinAlgError:
        return False

def nearPD_eig(A, eps=1e-8):
    """
    Eigen-decomposition based near positive-definite approximation:
    replace eigenvalues < eps with eps and reconstruct the matrix.
    """
    A = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals_clipped = np.clip(eigvals, eps, None)
    A_pd = (eigvecs * eigvals_clipped) @ eigvecs.T
    return (A_pd + A_pd.T) / 2.0

def make_corr_from_sym(A, ensure_pd=True, eps=1e-8, jitter=1e-10):
    """
    Project a symmetric matrix A to a valid correlation matrix:
      - symmetric
      - diagonal = 1
      - positive semi-definite (attempted)
    If ensure_pd is True, try nearPD and add tiny jitter if needed.
    """
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 1.0)

    if ensure_pd and not is_pos_def(A):
        A = nearPD_eig(A, eps=eps)
        # Normalize to correlation (make diagonals 1)
        diag = np.sqrt(np.diag(A))
        diag[diag < 1e-12] = 1e-12
        A = A / np.outer(diag, diag)
        A = (A + A.T) / 2.0
        np.fill_diagonal(A, 1.0)

        # If still not PD, try adding small jitter iteratively
        if not is_pos_def(A):
            tries = 0
            while tries < 5 and not is_pos_def(A):
                tries += 1
                A = A + np.eye(A.shape[0]) * jitter * (10 ** tries)
            if not is_pos_def(A):
                raise np.linalg.LinAlgError("Cannot make matrix positive definite even after jitter.")
    else:
        np.fill_diagonal(A, 1.0)

    # Clip numerical values into [-1, 1], symmetrize, and set diagonal to 1
    A = np.clip(A, -1.0, 1.0)
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 1.0)
    return A

def process_snp_data(snp_indices, snp_on_disk, fillna_with_mean=True):
    """
    Read specified SNP columns from a Bed object and return standardized numpy array.
    snp_indices: integer or list of integers (global SNP indices)
    Returns shape (n_individuals, n_snps)
    """
    if np.isscalar(snp_indices):
        snp_indices = [int(snp_indices)]
    snp_indices = list(map(int, snp_indices))

    arr = snp_on_disk[:, snp_indices].read().val  # shape (n_indiv, n_snps)
    arr = np.array(arr, dtype=float)

    # Fill missing values by column mean if requested
    if fillna_with_mean:
        col_means = np.nanmean(arr, axis=0)
        inds_nan = np.where(np.isnan(arr))
        if inds_nan[0].size > 0:
            arr[inds_nan] = np.take(col_means, inds_nan[1])

    # Standardize columns (avoid division by zero)
    col_mean = np.mean(arr, axis=0)
    col_std = np.std(arr, axis=0, ddof=0)
    col_std[col_std < 1e-12] = 1.0
    arr = (arr - col_mean) / col_std
    return arr

# --------------------------
# main
# --------------------------
os.chdir("/net/zootopia/disk1/chaon/WORK/fastgxe/example/")

# reproducibility
seed = 2025
rng = np.random.default_rng(seed)

num_envi = 40
bed_file = "test"
iid_lst = pd.read_csv(bed_file + ".fam", sep=r"\s+", header=None, dtype=str)[1].tolist()
sid_lst = pd.read_csv(bed_file + ".bim", sep=r"\s+", header=None, dtype=str)[1].tolist()
snp_on_disk = Bed(bed_file, count_A1=True)
num_iid = snp_on_disk.iid_count
num_snp = snp_on_disk.sid_count

# ---- environmental covariance (correlation) matrix ----
env_cov = np.eye(num_envi)
iu = np.triu_indices(num_envi, k=1)
# generate uniform random off-diagonal values between -0.2 and 0.2
offdiag = rng.uniform(-0.2, 0.2, size=iu[0].shape[0])
env_cov[iu] = offdiag
env_cov = env_cov + env_cov.T - np.diag(np.diag(env_cov))

print("Initial env_cov symmetric:", np.allclose(env_cov, env_cov.T))
print("Initial diag (should be 1):", np.diag(env_cov)[:5])

# If not PD, adjust and ensure it's a correlation matrix
env_cov = make_corr_from_sym(env_cov, ensure_pd=True, eps=1e-8, jitter=1e-12)
print("Post-processed env_cov is PD:", is_pos_def(env_cov))
# print minimum eigenvalue for checking
eigvals = np.linalg.eigvalsh(env_cov)
print("min eigenval env_cov:", eigvals.min())

# ----------------- simulate environmental data -----------------
env_data = rng.multivariate_normal(mean=np.zeros(num_envi), cov=env_cov, size=num_iid)
env_data = pd.DataFrame(env_data, columns=[f"E{i+1}" for i in range(num_envi)])
# standardize columns
env_data = (env_data - env_data.mean()) / env_data.std(ddof=0)

# ----------------- genetic/environment parameters -----------------
h2_env = 0.01
h2_g = 0.3
h2_gxe = 0.5
h2_gxe_focus = 0.4
num_envi_focus = 2
num_snp_g = 1000
num_snp_gxe = 500

# sample SNP indices reproducibly
snp_g_indices = rng.choice(num_snp, size=num_snp_g, replace=False)
snp_gxe_indices = rng.choice(num_snp, size=num_snp_gxe, replace=False)
snp_gxe_focus_index = rng.choice(num_snp, size=1, replace=True)

# environmental main effects
env_effect = rng.normal(0, np.sqrt(h2_env / num_envi), size=num_envi)
envs_effect = env_data.values @ env_effect  # shape (num_iid,)
envs_effect = envs_effect / np.std(envs_effect) * np.sqrt(h2_env)  # rescale to exact variance h2_env

# SNP main effects (using the selected snp_g_indices)
snp_g_mat = process_snp_data(snp_g_indices, snp_on_disk)  # (num_iid, num_snp_g)
snp_g_effect = rng.normal(0, np.sqrt(h2_g / num_snp_g), size=num_snp_g)
g_effect = snp_g_mat @ snp_g_effect  # (num_iid,)
g_effect = g_effect / np.std(g_effect) * np.sqrt(h2_g)  # rescale to exact variance h2_g

# --------------- construct GxE pairs ---------------
# weights: first 10 environments are more likely to be chosen
numbers = list(range(1, num_envi+1))
weights = np.array([10]*10 + [1]*30, dtype=float)

lst1 = []  # global SNP indices for pairs
lst2 = []  # environment indices aligned with lst1
for isnp in range(num_snp_gxe):
    # sample how many environments to pair with this SNP (1..num_envi)
    chosen_number = int(rng.choice(numbers, p=weights/weights.sum()))
    env_rnd_index = rng.choice(range(num_envi), chosen_number, replace=False)
    env_rnd_index = np.sort(env_rnd_index)
    lst1.extend([int(snp_gxe_indices[isnp])] * len(env_rnd_index))
    lst2.extend(list(env_rnd_index))

total_pairs = len(lst1)
print("Total GxE pairs:", total_pairs)

gxe_effect = np.zeros(num_iid)
if total_pairs > 0:
    snp_gxe_mat = process_snp_data(lst1, snp_on_disk)  # (num_iid, total_pairs)
    env_sub_mat = env_data.iloc[:, lst2].values       # (num_iid, total_pairs)
    # scale per-pair variance so that total GxE variance = h2_gxe - h2_gxe_focus
    var_per_pair = (h2_gxe - h2_gxe_focus) / total_pairs
    if var_per_pair < 0:
        raise ValueError("h2_gxe must be >= h2_gxe_focus")
    snp_gxe_effect = rng.normal(0, np.sqrt(var_per_pair), size=total_pairs)
    # elementwise multiply and sum across pairs to get individual's GxE contribution
    gxe_effect = np.dot(snp_gxe_mat * env_sub_mat, snp_gxe_effect)
    gxe_effect = gxe_effect / np.std(gxe_effect) * np.sqrt(h2_gxe)  # rescale to exact variance h2_gxe

# --------------- focused GxE on one SNP across multiple environments ---------------
print("Focused GxE SNP:", snp_gxe_focus_index[0], sid_lst[snp_gxe_focus_index[0]])
snp_focus_arr = process_snp_data(snp_gxe_focus_index, snp_on_disk)  # shape (num_iid, 1)
env_focus_indices = np.sort(rng.choice(range(num_envi), num_envi_focus, replace=False))
print("Focused GxE env indices:", env_focus_indices, [f"E{i+1}" for i in env_focus_indices])
snp_gxe_focus_effect = rng.normal(0, np.sqrt(h2_gxe_focus / num_envi_focus), size=num_envi_focus)
# broadcast multiply then dot to get focused GxE contribution
gxe_focus_effect = np.dot(snp_focus_arr * env_data.iloc[:, env_focus_indices].values, snp_gxe_focus_effect)
gxe_focus_effect = gxe_focus_effect / np.std(gxe_focus_effect) * np.sqrt(h2_gxe_focus)  # rescale to exact variance h2_gxe_focus

# ----------------- residual and phenotype -----------------
resid_var = 1.0 - h2_env - h2_g - h2_gxe
if resid_var < 0:
    raise ValueError("Total heritabilities > 1. Reduce h2_* params.")
residual = rng.normal(0, np.sqrt(resid_var), size=num_iid)
residual = residual / np.std(residual) * np.sqrt(resid_var)  # rescale to exact variance resid_var

pheno = envs_effect + g_effect + gxe_effect + gxe_focus_effect + residual
print("Phenotype variance components:", np.var(pheno), np.var(envs_effect), np.var(g_effect), np.var(gxe_effect), np.var(gxe_focus_effect), np.var(residual))
print(np.sum([np.var(envs_effect), np.var(g_effect), np.var(gxe_effect), np.var(gxe_focus_effect), np.var(residual)]))

# ----------------- assemble dataframe and save -----------------
dfP = pd.DataFrame({
    "iid": iid_lst,
    "pheno": pheno
})
dfP = pd.concat([dfP.reset_index(drop=True), env_data.reset_index(drop=True)], axis=1)

outfile = "test_simu_pheno.txt"
dfP.to_csv(outfile, sep=" ", index=False)
print("Saved simulated phenotype to:", outfile)

# quick diagnostics
print("pheno mean/std:", pheno.mean(), np.std(pheno, ddof=0))
print("env_data sample means (first 5):", env_data.iloc[:5, :5].mean().values)
print("env_cov diag (first 5):", np.diag(env_cov)[:5])
print("env_cov shape:", env_cov.shape)
