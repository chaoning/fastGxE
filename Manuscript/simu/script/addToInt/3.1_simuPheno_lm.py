import pandas as pd
import numpy as np
import sys
import os
from scipy.stats import pearsonr
import statsmodels.api as sm
from tqdm import tqdm


# Set working directory
os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_addToint_lm/")

# Read command line arguments
unknown_h2SNP_h2E = float(sys.argv[1]) / 100
R2 = float(sys.argv[2])  # squared correlation coefficient between G and unknown G, E and unknown E
sample_size = 100_000

res_lst = []
for rep in tqdm(range(1, 101)):
    np.random.seed(rep)
    # ---------------------------------------------------------
    # Construct covariance matrix for bivariate normal sampling
    # Each pair (true, proxy) is drawn from N(0, Σ) where
    # Var(true) = Var(proxy) = 1, Corr(true, proxy) = sqrt(R2)
    # ---------------------------------------------------------
    rho = np.sqrt(R2)
    cov_matrix = [[1, rho],
                  [rho, 1]]

    # ---------------------------------------------------------
    # Generate genetic values (G_true, G_proxy)
    # ---------------------------------------------------------
    G_joint = np.random.multivariate_normal(
        mean=[0, 0],      # Mean vector
        cov=cov_matrix,   # Covariance matrix
        size=sample_size  # Number of samples
    )
    G_true = G_joint[:, 0]
    G_proxy = G_joint[:, 1]
    R2_G = pearsonr(G_true, G_proxy)[0] ** 2

    # ---------------------------------------------------------
    # Generate environmental values (E_true, E_proxy)
    # ---------------------------------------------------------
    E_joint = np.random.multivariate_normal(
        mean=[0, 0],
        cov=cov_matrix,
        size=sample_size
    )
    E_true = E_joint[:, 0]
    E_proxy = E_joint[:, 1]
    R2_E = pearsonr(E_true, E_proxy)[0] ** 2

    y = np.sqrt(unknown_h2SNP_h2E) * G_true + np.sqrt(unknown_h2SNP_h2E) * E_true + np.random.normal(0, 1, sample_size) * np.sqrt(1 - 2 * unknown_h2SNP_h2E)

    X_true = np.concatenate([G_true.reshape(-1, 1), E_true.reshape(-1, 1), (G_true * E_true).reshape(-1, 1)], axis=1)
    X_proxy = np.concatenate([G_proxy.reshape(-1, 1), E_proxy.reshape(-1, 1), (G_proxy * E_proxy).reshape(-1, 1)], axis=1)

    # Fit linear models and calculate p-values of GxE term
    model = sm.OLS(y, X_true).fit()
    p_value_G_true = model.pvalues[0]
    p_value_trueE = model.pvalues[1]
    p_value_true = model.pvalues[2]

    model = sm.OLS(y, X_proxy).fit()
    p_value_G_proxy = model.pvalues[0]
    p_value_proxyE = model.pvalues[1]
    p_value_proxy = model.pvalues[2]

    GE_true = (G_true * E_true).reshape(-1, 1)
    model = sm.OLS(y, GE_true).fit()
    p_value_true_noMain = model.pvalues[0]
    
    GE_proxy = (G_proxy * E_proxy).reshape(-1, 1)
    model = sm.OLS(y, GE_proxy).fit()
    p_value_proxy_noMain = model.pvalues[0]

    res_lst.append((R2_G, R2_E, p_value_G_true, p_value_G_proxy, p_value_trueE, p_value_proxyE, p_value_true, p_value_proxy, p_value_true_noMain, p_value_proxy_noMain))

p_df = pd.DataFrame(res_lst, columns=['R2_G', 'R2_E', 'p_value_G_true', 'p_value_G_proxy', 'p_value_trueE', 'p_value_proxyE', 'p_value_true', 'p_value_proxy', 'p_value_true_noMain', 'p_value_proxy_noMain'])
p_df.to_csv(f"h2_{sys.argv[1]}_R2_{sys.argv[2]}.csv", index=False)
