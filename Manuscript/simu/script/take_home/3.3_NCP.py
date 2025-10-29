import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ncx2

def calculate_power(ncp, df, alpha=5e-8):
    """
    Calculate statistical power using a non-central chi-square distribution.

    Parameters:
    ncp : float
        Non-centrality parameter λ = n * h²
    df : int
        Degrees of freedom (e.g., number of GxE interaction terms)
    alpha : float
        Significance level (default 5e-8 for GWAS)

    Returns:
    power : float
        Probability of detecting the signal (i.e., power)
    """
    # Get critical chi-square value under H0
    critical_value = ncx2.isf(alpha, df, 0)
    # Compute power under H1
    power = ncx2.sf(critical_value, df, ncp)
    return power

# Parameters
df = 40  # degrees of freedom (e.g., 40 GxE terms)
alpha = 5e-8  # significance level
sample_sizes = [50000, 100000, 200000, 300000, 400000, 800000, 1200000, 1600000, 2000000]  # sample sizes

# SNP GxE heritability values: 0.02% to 0.14% (i.e., 0.0002 to 0.0014)
h2_gxe_values = np.linspace(0.02, 0.14, 7)

h2_gxe_values = [0.005, 0.01] + list(np.linspace(0.02, 0.14, 7))  # Adding 0.005 and 0.01 for more granularity

# Calculate power matrix: one line per sample size
sample_size_lst = []
h2_gxe_lst = []
power_lst = []
for n in sample_sizes:
    powers = [calculate_power(ncp=n * h2 / 100, df=df, alpha=alpha) for h2 in h2_gxe_values]
    sample_size_lst.extend([n] * len(h2_gxe_values))
    h2_gxe_lst.extend(h2_gxe_values)
    power_lst.extend(powers)

df = pd.DataFrame({
    'sample_size': sample_size_lst,
    'power_snp_h2_gxe': h2_gxe_lst,
    'power_val': power_lst
})

df.to_csv('/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_power/take_home_ncp.csv', index=False)
