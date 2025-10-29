import pandas as pd
import numpy as np
import sys
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from acat import acat
from saddle_point import pchisqsum
from scipy.stats import chi2

from sklearn.linear_model import LinearRegression

num_envi_power = int(sys.argv[1])
power_snp_h2_gxe = float(sys.argv[2])/100
sample_size = int(sys.argv[3])

dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi_new_2M.csv", header=0)
dfE = dfE.sample(n=sample_size, random_state=42)
dfE = dfE.reset_index(drop=True)
dataE = np.array(dfE)
corrE = dfE.corr().values
eigenvalues = np.linalg.eigvals(corrE)


def simulate_genotype_hwe(N=1000, p=0.3, seed=42):
    """
    Simulate genotypes (encoded as 0, 1, 2) under Hardy-Weinberg Equilibrium (HWE)
    
    Parameters:
    - N: Number of individuals
    - p: Allele frequency of the A allele (0 <= p <= 1)
    - seed: Random seed for reproducibility
    
    Returns:
    - A NumPy array of length N containing genotypes for each individual
      where:
        0 = homozygous for a (aa)
        1 = heterozygous (Aa)
        2 = homozygous for A (AA)
    """
    np.random.seed(seed)
    q = 1 - p  # Frequency of allele a
    genotype_probs = [q**2, 2*p*q, p**2]  # Probabilities for genotypes: aa, Aa, AA
    genotypes = np.random.choice([0, 1, 2], size=N, p=genotype_probs)
    genotypes = genotypes - np.mean(genotypes)  # Centering genotypes
    genotypes = genotypes / np.std(genotypes)  # Standardizing genotypes
    return genotypes

def simulate_phenotype(genotypes, environmental_factors, power_snp_h2_gxe, num_envi_power, seed=42):

    # GxE effects
    np.random.seed(seed)
    num_envi = environmental_factors.shape[1]
    power_envi_indices = np.random.choice(range(num_envi), num_envi_power, replace=False)
    power_envi_indices = np.sort(power_envi_indices)
    dataE_power = environmental_factors[:, power_envi_indices]
    eff = np.random.normal(size=num_envi_power)
    genotypes = genotypes.reshape(-1, 1)  # Ensure genotypes is a 2D array

    eff_var = np.sum(np.var(genotypes * dataE_power * eff.reshape(1, -1), axis=0))
    eff *= np.sqrt(power_snp_h2_gxe / eff_var)
    GxE_power_sum = np.dot(genotypes * dataE_power, eff)

    error = np.random.normal(size=environmental_factors.shape[0]) * np.sqrt(1 - power_snp_h2_gxe)
    y = GxE_power_sum + error

    return y

res_lst = []
for i in range(100):
    freq = np.random.uniform(0.05, 0.5)
    genotypes = simulate_genotype_hwe(N=sample_size, p=freq, seed=i)
    y = simulate_phenotype(genotypes, dataE, power_snp_h2_gxe, num_envi_power, seed=i)

    # GxE test for each environmental factor, linear regression
    p_lst = []
    z_lst = []
    for j in range(dataE.shape[1]):
        X = dataE[:, j] * genotypes
        effect = np.dot(X, y) / np.sum(X**2)
        se = 1 / np.sqrt(np.sum(X**2))
        t_stat = effect / se
        p_value = chi2.sf(t_stat ** 2, df=1)
        p_lst.append(p_value)
        z_lst.append(t_stat)
    _, p_acat = acat(p_lst)
    chi2_stat = np.sum(np.array(z_lst) ** 2)
    # saddle point approximation
    p_saddle = pchisqsum(np.array([chi2_stat]), df=[1] * len(z_lst), a=eigenvalues, lower_tail=True, method="saddlepoint")[0]
    _, p_gxe = acat([p_acat, p_saddle])
    res_lst.append([i, freq, p_acat, p_saddle, p_gxe])
    print(i)

dfR = pd.DataFrame(res_lst, columns=["rep", "freq", "p_acat", "p_saddle", "p_gxe"])
dfR.to_csv(f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_take_home_extend/{sample_size}.h2_{sys.argv[2]}.E_{num_envi_power}.csv", index=False)
