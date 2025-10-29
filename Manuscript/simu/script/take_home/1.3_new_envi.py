import numpy as np
import pandas as pd

# Step 1: Read the environmental variable file and compute the correlation matrix
# 'eid' is read as string to avoid parsing issues
dfE = pd.read_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", dtype={"eid": str})
env_data = dfE.loc[:, "age":"alcohol_frequency"]  # Select environmental variables
corr_mat = env_data.corr()  # Compute the correlation matrix

# Step 2: Set random seed for reproducibility
np.random.seed(42)

# Step 3: Define the number of samples and variables
n_samples = 2_000_000  # 2 million individuals
n_variables = corr_mat.shape[0]  # Number of environmental variables
mean_vec = np.zeros(n_variables)  # Zero mean for each variable

# Step 4: Sample from a multivariate normal distribution
# Use the correlation matrix as the covariance matrix since mean is 0
synthetic_env = np.random.multivariate_normal(mean=mean_vec, cov=corr_mat.values, size=n_samples)

# Step 5: Convert the result into a DataFrame
synthetic_env_df = pd.DataFrame(synthetic_env, columns=corr_mat.columns)

# Optional: Save the synthetic data to CSV
synthetic_env_df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi_new_2M.csv", index=False)

