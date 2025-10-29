import sys
import pandas as pd
import numpy as np

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file (space/tab separated with header)
df = pd.read_csv(input_file, sep=r"\s+", header=0)

# Select the second column as the phenotype/trait column
trait_col = df.columns[1]
y = df[trait_col].to_numpy()

# If there are negative values, stop the script
# because log transform is undefined for negative numbers
if np.nanmin(y) < 0:
    sys.exit(
        f"Error: Minimum value of {trait_col} is < 0. "
        f"Log transform is undefined. Consider Yeo-Johnson "
        f"or shifting the data by a larger positive constant."
    )

# If zeros are present, replace them with a small positive constant (epsilon)
# to avoid log(0). We use half of the minimum positive value as epsilon.
positive_mask = (y > 0)
if np.any(y == 0):
    if np.any(positive_mask):
        min_pos = np.nanmin(y[positive_mask])
        eps = min_pos / 2.0
        # Fallback in case eps is invalid (zero, negative, or NaN)
        if eps <= 0 or not np.isfinite(eps):
            eps = 1e-8
        y = np.where(y == 0, eps, y)
    else:
        # If the column only contains zeros (and/or NaN), log transform is meaningless
        sys.exit(
            f"Error: {trait_col} contains only zeros (and/or NaN). "
            f"Log transform is not meaningful."
        )

# Apply natural logarithm
# (use np.log10(y) for log base 10, or np.log1p(original_y) for log(x+1))
y_log = np.log(y)

# Replace the column with log-transformed values
df[trait_col] = y_log

# Save the transformed data to output file (tab-delimited with header)
df.to_csv(output_file, sep="\t", index=False, header=True)
