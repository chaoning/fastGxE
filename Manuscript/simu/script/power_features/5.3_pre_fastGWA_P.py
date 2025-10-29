import os

# Define the input and output folder paths
input_folder = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_baselineLD/"
output_folder = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_baselineLD_fastGWA/"

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop through all files in the input folder
for filename in os.listdir(input_folder):
    input_path = os.path.join(input_folder, filename)
    output_path = os.path.join(output_folder, filename)

    # Make sure it's a file (not a directory)
    if os.path.isfile(input_path):
        # Read all lines from the input file
        with open(input_path, "r") as infile:
            lines = infile.readlines()

        # Write to the output file, skipping the first line (header)
        with open(output_path, "w") as outfile:
            outfile.writelines(lines[1:])
