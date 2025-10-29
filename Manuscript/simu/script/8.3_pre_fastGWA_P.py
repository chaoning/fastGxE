import os
import glob

# Define the input and output folder paths
input_folder = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/"
output_folder = "/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_fastGWA/"

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Use glob to find files matching the pattern
# This matches files containing "rep_" followed by a number and ending with ".txt"
input_files = glob.glob(os.path.join(input_folder, "*rep_*[0-9].txt"))

# Debugging: Print the input files that were found
print("Matched files:", len(input_files), "files found.")

# Loop through all matched files
for input_path in input_files:
    filename = os.path.basename(input_path)  # Get the file name from the path
    output_path = os.path.join(output_folder, filename)

    # Make sure it's a file (not a directory)
    if os.path.isfile(input_path):
        # Read all lines from the input file
        with open(input_path, "r") as infile:
            lines = infile.readlines()

        # Write to the output file, skipping the first line (header)
        with open(output_path, "w") as outfile:
            outfile.writelines(lines[1:])
