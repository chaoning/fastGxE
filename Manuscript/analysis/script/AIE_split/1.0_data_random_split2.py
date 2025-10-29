import random

def split_file_randomly_with_header(filename, output1, output2):
    with open(filename, 'r') as f:
        lines = f.readlines()

    header = lines[0]
    data_lines = lines[1:]

    total_lines = len(data_lines)
    indices = list(range(total_lines))
    random.shuffle(indices)

    half = total_lines // 2
    idx1 = set(indices[:half])
    idx2 = set(indices[half:])

    lines1 = [data_lines[i] for i in sorted(idx1)]
    lines2 = [data_lines[i] for i in sorted(idx2)]

    with open(output1, 'w') as f1:
        f1.write(header)
        f1.writelines(lines1)

    with open(output2, 'w') as f2:
        f2.write(header)
        f2.writelines(lines2)

import sys

split_file_randomly_with_header(sys.argv[1], sys.argv[2], sys.argv[3])
