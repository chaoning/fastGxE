import numpy as np
import os

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_mmsusie/")

lst = []
for i in range(1, 11):
    for index in range(1, 139):
        try:
            file_path = f"{i}.{index}.mmsusie_out.pip.txt"
            arr = np.loadtxt(file_path, dtype=float)
            lst.append(np.sum(arr > 0.95))
        except Exception as e:
            pass


print(len(lst))
print(np.mean(lst))

lst = []
for i in range(1, 11):
    for index in range(1, 139):
        try:
            file_path = f"{i}.{index}.susie.pip"
            arr = np.loadtxt(file_path, dtype=float)
            lst.append(np.sum(arr > 0.95))
        except Exception as e:
            pass


print(len(lst))
print(np.mean(lst))
