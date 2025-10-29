import pandas as pd
import numpy as np
import os

k = 0
lst = []
for rep in range(100):
    for h2_gxe in [1, 5, 10]:
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_30.h2_gxe_{h2_gxe}.h2_nxe_15.sample_100000.rep_{rep+1}.var"
        k += 1
        if not os.path.exists(file):
            lst.append(str(k))
        
        file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_mom/h2_add_30.h2_gxe_{h2_gxe}.h2_nxe_15.sample_100000.rep_{rep+1}.noNxE.var"
        k += 1
        if not os.path.exists(file):
            lst.append(str(k))

print((len(lst)))
print(",".join(lst))
