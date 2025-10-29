import os

power_snp_index=-1
h2_add=30
h2_gxe=5
h2_nxe=15

os.chdir("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_take_home/")
k = 0
lst = []
for rep in range(1, 101):  # 1 to 100 inclusive
    for num_envi_power in [1, 2, 5, 10, 20, 30, 40]:
        for power_snp_h2_gxe in [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14]:
            for sample_size in [50000, 100000, 200000, 300000, 400000]:
                k += 1
                out_prefix = f"power_snp_h2_gxe_{power_snp_h2_gxe}.num_envi_{num_envi_power}.h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_{sample_size}.rep_{rep}"
                # check if output file already exists
                if not os.path.exists(f"{out_prefix}.txt"):
                    print(k, out_prefix)
                    lst.append(k)

print(",".join(map(str, lst)))
print(f"Total combinations checked: {len(lst)}")

with open("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/1.1_check.txt", "w") as f:
    f.write(",".join(map(str, lst)) + "\n")
