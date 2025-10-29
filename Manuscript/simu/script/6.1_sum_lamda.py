import pandas as pd

df_lst = []
h2_add = 30
h2_gxe = 5
for h2_nxe in [0, 5, 15, 25]:
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQ_plot/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_100000.csv"
    df = pd.read_csv(file, index_col=0)
    df = df.T
    df["h2_add"] = h2_add / 100
    df["h2_gxe"] = h2_gxe / 100
    df["h2_nxe"] = h2_nxe / 100
    df_lst.append(df)

h2_add = 30
h2_nxe = 15
for h2_add in [5, 30, 50]:
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQ_plot/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_100000.csv"
    df = pd.read_csv(file, index_col=0)
    df = df.T
    df["h2_add"] = h2_add / 100
    df["h2_gxe"] = h2_gxe / 100
    df["h2_nxe"] = h2_nxe / 100
    df_lst.append(df)


h2_add = 30
h2_nxe = 15
for h2_gxe in [1, 5, 10]:
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQ_plot/h2_add_{h2_add}.h2_gxe_{h2_gxe}.h2_nxe_{h2_nxe}.sample_100000.csv"
    df = pd.read_csv(file, index_col=0)
    df = df.T
    df["h2_add"] = h2_add / 100
    df["h2_gxe"] = h2_gxe / 100
    df["h2_nxe"] = h2_nxe / 100
    df_lst.append(df)

df = pd.concat(df_lst)

df.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.1_sum_lamda.csv")
