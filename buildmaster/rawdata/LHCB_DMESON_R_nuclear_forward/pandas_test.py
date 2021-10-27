import pandas as pd

df = pd.read_csv("./FwdRpPb2D.txt", sep=" ", header=None, dtype=float, names= ["pt","y","ratio","sys_corr","sys_uncorr"] )
df = df[(df[['ratio']] != 0).all(axis=1)] # drop values with ratio=0
sort_y_pt= df.sort_values(by=["y","pt"])
sort_y_pt.to_csv("./data_Fwd.csv", sep=" ", header=False, index=False, float_format='%.3f')
