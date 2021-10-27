import numpy as np
import math


f=open("./filtered_data_backward.dat", "a") 
data = []
stat = []
uncorr_sys = []
corr_sys = []
bins = [0,1,2,3]
for bin in bins:
    data_, stat_, uncorr_sys_, corr_sys_ = np.genfromtxt("./data_backward.dat", usecols=(bin*4+0,bin*4+1,bin*4+2,bin*4+3),
        delimiter='±', unpack=True, missing_values='-', skip_header=1)
    data.append(data_*1e3)
    stat.append(stat_*1e3)
    uncorr_sys.append(uncorr_sys_*1e3)
    corr_sys.append(corr_sys_*1e3)

data = np.concatenate(data)
stat = np.concatenate(stat)
uncorr_sys = np.concatenate(uncorr_sys)
corr_sys = np.concatenate(corr_sys)

for i in range(0,36):
        if not math.isnan(data[i]):
            f.write(str(data[i]) + '\t' + str(stat[i]) + '\t' + str(uncorr_sys[i])+ '\t' + str(corr_sys[i]) + '\n') 
