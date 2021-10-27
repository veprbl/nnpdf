import numpy as np
import math

f=open("./data_pp.dat", "a") 
data  = []
statp = []
statm = []
sysp = []
sysm = []

bins = [1,2,3,4]
for bin in bins:
    pt, ptmin, ptmax, data_, statp_, statm_, sysp_, sysm_ = np.genfromtxt("./data_pp/bin" + str(bin) + ".dat", usecols=(0,1,2,3,4,5,6,7),
        delimiter=',', unpack=True, missing_values='-')
    data.append(data_)
    statp.append(statp_)
    statm.append(statm_)
    sysp.append(sysp_)
    sysm.append(sysm_)

data = np.concatenate(data)
statp = np.concatenate(statp)
statm = np.concatenate(statm)
sysp = np.concatenate(sysp)
sysm = np.concatenate(sysm)


for i in range(0,36):
        if not math.isnan(data[i]):
            f.write(str(data[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i]) + '\t' + str(sysp[i]) + '\t' + str(sysm[i]) + '\n') 

