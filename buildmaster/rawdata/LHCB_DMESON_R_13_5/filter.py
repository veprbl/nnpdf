#filter hepdata tables and print data in a c++ friendly format
import numpy as np

f=open("5TeV_ratios.dat", "a")

pt, ptmin, ptmax, data, statp, statm, sysp, sysm = np.loadtxt(
    "/home/tommy/physics/Dmeson/HEPdata/ratio_5TeV/data.csv", delimiter=",", 
    skiprows=2, usecols=(0,1,2,3,4,5,6,7), unpack=True)

for i in range(0,pt.size):
    f.write(str(pt[i]) + '\t' + str(data[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i]) + '\t' + str(sysp[i]) + '\t' + str(sysm[i]) +  '\n')

f.close()

