import numpy as np
import math

other_bins = [0,1,3,4]
f=open("13TeV_ratios.dat", "a")

#D0 meson
pt_ref, ptmin, ptmax, data_ref, statp_ref, statm_ref, sysp_ref, sysm_ref = np.genfromtxt(
    "./D0_bin2.dat", usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')

for dd in other_bins:
    pt, ptmin, ptmax, data, statp, statm, sysp, sysm = np.genfromtxt("./D0_bin" + str(dd) + ".dat", 
        usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')
    ratio = data/data_ref
    statp = ratio*np.sqrt( statp**2/data**2 + statp_ref**2/data_ref**2 )
    sysp = ratio*np.sqrt( sysp**2/data**2 + sysp_ref**2/data_ref**2 )
    statm = -ratio*np.sqrt( statm**2/data**2 + statm_ref**2/data_ref**2 )
    sysm = -ratio*np.sqrt( sysm**2/data**2 + sysm_ref**2/data_ref**2 )

    for i in range(0,11):
        if not math.isnan(ratio[i]):
            f.write(str(pt[i]) + '\t' + str(ratio[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i])
                + '\t' + str(sysp[i]) + '\t' +  str(sysm[i]) + '\n')

#Dp meson
pt_ref, ptmin, ptmax, data_ref, statp_ref, statm_ref, sysp_ref, sysm_ref = np.genfromtxt(
    "./Dp_bin2.dat", usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')

for dd in other_bins:
    pt, ptmin, ptmax, data, statp, statm, sysp, sysm = np.genfromtxt("./Dp_bin" + str(dd) + ".dat", 
        usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')
    ratio = data/data_ref
    statp = ratio*np.sqrt( statp**2/data**2 + statp_ref**2/data_ref**2 )
    sysp = ratio*np.sqrt( sysp**2/data**2 + sysp_ref**2/data_ref**2 )
    statm = -ratio*np.sqrt( statm**2/data**2 + statm_ref**2/data_ref**2 )
    sysm = -ratio*np.sqrt( sysm**2/data**2 + sysm_ref**2/data_ref**2 )

    for i in range(0,11):
        if not math.isnan(ratio[i]):
            f.write(str(pt[i]) + '\t' + str(ratio[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i])
                + '\t' + str(sysp[i]) + '\t' +  str(sysm[i]) + '\n')


#Dps meson
pt_ref, ptmin, ptmax, data_ref, statp_ref, statm_ref, sysp_ref, sysm_ref = np.genfromtxt(
    "./Dps_bin2.dat", usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')

for dd in other_bins:
    pt, ptmin, ptmax, data, statp, statm, sysp, sysm = np.genfromtxt("./Dps_bin" + str(dd) + ".dat", 
        usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')
    ratio = data/data_ref
    statp = ratio*np.sqrt( statp**2/data**2 + statp_ref**2/data_ref**2 )
    sysp = ratio*np.sqrt( sysp**2/data**2 + sysp_ref**2/data_ref**2 )
    statm = -ratio*np.sqrt( statm**2/data**2 + statm_ref**2/data_ref**2 )
    sysm = -ratio*np.sqrt( sysm**2/data**2 + sysm_ref**2/data_ref**2 )

    for i in range(0,10):
        if not math.isnan(ratio[i]):
            f.write(str(pt[i]) + '\t' + str(ratio[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i])
                + '\t' + str(sysp[i]) + '\t' +  str(sysm[i]) + '\n')            
