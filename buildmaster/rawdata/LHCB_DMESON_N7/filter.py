import numpy as np

pt, ptmin, ptmax, data, statp, statm, sysp, sysm = np.genfromtxt(
    "/home/tommy/physics/nnpdfgit/nnpdf/buildmaster/rawdata/LHCB_DMESON_N7/data.csv",
     usecols=(0,1,2,3,4,5,6,7), unpack=True, delimiter=',', missing_values='-')

D0_bin2 = data[16:24]
D0_statp_bin2 = statp[16:24]
D0_statm_bin2 = statp[16:24]
D0_sysp_bin2 = sysp[16:24]
D0_sysm_bin2 = sysm[16:24]


Dp_bin2 = data[56:64]
Dps_bin2 = data[96:104]

D0 = np.concatenate((data[:8]/D0_bin2, data[8:16]/D0_bin2, data[24:32]/D0_bin2, data[32:40]/D0_bin2))
D0_statp = np.concatenate((statp[:8], statp[8:16], statp[24:32], statp[32:40]))
D0_statm = np.concatenate((statm[:8], statm[8:16], statm[24:32], statm[32:40]))
D0_sysp = np.concatenate((sysp[:8], sysp[8:16], sysp[24:32], sysp[32:40]))
D0_sysm = np.concatenate((sysm[:8], sysm[8:16], sysm[24:32], sysm[32:40]))


Dp = np.concatenate((data[40:48]/Dp_bin2, data[48:56]/Dp_bin2, data[64:72]/Dp_bin2, data[72:80]/Dp_bin2))
Dp_statp = np.concatenate((statp[40:48], statp[48:56], statp[64:72], statp[72:80]))
Dp_statm = np.concatenate((statm[40:48], statm[48:56], statm[64:72], statm[72:80]))
Dp_sysp = np.concatenate((sysp[40:48], sysp[48:56], sysp[64:72], sysp[72:80]))
Dp_sysm = np.concatenate((sysm[40:48], sysm[48:56], sysm[64:72], sysm[72:80]))


Dps = np.concatenate((data[80:88]/Dps_bin2, data[88:96]/Dps_bin2, data[104:112]/Dps_bin2, data[112:120]/Dps_bin2))
Dps_statp = np.concatenate((statp[80:88], statp[88:96], statp[104:112], statp[112:120]))
Dps_statm = np.concatenate((statm[80:88], statm[88:96], statm[104:112], statm[112:120]))
Dps_sysp = np.concatenate((sysp[80:88], sysp[88:96], sysp[104:112], sysp[112:120]))
Dps_sysm = np.concatenate((sysm[80:88], sysm[88:96], sysm[104:112], sysm[112:120]))


delete_bins = [30,31,32,62,63,64,72,79,80,86,87,88,89,92,93,94,95]

data = np.concatenate((D0,Dp,Dps))
data = np.delete(data,delete_bins)
print(data)

statp = np.concatenate((D0_statp, Dp_statp, Dps_statp))
statp = np.delete(statp,delete_bins)
statm = np.concatenate((D0_statm, Dp_statm, Dps_statm))
statm = np.delete(statm,delete_bins)
sysp = np.concatenate((D0_sysp, Dp_sysp, Dps_sysp))
sysp = np.delete(sysp,delete_bins)
sysm = np.concatenate((D0_sysm, Dp_sysm, Dps_sysm))
sysm = np.delete(sysm,delete_bins)


f=open("7TeV_ratios.dat", "a")
for i in range(0,data.size):
    f.write( str(data[i]) + '\t' + str(statp[i]) + '\t' + str(statm[i]) + '\t' + str(sysp[i]) + '\t' + str(sysm[i]) +  '\n')

f.close()
