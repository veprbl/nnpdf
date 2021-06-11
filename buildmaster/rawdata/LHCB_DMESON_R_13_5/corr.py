import numpy as np

def read_triangular_matrix(file, n):
    mat = np.empty((n,n),dtype=float)
    vec = np.genfromtxt(file, usecols=(4), unpack=True, delimiter=';', missing_values='-')
    vec = vec[~np.isnan(vec)]
    N = 0
    for i in range(0,n):
        for k in range(i,n):
            mat[i,k] = vec[N]
            mat[k,i] = vec[N]
            N+=1
    return(mat)

def read_matrix(file, n, m):
    mat = np.empty((n,m),dtype=float)
    vec = np.genfromtxt(file, usecols=(4), unpack=True, delimiter=';', missing_values='-')
    vec = vec[~np.isnan(vec)]
    N = 0
    for i in range(0,n):
        for k in range(0,m):
            mat[i,k] = vec[N]
            N+=1
    return(mat)

f = open("corr_tot.dat", "a")

npoints_D0 = 38
npoints_Dp = 37
npoints_Dps = 32

corr_D0 = read_triangular_matrix("./D0.dat", npoints_D0)
corr_Dp = read_triangular_matrix("./Dp.dat", npoints_Dp)
corr_Dps = read_triangular_matrix("./Dps.dat", npoints_Dps)

#eigenvalues, eigenvect = np.linalg.eig(corr_Dp)
#print(eigenvalues)

np.savetxt("corr_Dp.dat", corr_Dp)



corr_D0_Dp = read_matrix("./D0_Dp.dat", npoints_D0, npoints_Dp)
corr_D0_Dps = read_matrix("./D0_Dps.dat", npoints_D0, npoints_Dps)
corr_Dp_Dps = read_matrix("./Dp_Dps.dat", npoints_Dp, npoints_Dps)


l1 = np.concatenate((corr_D0,corr_D0_Dp, corr_D0_Dps), axis=1)
l2 = np.concatenate((np.transpose(corr_D0_Dp),corr_Dp,corr_Dp_Dps),axis=1)
l3 = np.concatenate((np.transpose(corr_D0_Dps),np.transpose(corr_Dp_Dps),corr_Dps),axis=1)
corr_tot = np.concatenate((l1,l2,l3),axis=0)

ll1 = np.concatenate((corr_D0,corr_D0_Dps), axis=1)
ll2 = np.concatenate((np.transpose(corr_D0_Dps),corr_Dps),axis=1)
corr_tot_test = np.concatenate((ll1,ll2),axis=0)

eigenvalues, eigenvect = np.linalg.eig(corr_tot_test)
print(eigenvalues)

for i in range(0,107):
    for j in range(0,107):
        f.write(str(corr_tot[i,j]) + '\t')
    f.write('\n')

f.close()

#np.savetxt("corr_tot.dat", corr_tot)





