# This script combines theory spectra Cl and noise spectra Nl into mock observed spectra Cl+Nl,
# including conversion to format of MontePython fiducial files
# Sabarish Venkartaramani, Julien lesgourgues, 10.22

import numpy as np
from math import pi

##############################################################################
#
# input files:
#fiducial_wo_noise = 'modelA/model_A_cl_lensed.dat'
fiducial_wo_noise = 'modelB/model_B_cl_lensed.dat'
noise = '../LB_Nldd_from_FuturCMB/noise-for-mcmc/noise_litebird_only_b50.dat'
#
# output files:
#fiducial_w_noise = 'modelA/cl_modelA_w_noise_b50.dat'
fiducial_w_noise = 'modelB/cl_modelB_w_noise_b50.dat'
fiducial_header = "# model B with unkown parameters (blind test) plus noise, in format of MP fiducial files"
#
##############################################################################

cl_wo_noise = np.loadtxt(fiducial_wo_noise)

# T_cmb in muK
T_cmb = 2.7255e6

twopi=2.*pi

clll = cl_wo_noise[:,0]
clTT = cl_wo_noise[:,1]
clEE = cl_wo_noise[:,2]
clTE = cl_wo_noise[:,3]
clBB = cl_wo_noise[:,4]
clpp = cl_wo_noise[:,5]
clTp = cl_wo_noise[:,6]
clEp = cl_wo_noise[:,7]

print ('-> fiducial input:',fiducial_wo_noise)
print ('   with: l_min=',clll[0],' l_max=',clll[-1],' number of ls=',len(clll))

l_size_fid = len(clll)

# convert Cl's of temp./pol. from dimensionless to [muK]**2
# convert to Cldd
clTT *= T_cmb**2 * twopi / (clll * (clll+1))
clEE *= T_cmb**2 * twopi / (clll * (clll+1))
clTE *= T_cmb**2 * twopi / (clll * (clll+1))
clBB *= T_cmb**2 * twopi / (clll * (clll+1))
cldd = clll * (clll+1) * clpp #TBC!!
clTd = T_cmb * np.sqrt(clll * (clll+1)) * clTp #TBC!

nl = np.loadtxt(noise)

nlll = nl[:,0]
nlTT = nl[:,1]
nlPP = nl[:,2]
nldd = nl[:,3]

# convert to Nldd
nldd *= twopi / nlll / (nlll+1.) #TBC!

print ('-> noise input:',noise)
print ('   with: l_min=',nlll[0],' l_max=',nlll[-1],' number of ls=',len(nlll))

l_size = len(nlll)

if (l_size_fid < l_size):
    print ("ERROR: the script assumes that l_max of Cl >= l_max of Nl")
else:
    print('-> the output Cl+Nl fiducial file will be written in:')
    print('  ',fiducial_w_noise)
    print('   up to l_max=',clll[l_size-1])
    print('   with header:')
    print('  ',fiducial_header)

cl_w_noise = np.zeros((l_size,6),'float64')

cl_w_noise[:,0] = clll[:l_size]
cl_w_noise[:,1] = clTT[:l_size] + nlTT[:]
cl_w_noise[:,2] = clEE[:l_size] + nlPP[:]
cl_w_noise[:,3] = clTE[:l_size]
cl_w_noise[:,4] = cldd[:l_size] + nldd[:]
cl_w_noise[:,5] = clTd[:l_size]

np.savetxt(fiducial_w_noise, cl_w_noise, header=fiducial_header, fmt='%d %.8e %.8e %.8e %.8e %.8e')
