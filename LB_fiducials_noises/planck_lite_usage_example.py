import numpy as np

#############################
### TO BE SET BY THE USER ###
#############################

# Read mock data vector from somewhere (model A or B in our case)
# Note: expected to be dimensioneless, no l(l+1) factor, ells=[0, 2508]
fid_cl_tt = ...
fid_cl_te = ...
fid_cl_ee = ...

# Choose the ell cut
ell_cut = ...

# Tcmb needed for conversion
T_cmb = ...


############################################################################
### CODE EXAMPLE: PART THAT NEED TO BE DONE ONLY ONCE, AT INITIALISATION ###
############################################################################

# Read Planck lite bin centers
av = np.loadtxt("noise_planck_bin_ref.dat")
ltt = av[:215]
lte = av[215:414]
lee = av[414:]

# Deduce Planck lite ell bins from centers, and keep only those above ell_cut
ixs_ells = []
# TT bins
ell_bins_tt = []
start_ell = 30
for i, mid_ell in enumerate(ltt):
    end_ell = 2 * int(mid_ell) - start_ell
    if start_ell >= ell_cut:
        ell_bins_tt.append(list(range(start_ell, end_ell + 1)))
        ixs_ells.append(i)
    start_ell = end_ell + 1
# TE bins
ell_bins_te = []
start_ell = 30
for i, mid_ell in enumerate(lte):
    end_ell = 2 * int(mid_ell) - start_ell
    if start_ell >= ell_cut:
        ell_bins_te.append(list(range(start_ell, end_ell + 1)))
        ixs_ells.append(len(ltt) + i)
    start_ell = end_ell + 1
# EE bins
ell_bins_ee = []
start_ell = 30
for i, mid_ell in enumerate(lee):
    end_ell = 2 * int(mid_ell) - start_ell
    if start_ell >= ell_cut:
        ell_bins_ee.append(list(range(start_ell, end_ell + 1)))
        ixs_ells.append(len(ltt) + len(lte) + i)
    start_ell = end_ell + 1

# Derive weight matrices for binning Cells
n_bins_tt, n_bins_te, n_bins_ee = len(ell_bins_tt), len(ell_bins_te), len(ell_bins_ee)
n_bins_tot = n_bins_tt + n_bins_te + n_bins_ee
weights_tt = np.zeros((2509, n_bins_tt))
weights_te = np.zeros((2509, n_bins_te))
weights_ee = np.zeros((2509, n_bins_ee))
for i, ell_bins in enumerate(ell_bins_tt):
    weights_tt[ell_bins, i] = 1. / len(ell_bins)
for i, ell_bins in enumerate(ell_bins_te):
    weights_te[ell_bins, i] = 1. / len(ell_bins)
for i, ell_bins in enumerate(ell_bins_ee):
    weights_ee[ell_bins, i] = 1. / len(ell_bins)

# Bin mock data vector
binned_fid_cl = np.hstack((
    fid_cl_tt @ weights_tt,
    fid_cl_te @ weights_te,
    fid_cl_ee @ weights_ee,
)) * 1e12 * T_cmb**2.

# Read, cut, and invert Planck lite covariance matrix
cov = np.loadtxt("noise_planck_lite.dat")
icov = np.linalg.inv(cov[ixs_ells, :][:, ixs_ells])


#############################################################
### CODE EXAMPLE: MCMC PART, TO BE DONE AT EACH MCMC STEP ###
#############################################################

# Theoretical Cells coming from CLASS or CAMB
# Note: expected to be dimensioneless, no l(l+1) factor, ells=[0, 2508]
cl_theo_tt = ...
cl_theo_te = ...
cl_theo_ee = ...

# Compute likelihood
binned_th_cl = np.hstack((
    cl_theo_tt @ weights_tt,
    cl_theo_te @ weights_te,
    cl_theo_ee @ weights_ee,
)) * 1e12 * T_cmb**2.
diff = binned_th_cl - binned_fid_cl
lnlike = -0.5 * diff @ icov @ diff
