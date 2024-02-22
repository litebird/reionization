import numpy as np
from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator, PchipInterpolator

# Constants
_c_ =  2.99792458e8
_Mpc_over_m_ = 3.085677581282e22
_G_ = 6.67428e-11
_m_H_ = 1.673575e-27
_sigma_ = 6.6524616e-29
_not4_ = 3.9715
_k_B_ = 1.3806504e-23
_h_P_ = 6.62606896e-34
sigma_B = 2.*(np.pi**5.)*(_k_B_**4.)/15./(_h_P_**3.)/(_c_**2)
T_cmb = 2.7255
def same(x):
    return x

# YHe interpolator, accurate at 1e-7 level for massless neutrinos,
# and 1e-4 for massive neutrinos (even for sum_mnu ~ 1 eV)
BBN_file = np.loadtxt("sBBN_2017.dat", skiprows=14)
omegab_BBN = BBN_file[:, 0][:48]
DeltaNeff_BBN = BBN_file[:, 1][::48]
YHe_BBN = BBN_file[:, 2].reshape(-1,48)
DeltaNeff_ref = 3.044 - 3.046
YHe_fct = RegularGridInterpolator((DeltaNeff_BBN, omegab_BBN),
                                  YHe_BBN, method="cubic")

# Create fHe function from YHe interpolator
def fHe_fct(omega_b):
    YHe = YHe_fct((DeltaNeff_ref, omega_b))
    return YHe/(_not4_ * (1. - YHe))

### Reionisation models
reio_he_z = 3.5
reio_he_dz = 0.5
reio_exp = 1.5

# Tanh + Delta z
def xe_tanhdz(z, omega_b, z_re, dz_re):
    fHe = fHe_fct(omega_b)
    arg = (((1.+z_re)**reio_exp - (1.+z[:, None])**reio_exp)
        / (reio_exp*((1.+z_re)**(reio_exp-1.))) / dz_re)
    xe = (1. + fHe) * (np.tanh(arg) + 1.) / 2.
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Tanh
def xe_tanh(z, omega_b, z_re):
    return xe_tanhdz(z, omega_b, z_re, 0.5*np.ones(z_re.shape))

# Flexknot with 1 knot, z_beg=30 and z_end=5
def xe_flex1(z, omega_b, z_f1, xe_f1, quiet=False):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    tqdm_or_not = same if quiet else tqdm
    for i in tqdm_or_not(range(len(omega_b))):
        x = [0., 4.9, 5., z_f1[i], 30., 100.]
        y = [1+fHe[i], 1+fHe[i], 1+fHe[i], xe_f1[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Flexknot with 2 knots, z_beg=30 and z_end=5
def xe_flex2(z, omega_b, z_f1, xe_f1, z_f2, xe_f2, quiet=False):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    tqdm_or_not = same if quiet else tqdm
    for i in tqdm_or_not(range(len(omega_b))):
        x = [0., 4.9, 5., z_f1[i], z_f2[i], 30., 100.]
        y = [1+fHe[i], 1+fHe[i], 1+fHe[i], xe_f1[i], xe_f2[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Flexknot with 1 knot, free z_beg and z_end
def xe_flexalt(z, omega_b, z_f1, z_f2, xe_f2, z_f3, quiet=False):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    tqdm_or_not = same if quiet else tqdm
    for i in tqdm_or_not(range(len(omega_b))):
        x = [0., 4.9, z_f1[i], z_f2[i], z_f3[i], 100.]
        y = [1+fHe[i], 1+fHe[i], 1+fHe[i], xe_f2[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Flexknot with 0 knot, free z_beg and z_end
def xe_flexaltalt(z, omega_b, z_f1, z_f2, quiet=False):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    tqdm_or_not = same if quiet else tqdm
    for i in tqdm_or_not(range(len(omega_b))):
        x = [0., z_f1[i], z_f2[i], 100.]
        y = [1+fHe[i], 1+fHe[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Compute tau(z) integrand
def tau_integ(z_arr, xe_arr, omega_b, H0, omega_cdm):
    h = H0/100
    YHe = YHe_fct((DeltaNeff_ref, omega_b))
    n_e = 3 * (H0 * 1e3 / _Mpc_over_m_)**2. * (omega_b / h**2.) / (8*np.pi*_G_*_m_H_) * (1-YHe)
    Omega_m = (omega_b + omega_cdm) / (H0/100)**2.
    Omega_g = (4.*sigma_B/_c_*(T_cmb**4.))/(3.*_c_*_c_*1e10*h*h/_Mpc_over_m_/_Mpc_over_m_/8./np.pi/_G_)
    Omega_ur = 3.044*7./8.*((4./11.)**(4./3.))*Omega_g
    Omega_r = Omega_g + Omega_ur
    Omega_l = 1. - Omega_m - Omega_g - Omega_r
    H_arr = H0 * 1e3 * np.sqrt(Omega_m * (1+z_arr)**3. + Omega_r*(1+z_arr)**4. + Omega_l) / _Mpc_over_m_
    comov = _c_ / H_arr
    return n_e * xe_arr * comov * (1+z_arr)**2. * _sigma_

# Compute tau(z) integrand
def tau_integ_bis(z_arr, xe_arr, omega_b, H0, omega_cdm):
    h = H0/100
    YHe = YHe_fct((DeltaNeff_ref, omega_b))
    n_e = 3 * (1e2 * 1e3)**2. / _Mpc_over_m_ * omega_b / (8*np.pi*_G_*_m_H_) * (1-YHe)
    Omega_m = (omega_b + omega_cdm) / h**2.
    Omega_g = (4.*sigma_B/_c_*(T_cmb**4.))/(3.*_c_*_c_*1e10*h*h/_Mpc_over_m_/_Mpc_over_m_/8./np.pi/_G_)
    Omega_ur = 3.044*7./8.*((4./11.)**(4./3.))*Omega_g
    Omega_r = Omega_g + Omega_ur
    Omega_l = 1. - Omega_m - Omega_g - Omega_r
    H_arr2 = (H0 * 1e3)**2. * (Omega_m * (1+z_arr)**3. + Omega_r*(1+z_arr)**4. + Omega_l)
    H_arr = np.sqrt(H_arr2)
    comov = _c_ / H_arr
    return n_e * xe_arr * (1+z_arr)**2. * _sigma_, comov
