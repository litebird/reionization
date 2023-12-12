import numpy as np
from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator, PchipInterpolator


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
_not4_ = 3.9715
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

# Flexknot with 1 knot
def xe_flex1(z, omega_b, z_f1, xe_f1):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    for i in tqdm(range(len(omega_b))):
        x = [0., 4.9, 5., z_f1[i], 30., 31.]
        y = [1+fHe[i], 1+fHe[i], 1+fHe[i], xe_f1[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe

# Flexknot with 1 knot
def xe_flex2(z, omega_b, z_f1, xe_f1, z_f2, xe_f2):
    fHe = fHe_fct(omega_b)
    # Deal with knots
    xe = np.zeros((len(z), len(omega_b)))
    for i in tqdm(range(len(omega_b))):
        x = [0., 4.9, 5., z_f1[i], z_f2[i], 30., 31.]
        y = [1+fHe[i], 1+fHe[i], 1+fHe[i], xe_f1[i], xe_f2[i], 0., 0.]
        interp = PchipInterpolator(x, y)
        xe[:, i] = interp(z)
    # Deal with HeII
    arg = (reio_he_z - z[:, None]) / reio_he_dz
    xe += fHe * (np.tanh(arg)+1.) / 2.
    return xe
