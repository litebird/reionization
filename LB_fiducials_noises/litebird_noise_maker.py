import numpy as np
import matplotlib.pyplot as plt

f = np.loadtxt("../LB_forecast_Aachen/data/noise_litebird_planck_v2.dat") # for the low-ell EE noise
E_noise_in_muKarcmin = 6.56
arcmin_to_radian = np.pi / 180. / 60
ells = np.arange(2, 2501)

# deconvolve noises with 20 arcmin
theta_in_arcmin = 20.
N_TT_20 = (arcmin_to_radian * E_noise_in_muKarcmin / np.sqrt(2.))**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_20 = (arcmin_to_radian * E_noise_in_muKarcmin)**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_20[:27] = f[:27, 2]
np.savetxt(
    "noise_litebird_b20.dat",
    np.column_stack((ells, N_TT_20, N_EE_20)),
    fmt=["%i", "%.10e", "%.10e"],
    header="All N_ell in µK^2\nell   |   N_ell_TT   |   N_ell_EE",
)

# deconvolve noises with 30 arcmin
theta_in_arcmin = 30.
N_TT_30 = (arcmin_to_radian * E_noise_in_muKarcmin / np.sqrt(2.))**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_30 = (arcmin_to_radian * E_noise_in_muKarcmin)**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_30[:27] = f[:27, 2]
np.savetxt(
    "noise_litebird_b30.dat",
    np.column_stack((ells, N_TT_30, N_EE_30)),
    fmt=["%i", "%.10e", "%.10e"],
    header="All N_ell in µK^2\nell   |   N_ell_TT   |   N_ell_EE",
)



# Add residuals
noise_30 = np.loadtxt("noise_litebird_b30.dat")
res_file = "out_cls_nilc_baseline_fwhm30acm_fsky80_maskP.dat"
noise_pP_30_f80 = np.loadtxt(f"litebird_residuals/{res_file}", skiprows=2)
noise_30[:1023, 2] = noise_pP_30_f80[:, 3]
np.savetxt(
    f"noise_litebird_nilc_b30_f80.dat",
    noise_30,
    header=f"ell   |   N_ell_TT   |   N_ell_EE\n{res_file}")


noise_30 = np.loadtxt("noise_litebird_b30.dat")
res_file = "out_cls_nilc_baseline_fwhm30acm_fsky70_maskP.dat"
noise_pP_30_f80 = np.loadtxt(f"litebird_residuals/{res_file}", skiprows=2)
noise_30[:1023, 2] = noise_pP_30_f80[:, 3]
np.savetxt(
    f"noise_litebird_nilc_b30_f70.dat",
    noise_30,
    header=f"ell   |   N_ell_TT   |   N_ell_EE\n{res_file}")
