import numpy as np
import matplotlib.pyplot as plt

f = np.loadtxt("../LB_forecast_Aachen/data/noise_litebird_planck_v2.dat")
E_noise_in_muKarcmin = 6.56
arcmin_to_radian = np.pi / 180. / 60
ells = np.arange(2, 1001)

# 30 arcmin.muK noises
theta_in_arcmin = 30.
N_TT_30 = (arcmin_to_radian * E_noise_in_muKarcmin / np.sqrt(2.))**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_30 = (arcmin_to_radian * E_noise_in_muKarcmin)**2. * (
    np.exp(ells * (ells+1) * (arcmin_to_radian * theta_in_arcmin)**2. / 8. / np.log(2.)))
N_EE_30[:27] = f[:27, 2]



