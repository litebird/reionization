This folder contains input and data files (fiducial spectra and noise) to run montepython forecast using TTTEEE assuming high noise in temperature to mimic forecast with EE only.

exact likelihood
TTTEEE
2<=ell<=30
fsky=1,0.7,0.5
fiducial spectra generated with Planck2018 TTTEEE+lowE+lensing bestfit (actually, mean of the posterior of) cosmological parameters -> same as EE only pipeline
noise and residual foregrounds from component separation for EE, N_ell=10^5 muK^2 for TT and dd
sampling over tau, with all params fixed to fiducial. Note we fix clamp=A_s exp(-2tau) instead of logA, so that A_s is recomputed at each mcmc step.

Note that the noise file has two columns repeated just to keep the expected format (TT,EE,dd).
Noise curves are downloaded from NERSC and correspond to the residual noise and foreground from the component separation analysis of the fg JSG (NERSC location: /project/projectdirs/litebird/results/compsep_products_from_PTEP_simulations/fgbuster_results/).

Note that noise curves are well below the cosmological signal at the scales we are looking at (2<=ell<=30) and the mcmc results are basically those of a CVL experiment.

Results:
fsky   sigma(tau)
0.5      0.0030
0.7      0.0025
1.0      0.0021
