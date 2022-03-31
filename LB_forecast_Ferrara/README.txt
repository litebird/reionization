Here is a list of the completed runs:

* litebird_planck_v1 likelihood and litebird_planck_v1_lcdm chains:

	- likelihood: LiteBird + Planck where:
        - the LiteBird noise NlEE comes from gluing 3 pieces: foregrounds residual from component separation at low l (deconvolved by a 80 arcmin beam), assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, and fsky = 0.57
        - the Planck noise NlTT and NlEE comes from the data/noise_fake_planck_realistic_two.dat file (reproducing fairly well the accuracy with which Planck constrains the cosmo parameters), which is used in the fake_planck_realistic likelihood with fsky = 0.57
        - the LiteBird noise NlEE replaces the Planck one as long as it is smaller, and lmax=3000
        - we kept the option of CMB lensing extraction, LensingExtraction = True, using the noise Nldd from fake_planck_realistic
        - further details can be checked in montepyhton/likelihoods/litebird_planck_v1

	- model: minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get sigma(tau_reio) = 0.0024
        - further details can be seen in input/litebird_planck_v1_lcdm.param (input) and in chains/litebird_planck_v1_lcdm/... (output)

* litebird_planck_v1 likelihood and litebird_planck_v1_mnu chains:

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - we get sigma(tau_reio) = 0.0024
        - we get no lower bound on m_ncdm, just an upper bound m_ncdm < 0.11 eV (95%CL)
        - this translates into Mnu < 0.33 eV (95%CL)

* litebird_w_lens_lcdm_fsky60_beam30 likelihood and litebird_w_lens_lcdm_fsky60_beam30 chains:

    - likelihood: we sketeched a litebird-alone likelihood in the following way:
        - the LiteBird noise NlEE comes from foregrounds residual from component separation at low l (deconvolved by a 80 arcmin beam), assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, l_max=1350, and fsky = 0.6
        - the LiteBird noise NLTT and Nldd comes from the paper 1808.05955, that is, from the file "noise_litebird.dat" of the public MontePython.
        - note that NlTT does not really matter since we are cosmic variance dominated at low l. Nldd does matter and we will do another test where we remove lensing extraction, for the sake of comparison.

    - model:  minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get sigma(tau_reio) = 0.0025 (from getdist,  we actually obtain sigma(tau_reio) = 0.0023)
        - further details can be seen in input/litebird_w_lcdm_fsky60_beam30.param (input) and in chains/litebird_w_lens_lcdm_fsky60_beam30/... (output)

* litebird_w_lens_mnu_fsky60_beam30 likelihood and litebird_w_lens_mnu_fsky60_beam30 chains:

    - likelihood: same likelihood as above, fiducial: data/litebird_w_lens_mnu_fsky60_beam30.param

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get sigma(tau_reio) = 0.0023 (like for fixed Mnu, analyzing with getdist)
        - the neutrino mass bound degrades a bit w.r.t Litebird+Planck, mnu < 0.129 --> Mnu < 0.39 eV (95%CL)
        - further details can be seen in input/litebird_w_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_w_lens_mnu_fsky60_beam30/... (output)


* litebird_w_lens_mnu_fsky60_beam30 likelihood and litebird_w_lens_mnu_fsky70_beam30_import_sampl chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get sigma(tau_reio) = 0.0022, (from getdist, sigma(tau_reio) = 0.0021, ~7% smaller than the fsky = 0.6 case)
        - the neutrino mass bound degrades a bit w.r.t Litebird+Planck, mnu < 0.12 --> Mnu < 0.36 eV (95%CL)
        - further details can be seen in input/litebird_w_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_w_lens_mnu_fsky70_beam30_import_sampl/... (output)
	The results on sigma(tau_reio) and Mnu are consistent with the run litebird_w_lens_mnu chains of the Aachen group (Litebird only, fsky = 0.7, beam = 30')

* litebird_w_lens_mnu_fsky60_beam30 likelihood and litebird_w_lens_mnu_fsky80_beam30_import_sampl chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get sigma(tau_reio) = 0.0020, ~13% smaller than the fsky = 0.6 case
        - the neutrino mass bound is mnu < 0.112 --> Mnu < 0.33 eV like in the case with Litebird+Planck (95%CL)
        - further details can be seen in input/litebird_w_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_w_lens_mnu_fsky80_beam30_import_sampl/... (output)

* litebird_wo_lens_mnu_fsky60_beam30 likelihood and litebird_wo_lens_mnu_fsky60_beam30 chains:
     - likelihood: the only difference w.r.t. litebird_w_lens_mnu_* is the fact that we turn off lensing extraction (LensingExtraction = False), fiducial: data/litebird_wo_lens_mnu_fsky60_beam30.param
     - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
         - we get sigma(tau_reio) = 0.0024 
         - the neutrino mass bound is mnu < 0.295 (95% CL) --> Mnu < 0.88
         - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_wo_lens_mnu_fsky60_beam30_import_sampl/... (output)

* litebird_wo_lens_mnu_fsky60_beam30 likelihood and litebird_wo_lens_mnu_fsky70_beam30_import_sampl chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above
    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 0.0022
        - mnu < 0.277 (95% CL) ---> Mnu < 0.83
        - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_wo_lens_mnu_fsky70_beam30_import_sampl/... (output)

     The results on sigma(tau_reio) is consistent with the run litebird_wo_lens_mnu chains of the Aachen group (Litebird only, fsky = 0.7, beam = 30'), while Mnu constraint is slightly larger (after rounding, I get mnu < 0.277 with montepython and mnu < 0.273 with getdist. By rounding as mnu < 0.277 I get Mnu < 0.83 instead of 0.81 like the Aachen run)

* litebird_wo_lens_mnu_fsky60_beam30 likelihood and litebird_wo_lens_mnu_fsky80_beam30_import_sampl chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above
    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 0.0021
        - mnu < 0.26 (95% CL) ---> Mnu < 0.78
        - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_wo_lens_mnu_fsky80_beam30_import_sampl/... (output)

sigma(tau_reio) constraints are slightly worse in the cases w/o lensing extraction with respect to the corresponding cases w/ lens., while the Mnu constraints degrade much more.


* litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30 and tau_prior likelihoods and litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30 chains:
     
    - likelihoods: same fiducial model as before, setting litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30.OnlyTT = True to select only temperature; also litebird_TT...beam30.LensingExtraction = True. The tau_prior likelihood sets a gaussian prior on tau with mean = fiducial value used for tau (0.0543) and as sigma the one from the run litebird_w_lens_mnu_fsky60_beam30 (0.0023287) 

!!! Notice that to work with both OnlyTT and LensingExtraction set to True, there is a small bug in montepython/likelihood_class.py to be corrected: in the block `elif self.LensingExtraction:` at line 1575 one needs to add a `if self.OnlyTT:` block with the correct computation of the covariance matrix and then `else:` the `Cov_obs` and `Cov_the` expressions as computed now
  
    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 0.00245 (with getdist, 68% limit is 0.0023)
        - mnu < 0.278 (95% CL) ---> Mnu < 0.83
        - further details can be seen in input/litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30.param (input) and in chains/litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30/... (output)
 

* litebird_TT_wo_lens_tau_prior_mnu_fsky60_beam30 and tau_prior likelihoods and litebird_TT_wo_lens_tau_prior_mnu_fsky60_beam30 chains:

    - likelihoods: same fiducial model as before, setting litebird_TT_wo_lens_tau_prior_mnu_fsky60_beam30.OnlyTT = True to select only temperature; also litebird_TT...beam30.LensingExtraction = False. The tau_prior likelihood sets a gaussian prior on tau with mean = fiducial value used for tau (0.0543) and as sigma the one from the run litebird_wo_lens_mnu_fsky60_beam30 (0.0024052)

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 0.0025 (with getdist, 68% limit is 0.0024)
        - mnu < 0.458 (95% CL) ---> Mnu < 1.37
        - further details can be seen in input/litebird_TT_wo_lens_tau_prior_mnu_fsky60_beam30.param (input) and in chains/litebird_TT_wo_lens_tau_prior_mnu_fsky60_beam30/... (output)

* litebird_TT_w_lens_tau_prior_lcdm_fsky60_beam30 and tau_prior likelihoods and litebird_TT_w_lens_tau_prior_lcdm_fsky60_beam30 chains:

    - likelihoods: same fiducial model as before, setting litebird_TT_w_lens_tau_prior_lcdm_fsky60_beam30.OnlyTT = True to select only temperature; also litebird_TT...beam30.LensingExtraction = True. The tau_prior likelihood sets a gaussian prior on tau with mean = fiducial value used for tau (0.0543) and as sigma the one from the run litebird_w_lens_lcdm_fsky60_beam30 (0.0024636)

    - model: minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV 
        - sigma(tau_reio) =  0.00249 (with getdist, 68% limit is 0.0024)
        - further details can be seen in input/litebird_TT_w_lens_tau_prior_lcdm_fsky60_beam30.param (input) and in chains/litebird_TT_w_lens_tau_prior_lcdm_fsky60_beam30/... (output)

* litebird_w_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_w_lens_mnu_fsky60_beam30_new_noise chains:

    - likelihood: we sketeched a litebird-alone likelihood in the following way:
        - the LiteBird noise NlEE comes from foregrounds residual from component separation at low l (deconvolved by a 80 arcmin beam), assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, l_max=1350, and fsky = 0.6
        - the LiteBird noise NLTT and Nldd comes from the new noise computed by Sabarish using FuturCMB, now in the file "noise_litebird_only_b30.dat"

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm        
        - we get sigma(tau_reio) = 2.3584e-03 ~ 0.0024
        - mnu < 1.3125e-01 ~ 0.13 (95% CL) --> Mnu < 0.39 eV (95%CL), worse than in the previous case with a wrong noise curve
        - further details can be seen in input/litebird_w_lcdm_fsky60_beam30_new_noise.param (input) and in chains/litebird_w_lens_lcdm_fsky60_beam30_new_noise/... (output)

* litebird_w_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_w_lens_mnu_fsky70_beam30_import_sampl_new_noise chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get sigma(tau_reio) = 0.0022
        - mnu < 0.12 --> Mnu < 0.36 eV (95%CL)
        - further details can be seen in input/litebird_w_lens_mnu_fsky60_beam30_new_noise.param (input) and in chains/litebird_w_lens_mnu_fsky70_beam30_import_sampl_new_noise/... (output)

* litebird_w_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_w_lens_mnu_fsky80_beam30_import_sampl_new_noise chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get sigma(tau_reio) = 2.0758e-03 ~ 0.0021
        - the neutrino mass bound is mnu < 0.113 --> Mnu < 0.34 eV 
        - further details can be seen in input/litebird_w_lens_mnu_fsky60_beam30_new_noise.param (input) and in chains/litebird_w_lens_mnu_fsky80_beam30_import_sampl_new_noise/... (output)

* litebird_wo_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_wo_lens_mnu_fsky60_beam30_new_noise chains:
     - likelihood: the only difference w.r.t. litebird_w_lens_mnu_* is the fact that we turn off lensing extraction (LensingExtraction = False), fiducial: data/litebird_wo_lens_mnu_fsky60_beam30_new_noise.param
     - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
         - we get sigma(tau_reio) = 2.4038e-03 ~ 0.0024
         - the neutrino mass bound is mnu < 0.287 (95% CL) --> Mnu < 0.86
         - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30_new_noise.param (input) and in chains/litebird_wo_lens_mnu_fsky60_beam30_import_sampl_new_noise/... (output)

* litebird_wo_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_wo_lens_mnu_fsky70_beam30_import_sampl_new_noise chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above
    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 2.2446e-03 ~ 0.0023
        - mnu < 0.27 (95% CL) ---> Mnu < 0.81
        - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30_new_noise.param (input) and in chains/litebird_wo_lens_mnu_fsky70_beam30_import_sampl_new_noise/... (output)

* litebird_wo_lens_mnu_fsky60_beam30_new_noise likelihood and litebird_wo_lens_mnu_fsky80_beam30_import_sampl_new_noise chains (obtained from the fsky = 0.6 case with importance sampling):

    - likelihood: same likelihood and fiducial model as above
    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 2.1187e-03 ~ 0.0021
        - mnu < 0.2566 ~ 0.26 (95% CL) ---> Mnu < 0.77
        - further details can be seen in input/litebird_wo_lens_mnu_fsky60_beam30.param (input) and in chains/litebird_wo_lens_mnu_fsky80_beam30_import_sampl/... (output)

* litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30_new_noise and tau_prior likelihoods and litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30_new_noise chains:

    - likelihoods: same fiducial model as before, setting litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30_new_noise.OnlyTT = True to select only temperature; also litebird_TT...beam30_new_noise.LensingExtraction = False. The tau_prior likelihood sets a gaussian prior on tau with mean = fiducial value used for tau (0.0543) and as sigma the one from the run litebird_w_lens_mnu_fsky60_beam30_new_noise (0.0023584)

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - sigma(tau_reio) = 2.4309e-03 ~ 0.0024
        - mnu < 2.9511e-01 ~ 0.30 (95% CL) ---> Mnu < 0.90
        - further details can be seen in input/litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30_new_noise.param (input) and in chains/litebird_TT_w_lens_tau_prior_mnu_fsky60_beam30_new_noise/... (output)




(I am also copypasting this comment of Julien on the neutrino parameters used in those runs):
the meaning of all CLASS input parameters is explained at length in the file explanatory.ini which comes with any distribution of CLASS. In this case: T_ncdm is the neutrino-to-photon temperature ratio. T_ncdm = 0.71611 is the value that is needed in order to get [M_nu / omega_nu] = 93.14 eV, which is the prediction of  the most precise studies of neutrino decoupling. N_ur means “number of ultra-relativistic species”, and you may see it just as the number of massless neutrino species. If the neutrino phase-space distributions predicted by neutrino decoupling were exactly thermal (Fermi-Dirac) distributions, we would just need to take deg_ncdm=3 species of neutrinos  with a temperature set by T_ncdm=0.71611 and with the same mass, then, N_ur=0, and we would be done. However, precise studies of neutrinos decoupling predict deviations from a pure thermal distribution. In order to mimick these distortions, it has been shown that we can introduce a small fake density of relativistic species on top  of the massive ones. The advantage of setting T_ncdm = 0.71611 for the 3 massive neutrinos, and on top of this, N_ur=0.00441, is that you will then get simultaneously [M_nu / omega_nu] = 93.14 eV AND N_eff=3.044, wich are precisely the predictions of the neutrino decoupling model.
