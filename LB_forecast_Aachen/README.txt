We post here the files thyat are necessary to reproduce (or compare with) our forecast runs.

* litebird_planck_v1 likelihood and litebird_planck_v1_lcdm chains (23.11.2021):

    - likelihood: LiteBird + Planck where:
        - the LiteBird noise NlEE comes from gluing 3 pieces: from simulations at low l, assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, and fsky = 0.57
        - the Planck noise NlTT and NlEE comes from the fake_planck_realistic likelihood with fsky = 0.57
        - the LiteBird noise NlEE replaces the Planck one as long as it is smaller, and lmax=3000
        - we kept the option of CMB lensing extraction, LensingExtraction = True, using the noise Nldd from fake_planck_realistic
        - further details can be checked in montepyhton/likelihoods/litebird_planck_v1

    - model: minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get sigma(tau_reio) = 0.0022
        - further details can be seen in input/litebird_planck_v1_lcdm.ini (input) and in chains/litebird_planck_v1_lcdm/... (output)

* litebird_planck_v1 likelihood and litebird_planck_v1_mnu chains (24.11.2021):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - we get the same sigma(tau_reio) = 0.0022 as for fixed Mnu
        - we get no lower bound on m_ncdm, just an upper bound m_ncdm < 0.11 eV (95%CL)
        - this translates into Mnu < 0.33 eV (95%CL)

* litebird_planck_v2 likelihood and litebird_planck_v2_lcdm chains (3.12.2021):

    - likelihood: LiteBird + Planck, differing from litebird_planck_v1 only through the fact that fsky is taken to be 0.70 in the LiteBird range 2 <= l <= 600 and 0.57 in the Planck range l > 600
                  (the l-dependent f_sky is written in the last column of the noise file and is read by the likelihood_mock_cmb code inside montepython/likelihood_class.py,
                  which required a small modification of this file)

    - model: minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get slightly smaller errors, but the impact is minor, and rounding at two digits we still get sigma(tau_reio) = 0.0022
        - further details can be seen in input/litebird_planck_v2_lcdm.ini (input) and in chains/litebird_planck_v2_lcdm/... (output)

* litebird_planck_v2 likelihood and litebird_planck_v2_mnu chains (3.12.2021):

    - likelihood: same likelihood and fiducial model as above with varying f_sky = 0.70 -> 0.57

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get the same sigma(tau_reio) = 0.0022 as for fixed Mnu
        - the neutrino mass bound is almost the same as in the v1 case, Mnu < 0.32 eV (95%CL)
        - further details can be seen in input/litebird_planck_v2_mnu.ini (input) and in chains/litebird_planck_v2_mnu/... (output)

* litebird_w_lens likelihood and litebird_w_lens_lcdm chains (updated 23.12.2021):

    - likelihood: we sketeched a litebird-alone likelihood in the following way:
        - the LiteBird noise NlEE comes from the simulations at low l, assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, l_max=1350, and fsky = 0.7
          (so this is like in the litebird_planck_v2 case but without the Planck noise at large l)
        - the LiteBird noise NLTT and Nldd comes from the paper 1808.05955, that is, from the file "noise_litebird.dat" of the public MontePython.
        - note that NlTT does not really matter since we are comsic variance dominated at low l. Nldd does matter and we will do another test where we remove lensing extraction, for the sake of comparison.

    - model:  minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get once more sigma(tau_reio) = 0.0022, which is identical to the results 1808.05955, despite of the fact that we have slightly degraded NlEE at very low l.
        - further details can be seen in input/litebird_w_lens_lcdm.ini (input) and in chains/litebird_w_lens_lcdm/... (output)

* litebird_w_lens likelihood and litebird_w_lens_mnu chains (23.12.2021):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get the same sigma(tau_reio) = 0.0022 as for fixed Mnu
        - the neutrino mass bound degrades a bit w.r.t Litebird+Planck, Mnu < 0.36 eV (95%CL)
        - further details can be seen in input/litebird_w_lens_mnu.ini (input) and in chains/litebird_planck_v2_mnu/... (output)

* litebird_wo_lens likelihood and litebird_wo_lens_lcdm chains (updated 23.12.2021):

    - likelihood: the only difference w.r.t. litebird_w_lens is the fact that we turn off lensing extraction (LensingExtraction = False)

    - model:  minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get once more sigma(tau_reio) = 0.0022, which suggests that this result comes entirely from the sensitivity to polarisation (+temperature) and is unaffected by lensing
        - further details can be seen in input/litebird_wo_lens_lcdm.ini (input) and in chains/litebird_wo_lens_lcdm/... (output)

* litebird_wo_lens likelihood and litebird_wo_lens_mnu chains (12.12.2021):

    - likelihood: same as above (litebird only without lensing extraction), wrong fiducial file, needs to be updated

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - after rounding, we get the sigma(tau_reio) = 0.0023, so lensing extraction has a tiny impact when Mnu fluctuates (with lensing extraction we had 0.022)
        - the neutrino mass bound, Mnu < 0.81 eV (95%CL), degrades considerably in absence of lensing extraction (it was 0.32 eV with lensing extraction)
        - further details can be seen in input/litebird_wo_lens_mnu.ini (input) and in chains/litebird_wo_lens_mnu/... (output)
