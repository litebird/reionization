We post here the files thyat are necessary to reproduce (or compare with) our forecast runs.

* litebird_planck_v1 likelihood and litebird_planck_v1_lcdm chains (23.11.2021):

    - likelihood: LiteBird + Planck where:
        - the LiteBird noise NlEE comes from gluing 3 pieces: from simulations at low l, assuming a sensitivity of 6.56 muK*arcmin at intermediate l, assuming a 30 arcmin beam at large l, and fsky = 0.57
        - the Planck noise NlTT and NlEE comes from the fake_planck_realistic likelihood with fsky = 0.57
        - the LiteBird noise NlEE replaces the Planck one as long as it is smaller, and lmax=3000
        - we kept the option of CMB lensing extraction, LensingExtraction = True, using the noise Nldd from fake_planck_realistic
        - further details can be checked in montepyhton/likelihoods/litebird_planck_v1

    -model: minimal LambdaCDM with the parameters discussed in previous telecons and Mnu fixed to 3*0.02eV
        - we get sigma(tau_reio) = 0.0022
        - further detials can be seen in input/litebird_planck_v1_lcdm.ini (input) and in chains/litebird_planck_v1_lcdm/... (output)

* litebird_planck_v1 likelihood and litebird_planck_v1_mnu chains (24.11.2021):

    - likelihood: same likelihood and fiducial model as above

    - model: LambdaCDM with a varying m_ncdm, and thus a total mass Mnu = 3 * m_ncdm
        - we get the same sigma(tau_reio) = 0.0022 as for fixed Mnu
        - we get no lower bound on m_ncdm, just an upper bound m_ncdm < 0.11 eV (95%CL)
        - this translates into Mnu < 0.33 eV (95%CL)
