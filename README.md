# Optical Depth, Reionization of the Universe and Neutrino Masses #

This is the repo for the "Optical Depth, Reionization of the Universe and Neutrino Masses" projet paper.


## Likelihood description

* ```Planck```:
  Planck data are spread into three different likelihoods following the original version of the Planck likelihoods.
  * ```Planck lowT```:
    likelihood based on the PR3 low-ℓ TT posterior distributions of measured $C_\ell$ fitted with a Wishart ($\alpha_\ell, \beta_\ell$). We then rescaled to the mock fiducial model with $\beta_\ell \longrightarrow \beta_\ell^{\rm mock} = C_{\ell,{\rm mock}}^{TT}\,(\alpha_\ell+1)$. Multipole range available $\ell \in [2,29]$.

  * ```Planck lowE```:
    Hamimeche&Lewis approximation using the Planck PR4 covariance matrix (LoLLiPoP) evaluated at the mock fiducial spectrum. Multipole range available $\ell \in [2,29]$.

  * ```Planck high-ℓ (TT,TE,EE)```:
    Gaussian likelihood for TT , TE and EE power spectra, with the full covariance matrix from the PR3 lite likelihood. Multipole range available $\ell \in [30,2508]$.

* ```LiteBIRD```:
  Wishart likelihood for TT , TE, and EE spectra, using beam-deconvolved noise spectra including foreground residuals after ILC and an effective sky fraction fsky ∈ {0.7, 0.8} from E-mode group. Multipole range available $\ell \in [2,2500]$.

* ```Simons-Observatory```:
  Wishart likelihood with noise curves including foreground residuals. Multipole range available $\ell \in [127,7979]$.

* ```CMB-S4```:
  Wishart likelihood with noise curves including foreground residuals. Multipole range available $\ell \in [127,10000]$.

  
