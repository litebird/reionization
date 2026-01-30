# Optical Depth, Reionization of the Universe and Neutrino Masses #

This is the repo for the "Optical Depth, Reionization of the Universe and Neutrino Masses" projet paper.

## Reionization history

We used two models as mock data (see `LB_fiducials_mocks`)

## Likelihood description

* ```Planck```:
  Planck data are spread into three different likelihoods following the original version of the Planck likelihoods.
  * ```Planck lowT```:
    likelihood based on the PR3 low-ℓ TT posterior distributions of measured $C_\ell$ fitted with a Wishart ($\alpha_\ell, \beta_\ell$). We then rescaled to the mock fiducial model with $\beta_\ell \longrightarrow \beta_\ell^{\rm mock} = (\alpha_\ell+1)C_{\ell,{\rm mock}}^{TT}$. Multipole range available $\ell \in [2,29]$.

  * ```Planck lowE```:
    Hamimeche&Lewis approximation using the Planck PR4 covariance matrix ([LoLLiPoP](https://github.com/planck-npipe/lollipop)) evaluated at the mock fiducial spectrum. Multipole range available $\ell \in [2,29]$.

  * ```Planck high-ℓ (TT,TE,EE)```:
    Gaussian likelihood for TT , TE and EE power spectra, with the full covariance matrix from the PR3 _Planck-lite_ likelihood (marginalised over foreground parameters). Multipole range available $\ell \in [30,2508]$.

* ```LiteBIRD```:
  Wishart likelihood for TT , TE, and EE spectra, using beam-deconvolved noise spectra including foreground residuals after ILC and an effective sky fraction fsky ∈ {0.7, 0.8} from E-mode group. Multipole range available (depending on the noise file).

* ```Simons-Observatory```:
  Wishart likelihood with noise curves including foreground residuals. Multipole range available (depending on the noise file).

* ```CMB-S4```:
  Wishart likelihood with noise curves including foreground residuals. Multipole range available (depending on the noise file).

## Datasets noise

Noises are described in `LB_fiducials_noises`

## Reionization models

We fit the data with 4 different models:
* `tanh`: tanh model with $z_{reio}$ ($\Delta z$ fixed)
* `tanhdz`:  tanh model with $z_{reio}$ and $\Delta z$
* `flex1`: FlexKnot model with 3 knots ( (z1=0,xe1=1.08), (z2,xe2=1.08), (z3,xe3=0) )
* `flex2`: FlexKnot model with 4 knots ( (z1=0,xe1=1.08), (z2,xe2=1.08), (z3,xe3), (z4,xe4=0) )

