REIO forecast
=============

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


``REIO forecast`` is a project for forecasting constraints on the reionization history for LiteBIRD

The likelihoods are centered on a fiducial model and use different approximations according to the dataset.

Likelihoods available are ``PLK``, ``PLK_lowlT``, ``PLK_lowlE``, ``LB``, ``SO``, ``S4`` for CMB, ``SOlens`` and ``S4lens`` for CMB lensing and ``bao`` for BAO.

It is interfaced with the ``cobaya`` MCMC sampler.

Install
-------

Then you can install the `REIO_forecast` likelihoods *via*

```shell
$ pip install -e .
```

Requirements
------------
* Python >= 3.5
* `numpy`
* `astropy`
* `cobaya`

