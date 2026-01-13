# Fiducial models for reionization studies

This folder contains the lensed and unlensed CMB spectra produced by CLASS (modified for reionization [CLASS_REIO](https://github.com/s-ilic/class_reio/tree/main)).

• `model_A`: tanh model

$x_e(y) = \frac{f}{2} \left[ 1 + \tanh \left( \frac{y(z_\text{reio}) - y}{\Delta_y} \right) \right]$

where $z_\text{reio}$ is the reionization mid-point, $y(z_\text{reio}) = (1 + z_\text{reio})^{3/2}$, $x_e = f/2$, and $\Delta y = 1.5 \sqrt{1 + z_{\text{reio}}} \Delta z$.

Default: $\Delta z = 0.5$, $\tau =0.06$

• `model_B`: exp asymmetric model

$x_e(z) = \frac{f}{2} \left[ (z_{beg} - z) / (z_{beg} - z_{end}) \right]^\alpha_{re}$

Default: $\alpha_{re} = 6$, $z_{end} = 5.$, $z_{beg} = 26.456663499594907$ so that $\tau = 0.06$.

