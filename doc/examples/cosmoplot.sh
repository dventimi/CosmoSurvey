#!/bin/sh
# -z1  (0.01)       Redshift start
# -z2  (2.50)       Redshift end
# -nz  (10)         Number of redshift bins
# -l1  (44.0)       log10(L1) luminosity start, [L1] = ergs/s
# -l2  (46.0)       log10(L2) luminosity end, [L2] = ergs/s
# -nl  (10)         Number of luminosity bins
# -do  (12.0)       Survey solid angle dOmega, [dOmega] = steradians
# -fl  (1.25e-13)   Survey flux limit, [fl] = ergs/s/cm^2
# -w   (-1.0)       Dark Energy equation of state parameter w
# -om  (0.30)       Matter density parameter Omega_M
# -ol  (0.70)       Dark Energy density parameter Omega_Lambda
# -s8  (0.80)       Power spectrum normalization sigma_8
# -g   (0.55)       Linder growth index gamma
cosmoplot \
    -s8 0.9
