#!/bin/sh
# -w   (-1.0)       Dark Energy equation of state parameter w
# -om  (0.30)       Matter density parameter Omega_M
# -ol  (0.70)       Dark Energy density parameter Omega_Lambda
# -s8  (0.80)       Power spectrum normalization sigma_8
# -g   (0.55)       Linder growth index gamma
cosmotest \
    -w -0.5 \
    -om 0.70 \
    -ol 0.30 \
    -s8 0.90 \
    -g 0.55