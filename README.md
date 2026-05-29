# polarimeter_recon_fit.cpp

A ROOT macro for proton polarimeter reconstruction and quality-control diagnostics for the GEp-V (GEP5) experiment at Jefferson Lab Hall A. It processes replayed data from the Super BigBite Spectrometer (SBS) GEM front tracker (FT) and focal-plane polarimeter (FPP) tracking detector systems, applies event selection cuts, and produces a suite of polarimeter reconstruction plots.

---

## Overview

The macro reads ROOT trees produced by the GEp-5 replay chain, selects clean elastic *ep* → *ep* events using a global cut string, and fills eight diagnostic histograms covering:

- FPP polar scattering angle (θ_FPP)
- Distance of closest approach (DOCA) between the FT and FPP tracks
- Reconstructed scattering vertex position along the beam (z_close)
- Angular differences between FT and FPP tracks (Δxp, Δyp)
- 2D correlation plots of Δxp vs. Δyp with and without the θ_FPP cut

Results are saved to `polarimeter_recon.pdf`.
