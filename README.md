# CoulexAngCorr
Imposes angular correlation of the projectile scattering angle and gamma-ray emission direction following Coulomb excitation

Maintainer: Aaron Chester

Durign Coulomb excitation, the projectile nucleus scattering angle and the emission direction of the subsequent gamma-ray are correlated and the gamma-ray emission direction is not isotropic. Further, this correlation is affected by the CM energy of the reaction, which varies as the projectile propagates through a thick reaction target. This code modifies Geant4-simulated spectra to include the effects of angular correlation using the Geant track weighting mechanism.

The program pre-computes a number of matrices which are used to calculate the weight imposed by the correlation function and then grabs the relavent information from the simulation output root files in order to determine the parameters of the correlation functoin on an event-by-event basis.

Useful References: K. Alder, A. Bohr, T. Huus, B. Mottelson, and A. Winther. Rev. Mod. Phys. 28 (1956) 432-542.; A. Chester. Ph. D. Thesis.
