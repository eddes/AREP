# AREP

In this folder you will find several python routines for comfort depicted below.

## VDI_PET_corrected.py
The corrected version of the PET thermal comfort index after the VDI norm (based on Djordje Spasic's python translation of the fortran code https://github.com/stgeorges)
 
## PET_fsolve_code.py
 a more compact and efficient version of the PET comfort index

  You are free to use the scripts above on the condition to cite the authors
     "E. Walther (AREP) and Q. Goestchel (ENS Paris-Saclay),  after D. Spasic's version".
     or alternately the article reference https://www.sciencedirect.com/science/article/pii/S0360132318301896

## steady_PET.py

A state-of-the art version of the steady PET index, as above. Major changes are 

[-] replacing the old school dichotomy method by scipy's fsolve (dividing by approx. 3 the function evaluations)
[-] writing the the docstring

## PET_transient.py

A work in progress version of the transient PET model.

## code_flux_integrale.py

A bit of code to compute the total flux on a standing invidiual, materialised as a cylinder, and including Perez diffuse/direct/circumsolar sky models.
