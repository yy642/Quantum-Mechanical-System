# Quantum-Mechanical-System
## Installation
required library:LAPACK/BLAS library
```
ifort -O3 -o compute_vdw.bin -qopenmp compute_vdw.f90
ifort -O3 -o response_function.bin -qopenmp response_function.f90
ifort -O3 -o RPA.bin -qopenmp RPA.f90  -lblas -llapack
ifort -O3 -o integral.bin -qopenmp integral.f90  -lblas -llapack
```
## For developers
**compute_vdw.f90** : compute dispersion energy given xyz file using openmp.
see "Tkatchenko, Alexandre, and Matthias Scheffler.Physical review letters 102, no. 7 (2009): 073005."

**response_function.f90** : compute the response function chi for harmonic oscillator 

**RPA.f90** : under random phase approximation, using the pre-computed response function, compute the eigenvalue for chi.T, where chi is the precomputed response function and T is the interaction tensor. 
Detials see supporting info of "Yang, Yan, Ka Un Lao, and Robert A. DiStasio Jr. Physical review letters 122, no. 2 (2019): 026001."

**integral.f90** : using casimir polder integral to compute the dispersion energy.
Detials see supporting info of "Yang, Yan, Ka Un Lao, and Robert A. DiStasio Jr. Physical review letters 122, no. 2 (2019): 026001."
