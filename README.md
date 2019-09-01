# Quantum-Mechanical-System
## Installation
```
ifort -O3 -o compute_vdw.bin -qopenmp compute_vdw.f90
ifort -O3 -o response_function.bin -qopenmp response_function.f90
ifort -O3 -o RPA.bin -qopenmp RPA.f90  -lblas -llapack
ifort -O3 -o integral.bin -qopenmp integral.f90  -lblas -llapack
```
## For developers

