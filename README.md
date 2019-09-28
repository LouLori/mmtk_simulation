# mmtk_simulation
example running script for MMTK simulations (PIMC/PIMD)

## nmv_prop
This is an asymmetric top density matrix calculator which is based on Fortran code.

### compilation
```
make clean
make
```
### calculate and generate density matrix
To optimize the calculation of density matrix, we run *submit_the.sh* script on server with *qsub* for a fixed theta value 181 times parallezation calculation as the range of theta angle is (0,180]:
```
for c in {0..180}; do qsub -v "ithe=$c" submit_the.sh; done
```

then, 181 files will be generated, we can use the *compile.x* script to combine all the files to get the final density matrix file:
```
./compile.x
```

final files are named as **rho.den\_rho, rho.den\_eng, rho.den\_esq**

## H2O Simulations
### Sample Forcefields


## CO2 Simulations
