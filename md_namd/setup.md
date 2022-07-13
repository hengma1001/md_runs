# NAMD3 simulations 

## Setup the container
```bash 
singularity pull namd3.sif docker://nvcr.io/hpc/namd:3.0-alpha11
```
## Parmed peeve
Parmed miswrites the number of REMARKs, default at 0, while there is 
no REMARK line written and it ends up skipping the NATOMS line, causing
namd error. 
