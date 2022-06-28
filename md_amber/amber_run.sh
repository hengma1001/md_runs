#!/bin/bash 

GPU_COUNT=4
# source /homes/ramanathana/amber20/amber.sh
source /homes/heng.ma/Pkgs/amber20/amber.sh

# inputs 
input_path=input
crd_file=${input_path}/ncd_dimers_0_2_10nm_em.inpcrd
top_file=${input_path}/ncd_dimers_0_2_10nm.prmtop

rm -r md
mkdir -p md 
cd md
cp ../md_npt.in ./
cp -r ../input ./

amber_in=md_npt.in
amber_out=md_npt.out

export CUDA_VISIBLE_DEVICES=0,1,2,3
mpirun -np ${GPU_COUNT} pmemd.cuda_SPFP.MPI -i ${amber_in} -o ${amber_out} \
    -p ${top_file} -c ${crd_file} -r md.nc 
    # -ref ${crd_file} \
