#!/bin/bash
set -e

# Script arguments
GPU_COUNT=1 # ${1:-1}
SIMG=/lambda_stor/homes/heng.ma/Pkgs/gmx/gromacs-2020_2.sif

# mdp file for gmx
MDP_PATH=/homes/heng.ma/Research/force_field/gmx_mdps
EM_MDP=${MDP_PATH}/minim.mdp

# inputs 
input_path=input
pdb_file=${input_path}/ncd_dimers_0_2_10nm.pdb
top_file=${input_path}/ncd_dimers_0_2_10nm.top

mkdir -p em
cd em
cp -r ../input ./

# Set number of OpenMP threads
# export OMP_NUM_THREADS=6 # ${OMP_NUM_THREADS:-1}

# Singularity will mount the host PWD to /host_pwd in the container
SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${SIMG}"

# Prepare benchmark data
${SINGULARITY} gmx grompp \
                -f ${EM_MDP} \
                -c ${pdb_file} \
                -p ${top_file}\
                -o em.tpr

# export GMX_GPU_DD_COMMS=true
# export GMX_GPU_PME_PP_COMMS=true
# export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export CUDA_VISIBLE_DEVICES=3

${SINGULARITY} gmx mdrun -v -deffnm em 
