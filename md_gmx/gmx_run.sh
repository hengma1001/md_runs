#!/bin/bash
set -e
# Usage: ./water.sh {GPU_COUNT} {IMG_NAME}

# Script arguments
GPU_COUNT=1 # ${1:-1}
SIMG=/lambda_stor/homes/heng.ma/Pkgs/gmx/gromacs-2020_2.sif

# mdp file for gmx
MDP_PATH=/homes/heng.ma/Research/force_field/gmx_mdps
mdp_file=${MDP_PATH}/nvt.mdp

# inputs 
input_path=input
pdb_file=${input_path}/ncd_dimers_0_2_10nm.pdb
top_file=${input_path}/ncd_dimers_0_2_10nm.top

mkdir -p md 
cd md
cp -r ../input ./
cp ../em/em.gro ./

# Set number of OpenMP threads
export OMP_NUM_THREADS=6 # ${OMP_NUM_THREADS:-1}


# Singularity will mount the host PWD to /host_pwd in the container
SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${SIMG}"
 
# Prepare benchmark data
${SINGULARITY} gmx grompp \
                -f ${mdp_file} \
                -c em.gro \
                -p ${top_file}\
                -o nvt.tpr

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5

# Run benchmark
${SINGULARITY} gmx mdrun -v -deffnm md \
                 -ntmpi ${GPU_COUNT} \
                 -nb gpu -pme gpu \
                 -bonded gpu \
                 -npme 1 \
                 -ntomp ${OMP_NUM_THREADS} \
                 -nstlist 400 \
                 -s nvt.tpr


