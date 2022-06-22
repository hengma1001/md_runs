#!/usr/bin/env bash
# Usage: ./singularity.sh <gpu count>
set -e; set -o pipefail

GPU_COUNT=${1:-1}
IMG="namd3.sif"

INPUT="run.namd"

SINGULARITY="singularity exec --nv -B $(pwd):/host_pwd --pwd /host_pwd ${IMG}"
NAMD2="namd3 +p 20 +devices 0,1,2,3 ${INPUT}"

echo "Running APOA1 example in ${SIMG} on ${GPU_COUNT} GPUS..."
${SINGULARITY}  ${NAMD2}
