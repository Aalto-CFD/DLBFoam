#!/bin/bash
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH -n 16

#set -euo pipefail

. /share/apps/spack/envs/fgci-centos7-generic/software/openfoam-org/7/beuz5lf/etc/bashrc
module load intel-mkl
module load hdf5/1.10.2



blockMesh
decomposePar -force

srun reactingFoam -parallel > log.reactingFoam




