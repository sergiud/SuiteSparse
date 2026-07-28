#!/bin/sh
#SBATCH --job-name=paru_and_umf
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --partition=bigmem
#SBATCH --mem=2900G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80

module load GCCcore/13.3.0
module load intel/2024a
module load imkl/2024.2.0
module load intel-compilers/2024.2.0
module load GMP/6.3.0
module load MPFR/4.2.1
module load CMake/3.29.3
module load Matlab/R2024a
module list

export OMP_MAX_ACTIVE_LEVELS=4
export OMP_NUM_THREADS=80,80
export OMP_DISPLAY_ENV=true
printenv | grep OMP

/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/rim.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b39_c30.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/twotone.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/std1_Jac2.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/mac_econ_fwd500.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/std1_Jac3.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b300_c1.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/lhr71.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/Zd_Jac2.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/av41092.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/Zd_Jac6.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/crashbasis.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/bbmat.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/Zd_Jac3.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/mc2depi.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b300_c2.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b678_c1.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b300_c3.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/Chebyshev4.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/Hamrle3.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/pre2.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b2052_c1.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/torso1.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b678_c2.mtx
/scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/TSOPF_RS_b2383.mtx

## too big:
## /scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/vas_stokes_4M.mtx
## /scratch/user/davis/SuiteSparse/ParU/build/paru_benchmark /scratch/user/davis/matrices/stokes.mtx

