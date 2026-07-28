#!/bin/sh
#SBATCH --job-name=mumps
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

echo "Matrix: rim.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/rim.mtx

echo "Matrix: TSOPF_RS_b39_c30.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b39_c30.mtx

echo "Matrix: twotone.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/twotone.mtx

echo "Matrix: std1_Jac2.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/std1_Jac2.mtx

echo "Matrix: mac_econ_fwd500.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/mac_econ_fwd500.mtx

echo "Matrix: std1_Jac3.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/std1_Jac3.mtx

echo "Matrix: TSOPF_RS_b300_c1.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b300_c1.mtx

echo "Matrix: lhr71.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/lhr71.mtx

echo "Matrix: Zd_Jac2.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/Zd_Jac2.mtx

echo "Matrix: av41092.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/av41092.mtx

echo "Matrix: Zd_Jac6.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/Zd_Jac6.mtx

echo "Matrix: crashbasis.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/crashbasis.mtx

echo "Matrix: bbmat.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/bbmat.mtx

echo "Matrix: Zd_Jac3.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/Zd_Jac3.mtx

echo "Matrix: mc2depi.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/mc2depi.mtx

echo "Matrix: TSOPF_RS_b300_c2.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b300_c2.mtx

echo "Matrix: TSOPF_RS_b678_c1.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b678_c1.mtx

echo "Matrix: TSOPF_RS_b300_c3.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b300_c3.mtx

echo "Matrix: Chebyshev4.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/Chebyshev4.mtx

echo "Matrix: Hamrle3.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/Hamrle3.mtx

echo "Matrix: pre2.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/pre2.mtx

echo "Matrix: TSOPF_RS_b2052_c1.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b2052_c1.mtx

echo "Matrix: torso1.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/torso1.mtx

echo "Matrix: TSOPF_RS_b678_c2.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b678_c2.mtx

echo "Matrix: TSOPF_RS_b2383.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/TSOPF_RS_b2383.mtx

echo "Matrix: vas_stokes_4M.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/vas_stokes_4M.mtx

echo "Matrix: stokes.mtx"
/scratch/user/davis/MUMPS_573_our/examples/c_mumps_paru /scratch/user/davis/matrices/stokes.mtx

