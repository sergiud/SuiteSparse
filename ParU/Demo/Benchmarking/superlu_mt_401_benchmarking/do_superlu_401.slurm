#!/bin/sh
#SBATCH --job-name=slu_401_v2
#SBATCH --time=08:00:00
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

echo "Matrix: rim,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/rim.rb

echo "Matrix: TSOPF_RS_b39_c30,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b39_c30.rb

echo "Matrix: twotone,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/twotone.rb

echo "Matrix: std1_Jac2,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/std1_Jac2.rb

echo "Matrix: mac_econ_fwd500,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/mac_econ_fwd500.rb

echo "Matrix: std1_Jac3,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/std1_Jac3.rb

echo "Matrix: TSOPF_RS_b300_c1,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b300_c1.rb

echo "Matrix: lhr71,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/lhr71.rb

echo "Matrix: Zd_Jac2,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/Zd_Jac2.rb

echo "Matrix: av41092,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/av41092.rb

echo "Matrix: Zd_Jac6,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/Zd_Jac6.rb

echo "Matrix: crashbasis,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/crashbasis.rb

echo "Matrix: bbmat,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/bbmat.rb

echo "Matrix: Zd_Jac3,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/Zd_Jac3.rb

echo "Matrix: mc2depi,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/mc2depi.rb

echo "Matrix: TSOPF_RS_b300_c2,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b300_c2.rb

echo "Matrix: TSOPF_RS_b678_c1,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b678_c1.rb

echo "Matrix: TSOPF_RS_b300_c3,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b300_c3.rb

echo "Matrix: Chebyshev4,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/Chebyshev4.rb

echo "Matrix: Hamrle3,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/Hamrle3.rb

echo "Matrix: pre2,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/pre2.rb

echo "Matrix: TSOPF_RS_b2052_c1,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b2052_c1.rb

echo "Matrix: torso1,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/torso1.rb

echo "Matrix: TSOPF_RS_b678_c2,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b678_c2.rb

echo "Matrix: TSOPF_RS_b2383,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/TSOPF_RS_b2383.rb

echo "Matrix: vas_stokes_4M,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/vas_stokes_4M.rb

echo "Matrix: stokes,"
/scratch/user/davis/superlu_mt_401_our/build/EXAMPLE/pdlinsolx_paru < /scratch/user/davis/RB/stokes.rb

