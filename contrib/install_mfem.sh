# bash file to install mfem
###
 # @Description: 
 # @version: 
 # @Author: luowen
 # @Date: 2023-10-19 18:03:19
 # @LastEditTime: 2023-10-19 22:00:49
 # @FilePath: /DC3DForward_CG/contrib/install_mfem.sh
### 

# set mpicxx
mpicxx=mpiicpc
#mpicxx=mpicxx # mpicxx

# # gslib
# cd mfem_package/gslib
# make clean
# make CC=mpicc MPI=1

# hypre-2.19.0
# cd ..
# rm -rf hypre
cd mfem_package
cd hypre-2.19.0/src
./configure --disable-fortran
make clean
make -j
cd ../..
ln -s hypre-2.19.0 hypre

# metis-5.1.0
cd metis-5.1.0
make config
make clean
make
sudo make install
# ln -s ../build/Linux-x86_64/libmetis/libmetis.a lib

# mfem-4.5
cd ../mfem-4.5
make clean
# make parallel -j 4  MPICXX=${mpicxx} MFEM_USE_MPI=YES MFEM_USE_METIS_5=YES METIS_DIR=@MFEM_DIR@/../metis-5.1.0 MFEM_USE_GSLIB=YES
make parallel -j 4  MPICXX=${mpicxx} MFEM_USE_MPI=YES MFEM_USE_METIS_5=YES METIS_DIR=@MFEM_DIR@/../metis-5.1.0 

