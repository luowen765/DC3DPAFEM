
# DC3DPAFEM user guide

DC3DPAFEM is a parallel adaptive finite-element algorithm for 3D direct current resistivity anisotropic forward. Details about the theory will be found in the article: _A scalable parallel finite-element algorithm using the algebraic multigrid solver for 3-D direct current resistivity modeling in anisotropic media_ currently in review.

## Prerequisites

The DC3DPAFEM was developed in C++ using MPI on the Linux operating system. Because we employ some compilers and libraries from the intel compiler and make, you should first install the intel compiler and make it. 

The Intel's **oneAPI toolkit** need to be installed including Base Toolkit and HPC Toolkit. They can be downloaded referring to https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html and https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html.

## Program description

This DC3DPAFEM is open source and freely available. It contains 4 folders:

* **contrib**: this folder shall contain needed third-party libraries. The open source finite element library MFEM, partitioning finite element meshes program METIS, and high performance preconditioners and solvers library HYPRE were used in this code. The MFEM can be available on the web: https://github.com/mfem/mfem, with the METIS available in https://github.com/KarypisLab/METIS , ans with the HYPRE available in https://github.com/hypre-space/hypre. Tested libraries and corresponding versions are __hypre-2.19.0__, __metis-5.1.0__, and __mfem-4.5__.You need to download these libraries and install them in sequence by *./install_mfem.sh*, in which:
  * *Installation for HYPRE*
    ```
    tar -zxvf hypre-2.19.0.tar.gz
    cd hypre-2.19.0/src/
    ./configure --disable-fortran
    make -j
    cd ../..
    ln -s hypre-2.19.0 hypre
    ```
  * *Installation for METIS*
    ``` 
    tar -zvxf metis-5.1.0.tar.gz
    cd metis-5.1.0
    make BUILDDIR=lib config
    make BUILDDIR=lib
    cp lib/libmetis/libmetis.a lib
    cd ..
    ```
  * *Installation for MFEM*
    ```
    tar zvxf mfem-4.5.tar.gz
    cd mfem-4.5/
    make parallel -j 4  MPICXX=${mpicxx} MFEM_USE_MPI=YES MFEM_USE_METIS_5=YES METIS_DIR=@MFEM_DIR@/../metis-5.1.0 
    ```

- **examples**: this folder contains two examples. The first example contains the Accuracy verification, Scalability test, and Robunstness test. The second example contains the Consistence test and Effect of topography and anisotropic test. Note that in the Consistency test of the second model, the model file mountainvalley.tar.xz need to be decompressed and restored as a .msh file. We show the configuration files in the Accuracy verification test, and the other tests are similar. It contains the following files:

  * *model.config*
    This file contains all the configuration options required by the DC3DPAFEM program. Commonly used ones are, for example, electrode configuration file name, resistivity parameter file name, number of adaptive refining, indication factor, grid file name, and so on.
  * *sigma.file*
    This file describes how many regions exist and the six separate parameters corresponding to the resistivity tensor for each region
  * *survey.input*
    This file describes the arrangement of the electrodes. The two numbers in the first line indicate the number of measurements and the pole-pole configuration, respectively. Each subsequent line denotes the coordinates of the source and measurement electrodes.
  * *twolay.msh*
    This file describes the computed 3D model.

  * *makefile*, *run.sh* 
    The former is used to compile the program to produce the executable file DC3D. The latter is used to run the program.

- **src**: this folder contains all sources file (including .cpp and .h).
Inside this file all source files have introductory comments giving the authors name, date of completion and IO etc.
- _DC3D.cpp_
  * This file is the main program of the package DC3DPAFEM, serving
 as the overall implementation framework for the algorithm. It includes reading
 input files, domain decomposition, solving finite element linear equation
 systems, and the mesh adaptation process. For more details, please refer to the
 article: _A scalable parallel finite-element algorithm using the algebraic
 multigrid solver for 3-D direct current resistivity modeling in anisotropic
 media_.
- _em.h em.cpp_
    * The namespace of the package DC3DPAFEM includes a variety of public functions and variables.
- _error\_estimators.h error\_estimators.cpp_
  * This class implements a posteriori error estimation based on continuity of current density (Ren et al., 2018b). _Ren, Z., QIU, L., Tang, J., ZHOU, F., CHEN, C., CHEN, H., HU, S., 2018b. 3D modeling of direct-current anisotropic resistivity using the adaptive finite-element method based on continuity of current density. Chinese Journal of Geophysics 61, 331–343_.
- _goafem.h goafem.cpp_
  * This class implements a framework for the goal-oriented finite element method, including finite element discretization, assembly and solving of linear equation systems, a posteriori error estimation, and mesh refinement. 
- _meshplus.h meshplus.cpp_
    * This class inherits MFEM's class in mesh/mesh.hpp. I do this for extending the print function. This file is part of the MFEM library. 
- _mfem\_coefficient.h mfem\_coefficient.cpp_
    * The following classes are the direct copies and modifies of MFEM's classes in fem/coefficient.hpp. I just did a minor modifications to compute different coeffients for assembling equations. These classes are called by class GOAFEM in goafem.h.
- _parahandler.h parahandler.cpp_
    * This class is used for post-processing the solution of FEM, converting the potential into apparent resistivity.
- _solversplus.h solversplus.cpp_
    * This class inherits MFEM's class in linalg/solvers.h. I do this for modify the print part of CGSolver::Mult function. 

- **figures**: this folder contains the figures presented in the article and corresponds to the data computed in the provided examples and tests.


## Running the tests

We show the Accuracy verification test in the two-layer model here.

```
cd DC3DPAFEM/examples/two-layer
```
Copy all files from folder Accuracy_verification to current directory
```
cp Accuracy_verification/* .
sh run.sh
```
The resulting files are saved in the **solutions** folder, where solution.x represents the potential and apparent resistivity results.

## License

This project is licensed under the terms of the BSD-3 license. See file LICENSE for details.



## Authors

* **Lewen Qiu**
School of Geosciences and Info-Physics, Central South University, Changsha, 410083, Hunan, China.
E-mail: luo_wen@csu.edu.cn

* **Zhengguang Liu**
School of Mathematics and Statistics, Central South University, Changsha, 410083, Hunan, China.
E-mail: zhengguang-liu@outlook.com

* **Hongbo Yao**
Macau lnstitute of Space Technology and Application, Taipa, Macao, 999078, China.
E-mail: hongbo.yao@outlook.com

* **Jingtian Tang**
School of Geosciences and Info-Physics, Central South University, Changsha, 410083, Hunan, China.
E-mail: jttang@csu.edu.cn




