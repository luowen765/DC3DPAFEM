
# DC3DPAFEM user guide

DC3DPAFEM is a parallel adaptive finite-element algorithm for 3D direct current resistivity anisotropic forward. Details about the theory will be found in the article: _A scalable parallel finite-element algorithm using the algebraic multigrid solver for 3-D direct current resistivity modeling in anisotropic media_ currently in review.

## Prerequisites

The DC3DPAFEM was developed using MPI on the Linux operating system. Because we employ some compliers and libraries from the intel complier and mkl. So you should firstly install the intel complier and mkl. 

The Intel's **oneAPI toolkit** need to be installed including Base Toolkit and HPC Toolkit. They can be downloaded referring to https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html and https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html.

## Program description

This DC3DPAFEM is open source and freely available. It contains 3 folders:

* _contrib_: this folder shall contain needed third-party libraries. The open source finite element library MFEM, partitioning finite element meshes program METIS, and high performance preconditioners and solvers library HYPRE were used in this code. The MFEM can be available on the web: https://github.com/mfem/mfem, with the METIS available in https://github.com/KarypisLab/METIS , ans with the HYPRE available in https://github.com/hypre-space/hypre. Tested libraries and corresponding versions are __hypre-2.19.0__, __metis-5.1.0__, and __mfem-4.5__.
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
- _examples_: this folder contains two examples. The first example contains the Accuracy verification, Scalability test, and Robunstness test. The second example contains the Consistence test and Effect of topography and anisotropic test. Note that in the Consistency test of the second model, the model file mountainvalley.tar.xz need to be decompressed and restored as a .msh file. We show the configuration files in the Accuracy verification test, and the other tests are similar. It contains the following files:

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

- _src_: this folder contains all sources file (including .cpp and .h).

- _figures_: this folder contains the figures presented in the article and corresponds to the data computed in the provided examples and tests.


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
hongbo.yao@outlook.com

* **Jingtian Tang**
School of Geosciences and Info-Physics, Central South University, Changsha, 410083, Hunan, China.
E-mail: jttang@csu.edu.cn




