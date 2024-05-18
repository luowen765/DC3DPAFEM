/*
 @Description: This file is the main program of the package DC3DPAFEM, serving
 as the overall implementation framework for the algorithm. It includes reading
 input files, domain decomposition, solving finite element linear equation
 systems, and the mesh adaptation process. For more details, please refer to the
 article: _A scalable parallel finite-element algorithm using the algebraic
 multigrid solver for 3-D direct current resistivity modeling in anisotropic
 media_.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with
// source code available in https://github.com/luowen765/DC3DPAFEM. You can
// redistribute it or modify it under the terms of the BSD-3 license. See file
// LICENSE for details.

#include "em.h"
#include "goafem.h"
// #include "memwatch.h"
#include "meshplus.h"
#include "mfem.hpp"
#include "parahandler.h"

#include "post.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;

int main(int argc, char *argv[]) {

  // MPI handler
  int n_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  double main_tic, main_toc;
  main_tic = MPI_Wtime();

  // Check program usage
  if (argc < 2) {
    if (myid == 0) {
      std::cout << "Usage: " << argv[0] << " input_model_filename\n";
    }
    return 1;
  }

  extern std::vector<double> memoryUsed;
  extern double tempMemory;

  ParaHandler para_handler(argv[1]);

  std::string result_path = "solutions/";
  if (access(result_path.c_str(), 0)) {
    std::cout << "folder not exists, create it ..." << std::endl;
    std::ostringstream command;
    command.str("");
    command << "mkdir " + result_path;
    int status = std::system(command.str().c_str());
    // assert(status == 0); /* normal exit*/
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Mesh *mesh = new Mesh(para_handler.mesh_file.c_str(), 1, 1);
  para_handler.preProcess();

  if (myid == 0)
    para_handler
        .printForwardParameters(); // print model.config file at terminal
  // test for Scalability
  int uniN = 0;
  for (int i = 0; i < para_handler.N_uniformRefine; i++) {
    mesh->UniformRefinement();
    uniN++;
    std::ostringstream msh_name;
    msh_name << para_handler.mesh_file << "." << std::to_string(i + 1);
    if (myid == 0) {
      mesh->Save(msh_name.str().c_str());
    }
    std::cout << msh_name.str().c_str() << " DOFs: " << mesh->GetNV()
              << std::endl;
  } // UniformRefinement
  std::vector<DenseMatrix> cond_att(mesh->GetNE());
  for (int i = 0; i < mesh->GetNE(); i++) {
    cond_att[i] = para_handler.get_elem_conductivity(mesh->GetAttribute(i));
    mesh->SetAttribute(i, i + 1);
  }
  mesh->SetAttributes();
  // Generate parallel mesh by a partitioning a serial mesh using METIS
  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  GOAFEM *goafem = new GOAFEM(cond_att, para_handler, pmesh);
  // AMR loop
  for (int iter = 1; iter <= para_handler.maxit; iter++) { // iter start from 1
    std::string amr = std::to_string(iter);
    goafem->initialize();
    if (goafem->get_problem_dofs() > para_handler.max_dofs) {
      if (myid == 0) {
        std::cout << "\nStop due to reached max number of dofs\n";
        std::cout << "Present dofs: " << goafem->get_problem_dofs() << "\n\n";
      }
      break;
    }

    if (goafem->get_problem_dofs() > 2147483648) {
      if (myid == 0) {
        std::cout << "\nmax_dofs has exited the range of int-type data, please "
                     "modify your code!";
      }
      break;
    }
    if (myid == 0) {
      std::cout << "*********************** AMR Loop " << iter
                << "**********************\n";
    }
    // Print number of elements, complex unknows and dofs
    goafem->print_problem_size();
    tempMemory = GetCurrentMemoryUsage();
    memoryUsed.push_back(tempMemory);
    goafem->solve_primal_problem();

    if (myid == 0) {
      std::cout << "Post processing...\n";
    }

    Post post(pmesh, goafem->pfes, goafem->up, &para_handler);

    // save sites solutions of Ut and compute apparent resistivity ...
    if (amr == "1" && uniN != 0) {
      post.main_post_process(myid, &para_handler,
                             "uni" + std::to_string(uniN) + "-" +
                                 std::to_string(n_procs),
                             goafem->sigma0, goafem->cond_att);
    } else {
      post.main_post_process(myid, &para_handler, amr, goafem->sigma0,
                             goafem->cond_att);
    }
    Meshplus pmeshplus = Meshplus(pmesh); // used to print mesh result

    goafem->solve_dual_problem();

    goafem->error_estimating();
    Vector relative_eta = goafem->local_err;
    double local_max_err = relative_eta.Max();
    double global_max_err;
    MPI_Allreduce(&local_max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX,
                  pmesh->GetComm());
    relative_eta *= (1.0 / global_max_err);

    if (para_handler.save_amr_mesh == 1) {
      post.main_save_local_mesh(amr, myid, &pmeshplus, relative_eta, cond_att);

    } else {
      if (myid == 0)
        std::cout << "do not save any mesh!" << std::endl;
    }
    goafem->refine_mesh();
    goafem->update();

    if (iter == 1 && memoryUsed.size() != 0) {
      if (myid == 0) {
        double max_element_error =
            *std::max_element(memoryUsed.begin(), memoryUsed.end());
        std::cout << "myid: " << myid << " need max memory in this code: \t"
                  << max_element_error << "\t\n";
      }
      std::cout << "myid: " << myid << " need memory for solving PDE: \t"
                << memoryUsed[1] - memoryUsed[0] << "\t\n";
    }
  }
  delete goafem;
  delete pmesh;

  main_toc = MPI_Wtime();
  if (myid == 0) {
    double local_time;
    double max_time = main_toc - main_tic;
    for (int n = 1; n < n_procs; n++) {
      MPI_Recv(&local_time, 1, MPI_DOUBLE, n, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      if (max_time < local_time) {
        max_time = local_time;
      }
    }
    std::cout << "Process unit: " << n_procs << "\n";
    std::cout << "the total calculation Time is:\t" << max_time << " seconds\n";
  } else {
    double local_time = main_toc - main_tic;
    MPI_Send(&local_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();

  return 0;
}
