
/*
 * @Description: This class implements a framework for the goal-oriented finite element method, including finite element discretization, assembly and solving of linear equation systems, a posteriori error estimation, and mesh refinement. 
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 * Input: Resistivity parameters, configuration settings, and a mesh file.
 * Output: Solution to the finite element linear equation systems.
 */





// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 


#ifndef _GOAFEM_H
#define _GOAFEM_H

#include "em.h"
#include "mfem.hpp"
#include "parahandler.h"

#include <mpi.h>
#include <vector>

using namespace std;
using namespace mfem;

class GOAFEM {
public:
  // MPI variables
  int myid;
  int n_procs;
  MPI_Comm comm;

  std::vector<DenseMatrix> cond_att;
  std::vector<DenseMatrix> sigma0;

  ParaHandler *para_handler;
  // Parallel mesh
  ParMesh *pmesh;

  FiniteElementCollection *pfec;

  // Parallel finite-element space
  ParFiniteElementSpace *pfes;
  BilinearFormIntegrator *integ; 

  ParBilinearForm *a;
  // Primal problem a(u,v) = f(v)
  std::vector<ParGridFunction *> up;
  // Dual problem a(w,v) = l(v)
  ParGridFunction *w;
  std::vector<ParLinearForm *> f;
  // dual problem linear form
  ParLinearForm *l;
  // Create "marker arrays" to define the portions of boundary associated
  //    with each type of boundary condition. These arrays have an entry
  //    corresponding to each boundary attribute.  Placing a '1' in entry i
  //    marks attribute i+1 as being active, '0' is inactive.
  Array<int> gamma0_bdr;
  Array<int> gamma1_bdr;

  MatrixCoefficient *LVolume_coef;
  Coefficient *LGamma1_coef;
  std::vector<VectorCoefficient *> RVolume_coef;
  std::vector<Coefficient *> RGamma0_coef;
  std::vector<Coefficient *> RGamma1_coef;
  OperatorHandle A;
  std::vector<Vector> U, F;
  Vector W, L;
  Vector local_err;
  Array<int> local_sources_tets;

  Array<int>
      global_initial_sources_tets; 

  Array<int> local_dual_tets; 

  int order; 

public:
  GOAFEM(std::vector<DenseMatrix> &att_cond, ParaHandler &para_handler_,
         ParMesh *pmesh_);
  ~GOAFEM();

public:
  void initialize();

  /*   Setup coefficient for left volume integral
    Setup coefficient for left distant surface integral
    Setup coefficient for right volume integral
    Setup coefficient for right distant surface integral
    Setup coefficient for air-earth surface integral */
  void setup_integral_coefficient();
  void set_bdr();

  // Problem size
  HYPRE_Int get_problem_dofs();
  void print_problem_size();
  void assemble_dual_linearform();
  void solve_with_pcg(OperatorHandle &A, Vector &X, Vector &B, std::string pre,
                      int maxit, std::string problem_type, int print_level,
                      int pre_print_level);

  // solve Boundary Value Problem with different solvers
  void solve(OperatorHandle &A, std::vector<Vector> &X, std::vector<Vector> &B,
             std::string problem_type);

  void solve_primal_problem();

  void solve_dual_problem();

  // Update all objects based on pmesh
  void update();
  void set_sigma0();
  // Estimating error
  void error_estimating();
  void refine_mesh();
};
#endif // _GOAFEM_H
