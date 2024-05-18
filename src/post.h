
/*
 * @Description: This class is used for post-processing the solution of FEM, converting the potential into apparent resistivity.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 * Input: configuration settings, FEM parameters, and a mesh file.
 * Output: Apparent resistivity.
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 


#ifndef _POST_H
#define _POST_H

#include <fstream>
#include <string>
#include <vector>
#include "em.h"
#include "meshplus.h"
#include "mfem.hpp"
#include "parahandler.h"

using namespace mfem;

class Post {
public:
  ParMesh *pmesh;
  ParFiniteElementSpace *pfes;
  std::vector<ParGridFunction *> up;
  ParaHandler *para_handler;
  std::vector<Array<int>> local_sites_tets; 
  // potential Up(i,m) of source i at measuring site m owned by my rank
  std::vector<Array<double>> local_Up; 
  std::vector<Array<double>> Ex;
  std::vector<Array<double>> Ey;
  std::vector<Array<double>> Ez;
  std::vector<Array<double>> Jx;
  std::vector<Array<double>> Jy;
  std::vector<Array<double>> Jz;
  std::vector<Array<double>> Ut;

public:
  Post(ParMesh *pmesh_, ParFiniteElementSpace *pfes_,
       std::vector<ParGridFunction *> &up_, ParaHandler *para_handler_
  );

  ~Post();

public:
  // Post-processing and save the solution
  void main_post_process(int rank, ParaHandler *para_handler, std::string amr,
                         std::vector<DenseMatrix> &sigma0,
                         std::vector<DenseMatrix> &cond_att);

  // Post-processing
  void post_process(std::vector<DenseMatrix> &cond_att,std::vector<DenseMatrix> &sigma0);
  // Save results of global sites in rank 0 (root rank)
  void save_as_one(std::vector<DenseMatrix> &sigma0);
  void main_save_local_mesh(std::string amr_, int rank, Meshplus *pmesh,
                            Vector &local_eta,
                            std::vector<DenseMatrix> &att_cond);
  // compute total potential Ut(c,m) and apparent resistivity rho, print them
  void read_compute_survery_mode(std::string survery_mode_file,
                                 std::string output_file);
};

#endif // _POST_H
