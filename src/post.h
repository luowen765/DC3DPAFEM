
// Copyright (c) 2023.
// This file is part of the 3DDCAF program. 3DDCAF is free software, you can redistribute it and/or modify it under the terms of the BSD-3 license. See file LICENSE for details.

/*
 * @Description:
  This class does post-processing for the solutions of FEM. The total potential can be obtained by adding the primary singular potential. Then apparent resistivity can be computed on different configrations. For more
information and source code availability, please visit https://github.com/luowen765/3DDCAF.
 * @Author: Lewen liu; Zhengguang liu; Hongbo Yao.
 */

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
