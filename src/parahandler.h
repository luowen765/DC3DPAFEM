
/*
 * @Description: This class is used to read various parameter settings from a configuration file.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details.  


#ifndef _PARAHANDLER_H
#define _PARAHANDLER_H

#include "em.h"
#include "mfem.hpp"
#include <fstream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

using namespace EM;
using namespace mfem;
class ParaHandler {
public:

  int source_number; 
  std::vector<Vertex> sources;
  std::vector<std::vector<Vertex>> sites;
  std::vector<Vertex> s_plus_m; 
  std::string survey_input;
  std::string sorted_survey;
  std::string model_parameters_file;
  int n_regions;          
  std::vector<int> marker_vec; 
  std::map<int, DenseMatrix> region2conductivity; 

  // FEM parameters
  int N_uniformRefine;     // the number for mesh UniformRefine function
  int maxit;               // max AMR iterations
  long int max_dofs;       // max dofs
  double beta;             // marking parameter
  // Mesh parameters
  int save_amr_mesh;
  std::string mesh_file;
  int pcg_maxit;
  double pcg_primal_tol;
  double pcg_dual_tol;
  int pcg_print_level;
  int amg_print_level;
  std::string if_compute_E_J;

public:
  ParaHandler(char *model_file);
  ~ParaHandler();

public:

  void skip_comments(std::istream &in, std::vector<std::string> &para_str_vec,
                     std::string comment_str = "#");
  // read model information from file
  void read_model_info(char *model_file);
  // calculate the tensor conductivity
  DenseMatrix cal_conductivity(Vector &main_cond, Vector &_angle);
  // get element conductivity according to marker (attribute)
  DenseMatrix get_elem_conductivity(int marker);
  // get center coordinates of sources
  Vector get_sources_center();

  void preProcess();
  void find_point_tets(ParMesh *pmesh, std::vector<Vertex> &points,
                       Array<int> &find_tets, std::string find_points_by);
  void printForwardParameters();
};
#endif // _PARAHANDLER_H
