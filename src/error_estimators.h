
/*
 * @Description: This class implements a posteriori error estimation based on continuity of current density (Ren et al., 2018b). Ren, Z., QIU, L., Tang, J., ZHOU, F., CHEN, C., CHEN, H., HU, S., 2018b. 3D modeling of direct-current anisotropic resistivity using the adaptive finite-element method based on continuity of current density. Chinese Journal of Geophysics 61, 331–343.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 * Input: resistivity parameters, configuration parameters, mesh file, and the finite element solution.
 * Output: the elemental posteriori error.
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 



#ifndef _ERROR_ESTIMATORS_H
#define _ERROR_ESTIMATORS_H

#include <cstdlib>

#include "mfem.hpp"
#include "parahandler.h"

using namespace mfem;

// Base error estimator class
class Estimator {
protected:
  // class ParaHandler for handlering geo-em model input parameter info
  ParaHandler *para_handler;
  ParMesh *pmesh;
  // parallel H(curl) finite-element space
  ParFiniteElementSpace *pfes;
  std::vector<ParGridFunction *> up;
  ParGridFunction *w;

public:
  Estimator(ParaHandler *para_handler_, ParMesh *pmesh_,
            ParFiniteElementSpace *pfes_, std::vector<ParGridFunction *> &up_,
            ParGridFunction *w_);
  virtual ~Estimator();

public:
  virtual void get_error_estimates(Vector &errors) = 0;
};

// Error estimator based on the face-jumps of the normal component of the
// current density nJ, sequential version is proposed by Zhengyong Ren et al.2018
// Here we implement a parallel version and call it nJEstimator.
class nJEstimator : public Estimator {
public:
  nJEstimator(ParaHandler *para_handler_, ParMesh *pmesh_,
              ParFiniteElementSpace *pfes_, std::vector<ParGridFunction *> &up_,
              ParGridFunction *w_, std::vector<DenseMatrix> &att_cond);
  ~nJEstimator();

public:
  std::vector<DenseMatrix> cond_att;
  double n_dot_J(bool if_loc1, bool if_shared_tet2,
                 FaceElementTransformations *face_elem_trans,
                 ElementTransformation *face_trans,
                 ElementTransformation *tet_trans, const FiniteElement *tet_fe,
                 int tet_id, const IntegrationPoint &ip, ParGridFunction *u);

  void compute_face_err(bool if_shared_tet2,
                        FaceElementTransformations *face_elem_trans,
                        double *face_err_u, double &face_err_w);

  void get_error_estimates(Vector &errors);
};

#endif // _ERROR_ESTIMATORS_H