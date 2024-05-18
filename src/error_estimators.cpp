

/*
 * @Description: This .cpp file contains the actual implementations of some functions declared in the error_estimators.h header file.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */


// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 


#include "error_estimators.h"
#include "em.h"

// Base error estimator class
Estimator::Estimator(ParaHandler *para_handler_, ParMesh *pmesh_,
                     ParFiniteElementSpace *pfes_,
                     std::vector<ParGridFunction *> &up_, ParGridFunction *w_)
    : para_handler(para_handler_), pmesh(pmesh_), pfes(pfes_), up(up_), w(w_) {}

Estimator::~Estimator() {}

// Error estimator based on the face-jumps of the normal component of the
// current density nJ, sequential version is proposed by Zhengyong Ren et al.
// Here we implement a parallel version and call it nJEstimator.
nJEstimator::nJEstimator(ParaHandler *para_handler_, ParMesh *pmesh_,
                         ParFiniteElementSpace *pfes_,
                         std::vector<ParGridFunction *> &up_,
                         ParGridFunction *w_,
                         std::vector<DenseMatrix> &att_cond)
    : Estimator(para_handler_, pmesh_, pfes_, up_, w_) {
  cond_att = att_cond;
}

nJEstimator::~nJEstimator() {}

double nJEstimator::n_dot_J(bool if_loc1, bool if_shared_tet2,
                            FaceElementTransformations *face_elem_trans,
                            ElementTransformation *face_trans,
                            ElementTransformation *tet_trans,
                            const FiniteElement *tet_fe, int tet_id,
                            const IntegrationPoint &ip, ParGridFunction *u) {
  int ndof = tet_fe->GetDof();
  int dim = tet_fe->GetDim();
  assert(dim == 3);

  // unit normal with Jacobian
  Vector normalJ(dim);
  normalJ = 0.0;
  CalcOrtho(face_trans->Jacobian(), normalJ);
  // unit normal
  Vector normal(dim);
  normal = 0.0;
  normal.Set(1.0 / normalJ.Norml2(), normalJ);
  IntegrationPoint tet_ip;
  if (if_loc1) {
    face_elem_trans->Loc1.Transform(ip, tet_ip);
  } else {
    face_elem_trans->Loc2.Transform(ip, tet_ip);
  }
  tet_trans->SetIntPoint(&tet_ip);
  DenseMatrix tet_shape(ndof, dim);              // 4*3
  tet_fe->CalcPhysDShape(*tet_trans, tet_shape); // delta U
  Array<int> tet_vdofs;
  Vector tU;
  if (!if_shared_tet2) {
    pfes->GetElementVDofs(tet_id, tet_vdofs);
    u->GetSubVector(tet_vdofs, tU);
  } else {
    pfes->GetFaceNbrElementVDofs(tet_id, tet_vdofs);
    u->FaceNbrData().GetSubVector(tet_vdofs, tU);
  }

  assert(tet_vdofs.Size() == ndof);
  assert(tU.Size() == ndof);
  Vector tempE(3);
  tempE = 0.0;
  for (int j = 0; j < ndof; j++) {
    for (int c = 0; c < 3; c++) {
      double shape = tet_shape(j, c);
      tempE(c) += tU[j] * shape;
    }
  }
  Vector tempJe(3);
  tempJe = 0.0;
  DenseMatrix tet_cond = cond_att[tet_trans->Attribute - 1];
  tet_cond.AddMult(tempE, tempJe);

  double nJ = 0.0;
  nJ = normal * tempJe;
  return nJ;
}

void nJEstimator::compute_face_err(bool if_shared_tet2,
                                   FaceElementTransformations *face_elem_trans,
                                   double *face_err_u, double &face_err_w) {
  ElementTransformation *face_trans = face_elem_trans->Face;
  int tet1_id;
  const FiniteElement *tet1_fe;
  ElementTransformation *tet1_trans;
  int tet2_id;
  const FiniteElement *tet2_fe;
  ElementTransformation *tet2_trans;

  if (!if_shared_tet2) {
    tet1_id = face_elem_trans->Elem1No;
    tet1_fe = pfes->GetFE(tet1_id);
    tet1_trans = face_elem_trans->Elem1;
    tet2_id = face_elem_trans->Elem2No;
    tet2_fe = pfes->GetFE(tet2_id);
    tet2_trans = face_elem_trans->Elem2;
  } else {
    tet1_id = face_elem_trans->Elem1No; // local element
    tet1_fe = pfes->GetFE(tet1_id);
    tet1_trans = face_elem_trans->Elem1;
    tet2_id = face_elem_trans->Elem2No - pmesh->GetNE(); // BE CAREFUL!
    tet2_fe = pfes->GetFaceNbrFE(tet2_id);               // BE CAREFUL!
    tet2_trans = face_elem_trans->Elem2;
  }
  const IntegrationRule *face_ir =
      &IntRules.Get(face_trans->GetGeometryType(), 2 * face_trans->Order() + 3);

  for (int q = 0; q < face_ir->GetNPoints(); ++q) {
    const IntegrationPoint &ip = face_ir->IntPoint(q);
    face_trans->SetIntPoint(&ip);
    double WJ = face_trans->Weight() * ip.weight;

    // calculate nJu and nJw in face's tet1
    double nJ_tet1_u = 0.0;
    double nJ_tet1_w = 0.0;
    // calculate nJu and nJw in face's tet2
    double nJ_tet2_u = 0.0;
    double nJ_tet2_w = 0.0;

    for (int s = 0; s < para_handler->source_number; s++) {

      // if_loc1 = true; if_shared_tet2 = false;
      nJ_tet1_u = this->n_dot_J(true, false, face_elem_trans, face_trans,
                                tet1_trans, tet1_fe, tet1_id, ip, up[s]);

      // if_loc1 = false;
      nJ_tet2_u =
          this->n_dot_J(false, if_shared_tet2, face_elem_trans, face_trans,
                        tet2_trans, tet2_fe, tet2_id, ip, up[s]);

      double nJu = std::abs(nJ_tet1_u) - std::abs(nJ_tet2_u);
      face_err_u[s] += 0.5 * nJu * nJu * WJ;
    }

    // if_loc1 = true; if_shared_tet2 = false;
    nJ_tet1_w = this->n_dot_J(true, false, face_elem_trans, face_trans,
                              tet1_trans, tet1_fe, tet1_id, ip, w);

    // if_loc1 = false;
    nJ_tet2_w = this->n_dot_J(false, if_shared_tet2, face_elem_trans,
                              face_trans, tet2_trans, tet2_fe, tet2_id, ip, w);

    double nJw = std::abs(nJ_tet1_w) - std::abs(nJ_tet2_w);
    face_err_w += 0.5 * nJw * nJw * WJ;
  } // integral point of  face done
}

void nJEstimator::get_error_estimates(Vector &errors) {

  for (int s = 0; s < para_handler->source_number; s++) {
    up[s]->ExchangeFaceNbrData();
  }
  w->ExchangeFaceNbrData();

  int n_faces = pmesh->GetNumFaces(); // 3D, return pmesh->GetNFaces()
  for (int f = 0; f < n_faces; f++) {
    // Interior faces (without boundary faces and shared faces)
    FaceElementTransformations *face_elem_trans =
        pmesh->GetInteriorFaceTransformations(f);
    if (face_elem_trans == NULL) // boundary faces or MPI shared faces
    {
      continue;
    }

    int tet1_id = face_elem_trans->Elem1No;
    int tet2_id = face_elem_trans->Elem2No;

    double face_err_w = 0.0;
    double face_err_u[para_handler->source_number];
    for (int i = 0; i < para_handler->source_number; i++) {
      face_err_u[i] = 0.0;
    }

    bool if_shared_tet2 = false;
    compute_face_err(if_shared_tet2, face_elem_trans, face_err_u, face_err_w);
    double elem_error_face = 0.0;
    for (int s = 0; s < para_handler->source_number; s++) {
      elem_error_face += std::sqrt(face_err_u[s]) * std::sqrt(face_err_w);
    }
    errors(tet1_id) += elem_error_face;
    errors(tet2_id) += elem_error_face;
  } // loop all faces done

  for (int f = 0; f < pmesh->GetNSharedFaces(); f++) {
    // In the returned object, 1 and 2 refer to the local and the neighbor
    // elements, respectively.
    FaceElementTransformations *face_elem_trans =
        pmesh->GetSharedFaceTransformations(f);

    int tet1_id = face_elem_trans->Elem1No; // local element

    double face_err_w = 0.0;
    double face_err_u[para_handler->source_number];
    for (int i = 0; i < para_handler->source_number; i++) {
      face_err_u[i] = 0.0;
    }

    bool if_shared_tet2 = true;
    compute_face_err(if_shared_tet2, face_elem_trans, face_err_u, face_err_w);

    double elem_error_face = 0.0;

    // MPI process will take care of tet2_id
    for (int s = 0; s < para_handler->source_number; s++) {
      elem_error_face += std::sqrt(face_err_u[s]) * std::sqrt(face_err_w);
    }
    errors(tet1_id) += elem_error_face;
  } // loop all shared faces done
}
