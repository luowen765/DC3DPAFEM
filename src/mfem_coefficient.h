


/*
 * @Description: The following classes are the direct copies and modifies of MFEM's classes in fem/coefficient.hpp. I just did a minor modifications to compute different coeffients for assembling equations. These classes are called by class GOAFEM in goafem.h.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 


#ifndef MFEM_COEFFICIENTS_H
#define MFEM_COEFFICIENTS_H

#include "em.h"
#include "mfem.hpp"
using namespace mfem;

class LVCoefficient : public MatrixCoefficient {
private:
  std::vector<DenseMatrix> mat;
  int temp_count;

public:
  /// Construct a coefficeint using matrix @a m for left side volume integral
  LVCoefficient(std::vector<DenseMatrix> &m)
      : MatrixCoefficient(m[0].Height(), m[0].Width()), mat(m) {}
  using MatrixCoefficient::Eval;
  /// Evaluate the matrix coefficient at @a ip.
  virtual void Eval(DenseMatrix &M, ElementTransformation &T,
                    const IntegrationPoint &ip);
};

class LGamma1Coefficient : public Coefficient {
private:
  std::vector<DenseMatrix> mat;
  Vector x_c;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a double coefficient for the surface integral in distant
  /// boundary face
  LGamma1Coefficient(std::vector<DenseMatrix> &m, Vector sources_center,
                     ParMesh *pmesh_)
      : mat(m), x_c(sources_center), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class RVCoefficient : public VectorCoefficient {
private:
  std::vector<DenseMatrix> mat;
  DenseMatrix sigma0;
  Vector x_i;
  // ParMesh *pmesh;

public:
  /// Construct a coeffient  using vector @a V for right side volume integral
  RVCoefficient(std::vector<DenseMatrix> &m, DenseMatrix sigma0_,
                Vector source_i)
      : VectorCoefficient(3), mat(m), sigma0(sigma0_), x_i(source_i)
  // ,pmesh(pmesh_)
  {}
  using VectorCoefficient::Eval;
  /// Evaluate the vector coefficient at @a ip.
  virtual void Eval(Vector &V, ElementTransformation &T,
                    const IntegrationPoint &ip);
};

class RGamma0Coefficient : public Coefficient {
private:
  DenseMatrix sigma0;
  Vector x_i;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a coeffient for right side volume integral
  RGamma0Coefficient(DenseMatrix sigma0_, Vector source_i, ParMesh *pmesh_)
      : sigma0(sigma0_), x_i(source_i), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the double coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

class RGamma1Coefficient : public Coefficient {
private:
  std::vector<DenseMatrix> mat;
  DenseMatrix sigma0;
  Vector x_i;
  ParMesh *pmesh;
  bool boundaryFlag = true;

public:
  /// Construct a coeffient for right side volume integral
  RGamma1Coefficient(std::vector<DenseMatrix> &m, DenseMatrix sigma0_,
                     Vector source_i, ParMesh *pmesh_)
      : mat(m), sigma0(sigma0_), x_i(source_i), pmesh(pmesh_) {}
  using Coefficient::Eval;
  /// Evaluate the double coefficient at @a ip.
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

#endif