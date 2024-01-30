
#include "mfem_coefficient.h"
#include "em.h"

#include <iostream>

void LVCoefficient::Eval(DenseMatrix &M, ElementTransformation &T,
                         const IntegrationPoint &ip) {
  M = mat[T.Attribute - 1];
}

double LGamma1Coefficient::Eval(ElementTransformation &T,
                                const IntegrationPoint &ip) {

  Vector normal(3);
  normal = 0.0;
  compute_normal(T, ip, normal);

  Vector x(3);
  x = 0.0;
  T.SetIntPoint(&ip);
  // Transform integration point from reference coordinates to physical
  // coordinates and store them in the vector
  T.Transform(ip, x); 
  int tet_id = -1;
  get_tetId_bdr(pmesh, T.ElementNo, tet_id);
  DenseMatrix sigma = mat[tet_id];
  double mass_coefficient = Beta_value(x, x_c, sigma, normal);
  return mass_coefficient;
}

void RVCoefficient::Eval(Vector &V, ElementTransformation &T,
                         const IntegrationPoint &ip) {
  T.SetIntPoint(&ip);
  Vector x(3);
  x = 0.0;
  // Transform integration point from reference coordinates to physical
  // coordinates and store them in the vector
  T.Transform(ip, x);

  Vector deltaUs(3);
  deltaUs = 0.0;
  gradient_U_i_s(x, x_i, sigma0, deltaUs);
  DenseMatrix temp = sigma0;
  temp -= mat[T.Attribute - 1];
  V.SetSize(3);
  V = 0.0;
  temp.AddMult(deltaUs, V);
}

double RGamma0Coefficient::Eval(ElementTransformation &T,
                                const IntegrationPoint &ip) {
  // compute outer noraml unit vector
  Vector normal(3); 
  normal = 0.0;
  compute_normal(T, ip, normal);

  Vector x(3);
  x = 0.0;
  T.SetIntPoint(&ip);
  // Transform integration point from reference coordinates to physical
  // coordinates and store them in the vector
  T.Transform(ip, x); 
  double Us = 0.0;
  Us = U_i_s(x, x_i, sigma0);
  double beta = Beta_value(x, x_i, sigma0, normal);
  double mass_coef = beta * Us;
  return mass_coef;
}

double RGamma1Coefficient::Eval(ElementTransformation &T,
                                const IntegrationPoint &ip) {

  //  std::cout<<T.Attribute<<" ";
  // compute outer noraml unit vector
  Vector normal(3); 
  normal = 0.0;
  compute_normal(T, ip, normal);

  Vector x(3);
  x = 0.0;
  T.SetIntPoint(&ip);
  // Transform integration point from reference coordinates to physical
  // coordinates and store them in the vector
  T.Transform(ip, x); 

  double Us = 0.0;
  Us = U_i_s(x, x_i, sigma0);

  int tet_id = -1;
  get_tetId_bdr(pmesh, T.ElementNo, tet_id);
  DenseMatrix sigma = mat[tet_id];
  double coef1 = Beta_value(x, x_i, sigma, normal) * Us;
  double coef2 = Beta_value(x, x_i, sigma0, normal) * Us;
  double mass_coef = coef2 - coef1;
  return mass_coef;
}

