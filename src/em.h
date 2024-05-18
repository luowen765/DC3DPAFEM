/*
 * @Description: The namespace of the package DC3DPAFEM includes a variety of public functions and variables.
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */


// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 




// 工具函数 和公共变量

#ifndef _EM_H
#define _EM_H

#include <complex>
#include <fstream>
#include <iostream>
#include <random>

#include <string>
#include <unistd.h> 
#include <sys/resource.h>
#include "mfem.hpp"
using namespace mfem;

namespace EM {
static const double PI = 3.1415926535897932384626433832795;
#define MAX_INT std::numeric_limits<int>::max()
// Tolerance used in our package.
static const double TOLERANCE = 1e-6;
} // namespace EM
using namespace EM;


inline void throwError(const std::string &errString) {
  std::cerr << errString << std::endl;
}

template <class T, class U> T max(const T &a, const U &b) {
  return std::max(a, T(b));
}

void print_time(double time_in_seconds);

// Calculate the length between two points: coord1 and coord2
double length_two_point(double *coord1, double *coord2);
double length_two_point(Vector &coord1, Vector &coord2);

// sorted by x，y，z
bool cmp(std::vector<double> a, std::vector<double> b);
// sorted by x1,y1,z1 x2,y2,z2
bool cmp2(std::vector<double> a, std::vector<double> b);

void remove_duplicate_vertex(std::vector<mfem::Vertex> &all_point);

double B_value(Vector p, Vector pi, DenseMatrix sigma);

double Beta_value(Vector p, Vector pi, DenseMatrix sigma, Vector normal);

double U_i_s(Vector p, Vector pi, DenseMatrix sigma0); 

void gradient_U_i_s(mfem::Vector p, mfem::Vector pi, mfem::DenseMatrix sigma0,
                    mfem::Vector &deltaUs); 

void compute_normal(ElementTransformation &T, const IntegrationPoint &ip,
                    Vector &normal);

void get_tetId_bdr(ParMesh *pmesh, int ElementNo, int &tet_id);

double GetCurrentMemoryUsage();

// whether a==b, note width == height for a and b
bool equalDM(DenseMatrix a, DenseMatrix b);
bool equalVector(Vector a, Vector b);

void Mult_DenseMatrix3(DenseMatrix &a, DenseMatrix &b, DenseMatrix &c);

#endif // _EM_H
