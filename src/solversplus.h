
/*
 * @Description: This class inherits MFEM's class in linalg/solvers.h. I do this for modify the print part of CGSolver::Mult function. 
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */

// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 


#ifndef _SOLVERSPLUS_h
#define _SOLVERSPLUS_h

#include "mfem.hpp"
using namespace mfem;

class SolversPlus : public CGSolver {

public:
  SolversPlus();
  SolversPlus(MPI_Comm comm_);
  void Mult(const Vector &b, Vector &x);
};

#endif // _SOLVERSPLUS_h