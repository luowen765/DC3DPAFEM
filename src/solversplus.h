
// Copyright (c) 2023.
// This file is part of the 3DDCAF program. 3DDCAF is free software, you can redistribute it and/or modify it under the terms of the BSD-3 license. See file LICENSE for details.

/*
 * @Description:
  This class inherits MFEM's class in linalg/solvers.h. I do this for modify
  the print part of CGSolver::Mult function. For more
information and source code availability, please visit https://github.com/luowen765/3DDCAF.
 * @Author: Lewen liu; Zhengguang liu; Hongbo Yao.
 */

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