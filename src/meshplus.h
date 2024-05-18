
/*
 * @Description: This class inherits MFEM's class in mesh/mesh.hpp. I do this for extending the print function. This file is part of the MFEM library. 
 * @Author: Lewen liu, Zhengguang liu, Hongbo Yao and Jingtian Tang.
 * @Date: 2023-12-19 
 */


// Copyright (c) 2023.
// This file is part of the DC3DPAFEM program. DC3DPAFEM is free software with source code available in https://github.com/luowen765/DC3DPAFEM. You can redistribute it or modify it under the terms of the BSD-3 license. See file LICENSE for details. 



#ifndef _MESHPLUS_H
#define _MESHPLUS_H

#include "mfem.hpp"
using namespace mfem;

class Meshplus : public ParMesh {

public:
  Meshplus();
  Meshplus(ParMesh *pmesh_);
  void PrintVTK_eta(std::ostream &out, Vector &other_attri);
};

#endif // _MESHPLUS_H