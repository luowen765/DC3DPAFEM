
// Copyright (c) 2023.
// This file is part of the 3DDCAF program. 3DDCAF is free software, you can redistribute it and/or modify it under the terms of the BSD-3 license. See file LICENSE for details.

/*
 * @Description:
This class inherits MFEM's class in mesh/mesh.hpp. I do this for extending
the print function This file is part of the MFEM library. For more
information and source code availability, please visit https://github.com/luowen765/3DDCAF.
 * @Author: Lewen liu; Zhengguang liu; Hongbo Yao.
 */

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