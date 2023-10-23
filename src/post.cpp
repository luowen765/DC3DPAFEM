
#include "post.h"
#include "em.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>
#include <unistd.h>

Post::Post(ParMesh *pmesh_, ParFiniteElementSpace *pfes_,
           std::vector<ParGridFunction *> &up_, ParaHandler *para_handler_)
    : pmesh(pmesh_), pfes(pfes_), up(up_), para_handler(para_handler_) {}

Post::~Post() {}

// Post-processing and save the solution
void Post::main_post_process(int rank, ParaHandler *para_handler,
                             std::string amr, std::vector<DenseMatrix> &sigma0,
                             std::vector<DenseMatrix> &cond_att) {

  double tic = MPI_Wtime();
  local_sites_tets.resize(para_handler->source_number);

  for (int i = 0; i < para_handler->source_number; i++) {
    // Find elements which contain sites
    para_handler->find_point_tets(pmesh, para_handler->sites[i],
                                  local_sites_tets[i],
                                  para_handler->find_points_by);
  }
  // compute Up of measurements
  this->post_process(cond_att, sigma0);
  // compute Ut
  this->save_as_one(sigma0);
  // compute and output rho
  if (rank == 0) {
    std::string output_file = "solutions/solution" + std::string(".") + amr;
    this->read_compute_survery_mode(para_handler->sorted_survey, output_file);
  }

  double toc = MPI_Wtime();
  if (rank == 0) {
    std::cout << "The time cost of post_process is: " << (toc - tic) << " s\n";
  }
}

void Post::post_process(std::vector<DenseMatrix> &cond_att,
                        std::vector<DenseMatrix> &sigma0) {

  MPI_Comm mycomm = pmesh->GetComm();
  int myid, n_procs;
  MPI_Comm_rank(mycomm, &myid);
  int Ns = para_handler->source_number;
  local_Up.resize(Ns);
  if (para_handler->if_compute_E_J == "true") {
    if (myid == 0) {
      Ex.resize(Ns);
      Ey.resize(Ns);
      Ez.resize(Ns);
      Jx.resize(Ns);
      Jy.resize(Ns);
      Jz.resize(Ns);
    } else {
      std::cout << "current density J cannot support computing in parallel!"
                << std::endl;
      std::abort();
    }
  }

  for (int k = 0; k < para_handler->source_number; k++) {
    int nsite = para_handler->sites[k].size();
    myassert(sigma0.size() == para_handler->source_number);
    local_Up[k].SetSize(nsite);
    local_Up[k] = 0.0;
    if (para_handler->if_compute_E_J == "true") {
      Ex[k].SetSize(nsite);
      Ey[k].SetSize(nsite);
      Ez[k].SetSize(nsite);
      Jx[k].SetSize(nsite);
      Jy[k].SetSize(nsite);
      Jz[k].SetSize(nsite);
    }

    for (int i = 0; i < nsite; i++) {
      int tet_id = local_sites_tets[k][i];
      if (tet_id != -2) { // find in this process

        const FiniteElement *fe = pfes->GetFE(tet_id);
        int ndof = fe->GetDof();
        Vector pt(3);
        pt = 0.0;
        pt[0] = para_handler->sites[k][i](0);
        pt[1] = para_handler->sites[k][i](1);
        pt[2] = para_handler->sites[k][i](2);
        IntegrationPoint ip;
        ElementTransformation *tr = pfes->GetElementTransformation(tet_id);
        tr->TransformBack(pt, ip);
        tr->SetIntPoint(&ip);
        // Compute shape function
        Vector shape(ndof);
        shape = 0.0;
        fe->CalcPhysShape(*tr, shape);
        Array<int> vdofs;
        pfes->GetElementDofs(tet_id, vdofs);
        Vector tempUp;
        tempUp = 0.0;
        up[k]->GetSubVector(vdofs, tempUp);

        double potential = 0.0;
        for (int j = 0; j < ndof; j++) {
          potential += shape[j] * tempUp[j];
        }
        local_Up[k][i] = potential; // local Up

        if (para_handler->if_compute_E_J == "true") {
          // compute E and J
          DenseMatrix Dshape(ndof, 3);
          Dshape = 0.0;
          fe->CalcPhysDShape(*tr, Dshape);
          Vector tempE(3);
          tempE = 0.0;
          Dshape.AddMultTranspose(tempUp, tempE);
          Vector deltaUs(3);
          deltaUs = 0.0;
          Vector pi(3);
          pi = 0.0;
          for (int kk = 0; kk < 3; kk++) {
            pi[kk] = para_handler->sources[k](kk);
          }
          if (equalVector(pt, pi)) {
            std::cout << "pt: " << std::endl;
            pt.Print();
            std::cout << "pi: " << std::endl;
            pi.Print();
          }
          myassert(!equalVector(pt, pi));
          gradient_U_i_s(pt, pi, sigma0[k], deltaUs);
          myassert(Ex[k].Size() != 0);
          Ex[k][i] = -1.0 * (deltaUs(0) + tempE(0));
          Ey[k][i] = -1.0 * (deltaUs(1) + tempE(1));
          Ez[k][i] = -1.0 * (deltaUs(2) + tempE(2));
          Vector totalE(3), totalJ(3);
          totalE = 0.0;
          totalJ = 0.0;
          totalE[0] = Ex[k][i];
          totalE[1] = Ey[k][i];
          totalE[2] = Ez[k][i];
          DenseMatrix tet_cond = cond_att[tr->Attribute - 1];
          tet_cond.AddMult(totalE, totalJ);
          Jx[k][i] = totalJ(0);
          Jy[k][i] = totalJ(1);
          Jz[k][i] = totalJ(2);
        }

      } else {
        myassert(local_Up[k][i] == 0.0);
      }
    }
  }
}

void Post::save_as_one(std::vector<DenseMatrix> &sigma0) {

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status status;
  MPI_Comm mycomm = pmesh->GetComm();
  int myid, n_procs;
  MPI_Comm_size(mycomm, &n_procs);
  MPI_Comm_rank(mycomm, &myid);

  int count = -1;
  double *msg_Up;
  int Ns = para_handler->source_number;
  // Up.resize(Ns);
  Ut.resize(Ns);

  for (int k = 0; k < Ns; k++) {
    if (myid == 0) {
      for (int p = 1; p < n_procs; p++) {
        MPI_Recv(&count, 1, MPI_INT, p, 101, mycomm, &status);
        msg_Up = new double[count];
        MPI_Recv(msg_Up, count, MPI_DOUBLE, p, 105, mycomm, &status);
        for (int i = 0; i < count; i++) {
          local_Up[k][i] += msg_Up[i];
        }
        delete msg_Up;
      }

      // compute total potential Ut(sources,measures)
      int sites_number = para_handler->sites[k].size();
      Ut[k].SetSize(sites_number);
      Vector pi(3);
      pi = 0.0;
      for (int kk = 0; kk < 3; kk++) {
        pi[kk] = para_handler->sources[k](kk);
      }
      for (int h = 0; h < sites_number; h++) {
        Vector pj(3);
        pj = 0.0;
        for (int gg = 0; gg < 3; gg++) {
          pj[gg] = para_handler->sites[k][h](gg);
        }
        if (pj(0) == pi(0) && pj(1) == pi(1) && pj(2) == pi(2)) {
          Ut[k][h] = 1E10;
        } else {
          Ut[k][h] = local_Up[k][h] + U_i_s(pj, pi, sigma0[k]);
        }
      }
    } else {
      count = local_Up[k].Size();
      msg_Up = local_Up[k].GetData();
      MPI_Send(&count, 1, MPI_INT, 0, 101, mycomm);
      MPI_Send(msg_Up, count, MPI_DOUBLE, 0, 105, mycomm);
    }
  } // source loop
}
// Save local mesh
void Post::main_save_local_mesh(std::string amr_, int rank, Meshplus *pmesh,
                                Vector &local_eta,
                                std::vector<DenseMatrix> &att_cond) {

  MPI_Barrier(MPI_COMM_WORLD);
  std::ostringstream vtk_name2;
  vtk_name2 << "solutions/"
            << "err" << std::setfill('0') << std::setw(2) << rank << "." << amr_
            << ".vtk";
  std::ofstream vtk_ofs2(vtk_name2.str().c_str());
  vtk_ofs2.precision(8);
  pmesh->PrintVTK_eta(vtk_ofs2, local_eta);
}
// checked
void Post::read_compute_survery_mode(std::string survery_mode_file,
                                     std::string output_file) {

  unsigned int n_data = 0;
  unsigned int mode = 0;
  int source;

  std::ifstream mode_stream(survery_mode_file.c_str());
  myassert(mode_stream.good());

  std::ofstream out_stream(output_file.c_str());
  myassert(out_stream.good());

  mode_stream >> n_data >> mode;

  if (mode == 11) { // pole-pole configuration
    out_stream << n_data << "\t" << mode << "\n";
    int ns = para_handler->sources.size();
    int surveys = 0;
    for (int i = 0; i < ns; i++) {
      int nsite = para_handler->sites[i].size();
      surveys += nsite;
      for (int j = 0; j < nsite; j++) {
        int C_local_id = i;
        int P_local_id = j;

        Vertex C_node = para_handler->sources[C_local_id];
        Vertex P_node = para_handler->sites[C_local_id][P_local_id];
        double a = length_two_point(C_node.operator()(), P_node.operator()());
        double k = 2.0 * PI * a;
        double U = this->Ut[C_local_id][P_local_id];
        double rho = k * U;
        out_stream << C_local_id << "\t" << (C_node)(0) << "\t" << (C_node)(1)
                   << "\t" << (C_node)(2) << "\t" << P_local_id << "\t"
                   << (P_node)(0) << "\t" << (P_node)(1) << "\t" << (P_node)(2)
                   << "\t" << k << "\t" << U << "\t" << rho << "\n";
      }
    }
    myassert(n_data == surveys);

    if (para_handler->if_compute_E_J == "true") {
      std::string out_JU = "EJU.txt";
      // compute JU
      std::ofstream out(out_JU);
      myassert(out.good());
      out.precision(8);
      out << "x"
          << "\t"
          << "y"
          << "\t"
          << "z"
          << "\t"
          << "Jx"
          << "\t"
          << "Jy"
          << "\t"
          << "Jz"
          << "\t"
          << "U\n";
      int Ns = para_handler->source_number;
      for (int i = 0; i < Ns; i++) {
        int Nsite = para_handler->sites[i].size();
        for (int j = 0; j < Nsite; j++) {
          out << para_handler->sites[i][j](0) << "\t"
              << para_handler->sites[i][j](1) << "\t"
              << para_handler->sites[i][j](2) << "\t" << Jx[i][j] << "\t"
              << Jy[i][j] << "\t" << Jz[i][j] << "\t" << Ut[i][j] << std::endl;
        }
      }
    }
  } else {
    std::cout << "Other modes are not supported temporarily!";
  }

  mode_stream.close();
  out_stream.close();

  return;
}
