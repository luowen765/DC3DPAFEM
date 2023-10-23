

#include "goafem.h"
#include "em.h"
#include "error_estimators.h"
#include "mfem_coefficient.h"
#include "solversplus.h"

#include <algorithm>
#include <cassert>
#include <cmath>

std::vector<double> memoryUsed;
double tempMemory;

GOAFEM::GOAFEM(std::vector<DenseMatrix> &att_cond, ParaHandler &para_handler_,
               ParMesh *pmesh_)
    : cond_att(att_cond), myid(0), n_procs(1), para_handler(&para_handler_),
      pmesh(pmesh_), pfes(NULL), pfec(NULL), w(NULL), a(NULL), l(NULL),
      LVolume_coef(NULL), LGamma1_coef(NULL), integ(NULL) {
  int sn = para_handler->source_number;
  U.resize(sn);
  F.resize(sn);
  up.resize(sn);
  f.resize(sn);
  RVolume_coef.resize(sn), RGamma0_coef.resize(sn), RGamma1_coef.resize(sn);
  for (int i = 0; i < sn; i++) {
    up[i] = NULL;
    f[i] = NULL;
    RVolume_coef[i] = NULL;
    RGamma0_coef[i] = NULL;
    RGamma1_coef[i] = NULL;
  }
}

GOAFEM::~GOAFEM() {
  delete pfes;
  delete pfec;
  delete a;
  delete w;
  delete l;
  delete LVolume_coef;
  delete LGamma1_coef;
  for (int i = 0; i < para_handler->source_number; i++) {
    delete up[i];
    delete f[i];
    delete RVolume_coef[i];
    delete RGamma0_coef[i];
    delete RGamma1_coef[i];
  }
}

void GOAFEM::initialize() {
  // MPI variables
  comm = pmesh->GetComm();
  MPI_Comm_size(comm, &n_procs);
  MPI_Comm_rank(comm, &myid);

  local_sources_tets.DeleteAll();
  // find sources
  para_handler->find_point_tets(pmesh, para_handler->sources,
                                local_sources_tets,
                                para_handler->find_points_by);
  this->set_sigma0();
  // Parallel  finite-element space
  int dim = pmesh->Dimension(); // The dimention of model.
  order = 1;
  pfec = new H1_FECollection(order, dim);
  pfes = new ParFiniteElementSpace(pmesh, pfec);

  // Bilinear form a(u,v), a(w,v)
  a = new ParBilinearForm(pfes);
  for (int i = 0; i < para_handler->source_number; i++) {
    // Gridfunction of primary problem
    up[i] = new ParGridFunction(pfes);
    *(up[i]) = 0.0;
    f[i] = new ParLinearForm(pfes);
    f[i]->Vector::operator=(0.0);
  }
  // Gridfunction of dual problem
  w = new ParGridFunction(pfes);
  *w = 0.0;
  // Linear form of dual problem
  l = new ParLinearForm(pfes);
  l->Vector::operator=(0.0); // VERY IMPORTANT!
  // number of elements
  local_err.SetSize(pmesh->GetNE());
  local_err = 0.0;
}

void GOAFEM::set_bdr() {
  gamma0_bdr.SetSize(pmesh->bdr_attributes.Size());
  gamma1_bdr.SetSize(pmesh->bdr_attributes.Size());
  gamma0_bdr = 0;
  gamma1_bdr = 0;
  gamma0_bdr[0] = 1;
  gamma1_bdr[1] = 1;
}

void GOAFEM::setup_integral_coefficient() {
  LVolume_coef = new LVCoefficient(cond_att);
  Vector sources_center = para_handler->get_sources_center();
  LGamma1_coef = new LGamma1Coefficient(cond_att, sources_center, pmesh);

  // Setup coefficient for right volume integral
  for (int i = 0; i < para_handler->source_number; i++) {
    Vector source_i(3);
    source_i = 0.0;
    for (int j = 0; j < 3; j++) {
      source_i[j] = para_handler->sources[i](j);
    }
    RVolume_coef[i] = new RVCoefficient(cond_att, sigma0[i], source_i);
    RGamma0_coef[i] = new RGamma0Coefficient(sigma0[i], source_i, pmesh);
    RGamma1_coef[i] =
        new RGamma1Coefficient(cond_att, sigma0[i], source_i, pmesh);
  }
}

HYPRE_Int GOAFEM::get_problem_dofs() { return pfes->GlobalTrueVSize(); }
void GOAFEM::print_problem_size() {
  HYPRE_Int size = pfes->GlobalTrueVSize();
  long n_global_elems = pmesh->GetGlobalNE();
  if (myid == 0) {
    std::cout << "Number of tetrahedral elements: " << n_global_elems << "\n";
    std::cout << "Degrees of Freedom (DOFs): " << size << "\n";
  }
}

void GOAFEM::solve_primal_problem() {
  double start, finish;
  start = MPI_Wtime();

  if (myid == 0)
    std::cout << "start assemble primal problem\n";
  /*   Setup coefficient for left volume integral
    Setup coefficient for left distant surface integral
    Setup coefficient for right volume integral
    Setup coefficient for right distant surface integral
    Setup coefficient for air-earth surface integral */
  this->setup_integral_coefficient();
  this->set_bdr();
  integ = new DiffusionIntegrator(*LVolume_coef);
  a->AddDomainIntegrator(integ);
  a->AddBoundaryIntegrator(new MassIntegrator(*LGamma1_coef), gamma1_bdr);
  a->Assemble();
  a->Finalize();

  for (int i = 0; i < para_handler->source_number; i++) {
    f[i]->AddDomainIntegrator(new DomainLFGradIntegrator(*(RVolume_coef[i])));
    f[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(*(RGamma0_coef[i])),
                                gamma0_bdr);
    f[i]->AddBoundaryIntegrator(new BoundaryLFIntegrator(*(RGamma1_coef[i])),
                                gamma1_bdr);
    f[i]->Assemble();
    Array<int> ess_tdof_list(0);
    // Linear system of primal problem au=f
    a->FormLinearSystem(ess_tdof_list, *(up[i]), *(f[i]), A, U[i], F[i]);
  }

  double finish_assemble_PDE = MPI_Wtime();

  if (myid == 0)
    std::cout << "The time of assembling PDE is " << finish_assemble_PDE - start
              << " s\n";

  double start_solve_primalPDE, finish_solve_primalPDE;
  start_solve_primalPDE = MPI_Wtime();

  std::string problem_type = "primal";

  if (myid == 0)
    std::cout << "start solve primal problem\n";
  this->solve(A, U, F, problem_type);

  finish_solve_primalPDE = MPI_Wtime();
  if (myid == 0)
    std::cout << "The time of solving PDE is "
              << finish_solve_primalPDE - start_solve_primalPDE << " s\n";

  for (int i = 0; i < para_handler->source_number; i++) {
    a->RecoverFEMSolution(U[i], *(f[i]), *(up[i]));
  }

  finish = MPI_Wtime();
  if (myid == 0)
    std::cout << "The elapsed time of solving primal problem is "
              << finish - start << " s\n";
}

void GOAFEM::solve_dual_problem() {
  double start, finish;
  start = MPI_Wtime();

  if (myid == 0)
    std::cout << "start assemble dual problem "
              << " \n";
  local_dual_tets.DeleteAll();
  para_handler->find_point_tets(pmesh, para_handler->s_plus_m, local_dual_tets,
                                para_handler->find_points_by);
  assemble_dual_linearform();
  Array<int> ess_tdof_list(0);
  myassert(a != NULL);
  // Linear system of dual problem aw=l
  a->FormLinearSystem(ess_tdof_list, *w, *l, A, W, L);

  std::string problem_type = "dual";
  std::vector<Vector> W_, L_;
  W_.push_back(W);
  L_.push_back(L);
  this->solve(A, W_, L_, problem_type);
  a->RecoverFEMSolution(W_[0], *l, *w);
  finish = MPI_Wtime();
  if (myid == 0)
    std::cout << "The elapsed time of solving dual problem is "
              << finish - start << " s\n";
}

void GOAFEM::assemble_dual_linearform() {
  // Compute non-zero element vector and assemble
  for (int i = 0; i < local_dual_tets.Size(); i++) {
    int tet_id = local_dual_tets[i];
    if (tet_id != -2) { // find it
      const FiniteElement *fe = pfes->GetFE(tet_id);
      const IntegrationRule *ir =
          &IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder() + 3);
      ElementTransformation *eltrans = pfes->GetElementTransformation(tet_id);
      double elem_volume = pmesh->GetElementVolume(tet_id);
      const int ndof = fe->GetDof();
      Vector shape(ndof);
      shape = 0.0;
      Vector elemvec(ndof);
      elemvec = 0.0;

      for (int q = 0; q < ir->GetNPoints(); ++q) {
        const IntegrationPoint &ip = ir->IntPoint(q);
        eltrans->SetIntPoint(&ip);
        double WJ = eltrans->Weight() * ip.weight;
        // fe->CalcPhysVShape(*eltrans, shape);
        fe->CalcPhysShape(*eltrans, shape);
        for (int k = 0; k < ndof; k++) // loop all dofs
        {
          elemvec[k] += WJ * (1.0 / elem_volume) * shape[k];
        }
      } // q loop
      Array<int> vdofs;
      pfes->GetElementDofs(tet_id, vdofs);
      l->AddElementVector(vdofs, elemvec);
    } // find it
  }
}

void GOAFEM::solve_with_pcg(OperatorHandle &A, Vector &X, Vector &B,
                            std::string pre, int maxit,
                            std::string problem_type, int print_level,
                            int pre_print_level) {
  double tol = 0.0;
  if (problem_type == "primal") {
    tol = para_handler->pcg_primal_tol;
  } else if (problem_type == "dual") {
    tol = para_handler->pcg_dual_tol;
  } else {
    std::cout << "please check your problem_type! primal or dual?" << std::endl;
    std::abort();
  }

  Solver *prec = NULL;

  if (pre == "amg") {
    HypreBoomerAMG *amg = new HypreBoomerAMG;
    amg->SetPrintLevel(pre_print_level);
    prec = amg;
  } else {
    std::cout << "Unsupported input linear solver: "
              << para_handler->linear_solver << "\n";
    std::abort();
  }

  SolversPlus cg(MPI_COMM_WORLD);
  cg.SetMaxIter(maxit);
  cg.SetRelTol(tol);
  cg.SetPrintLevel(print_level);
  // only
  if (prec) {
    cg.SetPreconditioner(*prec);
  }
  cg.SetOperator(*A);
  tempMemory = GetCurrentMemoryUsage();
  memoryUsed.push_back(tempMemory);

  cg.Mult(B, X);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == 0)
    std::cout << problem_type + " iterations: " << cg.GetNumIterations()
              << std::endl;
  delete prec;
}

void GOAFEM::solve(OperatorHandle &A, std::vector<Vector> &X,
                   std::vector<Vector> &B, std::string problem_type) {
  std::string linear_solver = para_handler->linear_solver;
  for (int i = 0; i < B.size(); i++) {
    solve_with_pcg(A, X[i], B[i], linear_solver, para_handler->pcg_maxit,
                   problem_type, para_handler->pcg_print_level,
                   para_handler->amg_print_level);
  }
}

void GOAFEM::update() {
  pfes->Update();
  for (int i = 0; i < para_handler->source_number; i++) {
    up[i]->Update();
  }
  w->Update();

  a->Update();
  for (int i = 0; i < para_handler->source_number; i++) {
    f[i]->Update();
  }
  l->Update();
}

void GOAFEM::set_sigma0() {
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status status;
  MPI_Comm mycomm = pmesh->GetComm();
  int myid, n_procs;
  MPI_Comm_size(mycomm, &n_procs);
  MPI_Comm_rank(mycomm, &myid);

  int count = -1;
  int *msg_tets;

  int ns = local_sources_tets.Size();
  myassert(ns == para_handler->source_number);
  global_initial_sources_tets.SetSize(ns);

  if (myid == 0) {
    for (int i = 0; i < ns; i++) {
      int tet_id = local_sources_tets[i];
      if (tet_id != -2) {
        global_initial_sources_tets[i] = pmesh->GetAttribute(tet_id) - 1;
      } else {
        global_initial_sources_tets[i] = 0;
      }
    }

    for (int p = 1; p < n_procs; p++) {
      MPI_Recv(&count, 1, MPI_INT, p, 101, mycomm, &status);
      msg_tets = new int[count];
      MPI_Recv(msg_tets, count, MPI_INT, p, 105, mycomm, &status);

      // Append the data of other threads
      for (int i = 0; i < count; i++) {
        global_initial_sources_tets[i] += msg_tets[i];
      }
      delete msg_tets;
    }
  } else {
    count = local_sources_tets.Size();
    msg_tets = new int[count];
    for (int i = 0; i < count; i++) {
      int tet_id = local_sources_tets[i];
      if (tet_id != -2) {
        msg_tets[i] = pmesh->GetAttribute(tet_id) - 1;
      } else {
        msg_tets[i] = 0;
      }
    }
    MPI_Send(&count, 1, MPI_INT, 0, 101, mycomm);
    MPI_Send(msg_tets, count, MPI_INT, 0, 105, mycomm);
    delete msg_tets;
  }

  int *sources_tets_g;
  sources_tets_g = new int[ns];

  if (myid == 0) {
    for (int j = 0; j < ns; j++) {
      sources_tets_g[j] = global_initial_sources_tets[j];
    }
    for (int p = 1; p < n_procs; p++) {
      MPI_Send(sources_tets_g, ns, MPI_INT, p, 6, mycomm);
    }
  } else {
    MPI_Recv(sources_tets_g, ns, MPI_INT, 0, 6, mycomm, &status);
    global_initial_sources_tets.SetSize(ns);
    for (int j = 0; j < ns; j++) {
      global_initial_sources_tets[j] = sources_tets_g[j];
    }
  }
  delete sources_tets_g;

  sigma0.clear();
  sigma0.resize(para_handler->source_number);

  for (int k = 0; k < para_handler->source_number; k++) {
    sigma0[k] = cond_att[global_initial_sources_tets[k]];
  }
}

void GOAFEM::error_estimating(int iter) {
  // Error estiamating
  if (para_handler->error_esti_way == "jn") {
    nJEstimator nJ_estimator(para_handler, pmesh, pfes, up, w, cond_att);
    nJ_estimator.get_error_estimates(local_err);
  } else {
    std::cout << "unsupported estimator Temporarily!" << std::endl;
  }
}

void GOAFEM::refine_mesh() {

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status status;
  MPI_Comm mycomm = pmesh->GetComm();
  int myid, n_procs;
  MPI_Comm_size(mycomm, &n_procs);
  MPI_Comm_rank(mycomm, &myid);

  double threshold;
  double local_max_err = local_err.Max();
  double global_max_err;
  MPI_Allreduce(&local_max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX,
                pmesh->GetComm());
  threshold = para_handler->beta * global_max_err;
  pmesh->RefineByError(local_err, threshold);
}
