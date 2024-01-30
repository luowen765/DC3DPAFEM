
#include "solversplus.h"
#include <iomanip>
#include <iostream>

using namespace std;

SolversPlus::SolversPlus() : CGSolver() {}
SolversPlus::SolversPlus(MPI_Comm comm_) : CGSolver(comm_) {}
void SolversPlus::Mult(const Vector &b, Vector &x) {
  int i;
  double r0, den, nom, nom0, betanom, alpha, beta;
  double initial_r;

  x.UseDevice(true);
  if (iterative_mode) {
    oper->Mult(x, r);
    subtract(b, r, r);
  } else {
    r = b;
    x = 0.0;
  }

  if (prec) {
    prec->Mult(r, z);
    d = z;
  } else {
    d = r;
  }
  nom0 = nom = Dot(d, r);

  MFEM_ASSERT(IsFinite(nom), "nom = " << nom);
  if (print_level == 1 || print_level == 3) {
    initial_r = sqrt(nom);
    mfem::out << "   Iteration : " << setw(3) << 0 << "  ||r|| = \t"
              << sqrt(nom) / initial_r << (print_level == 3 ? " ...\n" : "\n");
  }
  Monitor(0, nom, r, x);

  if (nom < 0.0) {
    if (print_level >= 0) {
      mfem::out
          << "PCG: The preconditioner is not positive definite. (Br, r) = "
          << nom << '\n';
    }
    converged = 0;
    final_iter = 0;
    final_norm = nom;
    return;
  }
  r0 = std::max(nom * rel_tol * rel_tol, abs_tol * abs_tol);

  if (nom <= r0) {
    converged = 1;
    final_iter = 0;
    final_norm = sqrt(nom);
    return;
  }

  oper->Mult(d, z);
  den = Dot(z, d);
  MFEM_ASSERT(IsFinite(den), "den = " << den);
  if (den <= 0.0) {
    if (Dot(d, d) > 0.0 && print_level >= 0) {
      mfem::out << "PCG: The operator is not positive definite. (Ad, d) = "
                << den << '\n';
    }
    if (den == 0.0) {
      converged = 0;
      final_iter = 0;
      final_norm = sqrt(nom);
      return;
    }
  }

  // start iteration
  converged = 0;
  final_iter = max_iter;
  for (i = 1; true;) {
    alpha = nom / den;
    add(x, alpha, d, x);
    add(r, -alpha, z, r);
    if (prec) {
      prec->Mult(r, z);
      betanom = Dot(r, z);
    } else {
      betanom = Dot(r, r);
    }
    MFEM_ASSERT(IsFinite(betanom), "betanom = " << betanom);
    if (betanom < 0.0) {
      if (print_level >= 0) {
        mfem::out
            << "PCG: The preconditioner is not positive definite. (Br, r) = "
            << betanom << '\n';
      }
      converged = 0;
      final_iter = i;
      break;
    }

    if (print_level == 1) {
      mfem::out << "   Iteration : " << setw(3) << i << "  ||r|| = \t"
                << sqrt(betanom) / initial_r << '\n';
    }

    Monitor(i, betanom, r, x);

    if (betanom <= r0) {
      if (print_level == 2) {
        mfem::out << "Number of PCG iterations: " << i << '\n';
      } else if (print_level == 3) {
        mfem::out << "   Iteration : " << setw(3) << i
                  << "  (B r, r) = " << betanom << '\n';
      }
      converged = 1;
      final_iter = i;
      break;
    }

    if (++i > max_iter) {
      break;
    }

    beta = betanom / nom;
    if (prec) {
      add(z, beta, d, d);
    } else {
      add(r, beta, d, d);
    }
    oper->Mult(d, z);
    den = Dot(d, z);
    MFEM_ASSERT(IsFinite(den), "den = " << den);
    if (den <= 0.0) {
      if (Dot(d, d) > 0.0 && print_level >= 0) {
        mfem::out << "PCG: The operator is not positive definite. (Ad, d) = "
                  << den << '\n';
      }
      if (den == 0.0) {
        final_iter = i;
        break;
      }
    }
    nom = betanom;
  }
  if (print_level >= 0 && !converged) {
    if (print_level != 1) {
      if (print_level != 3) {
        mfem::out << "   Iteration : " << setw(3) << 0
                  << "  (B r, r) = " << nom0 << " ...\n";
      }

      mfem::out << "   Iteration : " << setw(3) << final_iter
                << "  (B r, r) = " << betanom << '\n';
    }
    mfem::out << "PCG: No convergence!" << '\n';
  }
  if (print_level >= 1 || (print_level >= 0 && !converged)) {
    mfem::out << "Average reduction factor = "
              << pow(betanom / nom0, 0.5 / final_iter) << '\n';
  }
  final_norm = sqrt(betanom);

  Monitor(final_iter, final_norm, r, x, true);
}
