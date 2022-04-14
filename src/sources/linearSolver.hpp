#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <iostream>
#include <cstdlib>
#include "fvm/grid.hpp"
#include "fvm/cellField.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

using namespace std;

template <typename var>
class LinearSolver {
  bool activated;

public:
  LinearSolver(): activated(false) {};
  LinearSolver(const int& type, const Grid& g) {

    nmb = var::nVars * g.Mvol() * g.Nvol();
    switch (type) {
    case 1: activated = false;
      break;
    case 2: {
      activated = true;
      
      VecCreate(PETSC_COMM_WORLD, &b);
      VecSetSizes(b, PETSC_DECIDE, nmb);
      VecSetFromOptions(b);
      VecSet(b, 0.);

      VecDuplicate(b, &x);
      VecSet(x, 0.);

      MatCreate(PETSC_COMM_WORLD, &A);
      MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nmb, nmb);
      MatSetFromOptions(A);
      MatSetUp(A);
      MatZeroEntries(A);

      KSPCreate(PETSC_COMM_WORLD, &solver);
      KSPSetOperators(solver, A, A);

      KSPSetFromOptions(solver);
      KSPSetUp(solver);
    }
      break;
    default:
      cout << "No such possibility for a linear solver!" << endl;
      exit(31);
    }
    
  };
  ~LinearSolver() {
    if (activated) {
      KSPDestroy(&solver);
      MatDestroy(&A);
      VecDestroy(&x);
      VecDestroy(&b);
    }
  };

  int nmb;
  
  Vec x, b;
  Mat A;
  KSP solver;

  void solve() {
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    
    KSPSetOperators(solver, A, A);
    KSPSolve(solver, b, x);
  };

  void getResults(CellField<var>& res) {
    const PetscScalar *px;

    VecGetArrayRead(x, &px);

    for (int i=0; i<res.M(); i++) {
      for (int j=0; j<res.N(); j++) {
	for (int k=0; k<var::nVars; k++) {
	  res[i][j][k] = px[var::nVars * (j*res.M() + i) + k];
	}
      }
    }

    VecRestoreArrayRead(x, &px);
  };

  void reset() {
    MatZeroEntries(A);
    VecSet(b, 0.);
  };
};

#endif
