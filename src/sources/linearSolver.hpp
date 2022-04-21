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
      VecSet(x, 1.);

      MatCreate(PETSC_COMM_WORLD, &A);
      MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nmb, nmb);
      MatSetType(A, MATAIJ);
      MatSeqAIJSetPreallocation(A, 20, NULL);
      MatMPIAIJSetPreallocation(A, 20, NULL, 20, NULL);
      MatSetFromOptions(A);
      MatSetUp(A);

      KSPCreate(PETSC_COMM_WORLD, &solver);
      //KSPSetType(solver, KSPMINRES);

      // PC prec;
      // KSPGetPC(solver,&prec);
      // PCSetType(prec,PCJACOBI);
      // //KSPSetOperators(solver, A, A);

      // KSPSetFromOptions(solver);
      // KSPSetUp(solver);
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
  KSPConvergedReason reason;
  PetscInt  its;
  PetscReal rnorm;
  
  void solve() {
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    
    KSPSetOperators(solver, A, A);
    KSPSolve(solver, b, x);
      
    PetscInt  restart, maxits;
    PetscReal rtol, abstol;
    KSPGetTolerances(solver, &rtol, &abstol, NULL, &maxits);
    KSPGMRESGetRestart(solver, &restart);
      
    cout<<"Initial settings "<<endl;
    cout<<"Relative tolerance : "<<rtol<<endl;
    cout<<"Absolute tolerance : "<<abstol<<endl;
    cout<<"Maximal number of iterations : "<<maxits<<endl;
    cout<<"Number of Krylov vectors before restarting : "<<restart<<endl<<endl;
      
    cout<<"Convergence reason     KSP iterations      Residual norm    "<<endl;

    KSPGetConvergedReason(solver, &reason);
    KSPGetIterationNumber(solver, &its);
    KSPGetResidualNorm(solver, &rnorm);
    
        
    cout << reason << " " << its << " " << rnorm << " " << endl;
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

  void free() {
    if (activated) {
      KSPDestroy(&solver);
      MatDestroy(&A);
      VecDestroy(&x);
      VecDestroy(&b);
    }
    activated = false;
  };
};

#endif
