#include <iostream>
#include <cstdlib>
#include "fvm/cellField.hpp"
#include "fvm/computeResidue.hpp"
#include "sources/setBoundaryConditions.hpp"
#include "sources/settings.hpp"
#include "sources/initialisation.hpp"
#include "sources/timeStep.hpp"
#include "sources/setGrid.hpp"
#include "sources/bcList.hpp"
#include "sources/typedefs.hpp"
#include "saving/storeNormResidue.hpp"
#include "saving/saveResults.hpp"
#include "sources/step.hpp"
#include "sources/linearSolver.hpp"
#include "compressible.hpp"

using namespace std;

int main(int argc,char **args) {

  PetscInitialize( &argc , &args , (char *)0 , 0 );
  
  Settings setting("starter.txt");
  
  Grid g;
  setGrid(g, setting);

  LinearSolver<Compressible> linSolver(setting.solver, g);

  map<string, bcWithJacobian> BC;

  for (auto it = setting.usedBC.begin(); it != setting.usedBC.end(); it++) {
    auto itbcL = bcList.find(it->second);
    if (itbcL != bcList.end()) {
      BC[it->first] = bcList[it->second];
    }
    else {
      cout << "No a such boundary condition type!" << endl;
      cout << "Possible boundary condition types: ";
      for (itbcL = bcList.begin(); itbcL != bcList.end(); itbcL++) {
	cout << itbcL->first << ", ";
      }
      cout << endl;
      exit(15);
    } 
  }
  
  CellField<Compressible> w(g), wOld(g), res(g);

  double dt;
  
  cout << "Start" << endl;

  initialisation(w, g, setting);
  wOld = w;

  for (int i=1; i<=setting.stop; i++) {
    dt = timeStep(w, g, setting);
    
    step<Compressible>(w, wOld, res, g, dt, BC, linSolver, setting);

    if (i%5 == 0) {
      cout << "iterace: " << i << ", dt = " << dt << endl;
      storeNormResidue(res, g, i);
    }
  }

  cout << "Saving results..." << endl;
  setBoundaryConditions(w, g, setting, BC);
  saveResults(w, g);

  cout << "Good bye!" << endl;

  linSolver.free();

  PetscFinalize();

  return 0;
}
