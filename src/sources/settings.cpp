#include "settings.hpp"

Settings::Settings(const string& name) {
  cout << "Reading data from: " << name << endl;

  set<string> sections;
  sections.insert("GRID"); sections.insert("INITIAL_CONDITIONS");
  sections.insert("BOUNDARY_CONDITIONS");sections.insert("FLUX_SPLITTER");
  sections.insert("TIME"); sections.insert("PHYSICAL_VALUES");
  sections.insert("SAVING");

  map<string, vector<string> > dataFile;
  loadDataFile(name, sections, dataFile);

  // Reading information about grid
  mCells = 0;
  nCells = 0;
  string section = "GRID";
  findSection(dataFile, "grid_type", section, grid_type);
  switch (grid_type) {
  case 1:
    findSection(dataFile, "mCells", section, mCells);
    findSection(dataFile, "nCells", section, nCells);
    break;
  default:
    cout << "Loading starter.txt: no a such grid type!" << endl;
    exit(19);
  }
  findSection(dataFile, "ghostCells", section, ghostCells);


  // Reading information about initial conditions
  section = "INITIAL_CONDITIONS";
  findSection(dataFile, "rhoInit", section, rhoInit);
  double u, v;
  findSection(dataFile, "uInit", section, u);
  findSection(dataFile, "vInit", section, v);
  uInit = Vector2d(u, v);
  findSection(dataFile, "pInit", section, pInit);

  // Reading information about boundary conditions
  section = "BOUNDARY_CONDITIONS";
  int numOfBoundaries;
  findSection(dataFile, "numOfBoundaries", section, numOfBoundaries);

  for (int i=1; i<=numOfBoundaries; i++) {
    stringstream number;
    number << i;
    string numberString;
    number >> numberString;

    string boundary = "boundary" + numberString;
    string bcType = "bcType" + numberString;

    string boundaryName, bcTypeName;
    
    findSection(dataFile, boundary, section, boundaryName);
    findSection(dataFile, bcType, section, bcTypeName);

    usedBC[boundaryName] = bcTypeName;
  }
  
  findSection(dataFile, "alpha", section, alpha);
  findSection(dataFile, "Ma2is", section, Ma2is);
  
  // Reading information about flux splitter
  section = "FLUX_SPLITTER";
  findSection(dataFile, "flux", section, flux);

  // Reading information about time
  section = "TIME";
  findSection(dataFile, "CFL", section, CFL);

  // Reading information about physical values
  section = "PHYSICAL_VALUES";
  findSection(dataFile, "kappa", section, kappa);
  findSection(dataFile, "rho0", section, rho0);
  findSection(dataFile, "p0", section, p0);
  findSection(dataFile, "R", section, R);
  findSection(dataFile, "Pr", section, Pr);

  // Reading information about stop criteria
  section = "SAVING";
  findSection(dataFile, "stop", section, stop);
}
