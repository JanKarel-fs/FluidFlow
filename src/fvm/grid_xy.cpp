#include "grid_xy.hpp"

void Grid_xy::loadXYFiles(const string& name1, const string& name2,
		   vector<double>& xCoord, vector<double>& yCoord) {

  ifstream fin(name1.c_str());
  if (!fin) {
    cout << "File \"" << name1 << "\" with x coordinates can not be opened!" << endl;
    exit(20);
  }

  int nb;
  fin >> nb;

  xCoord.resize(nb);
  for (int i=0; i<nb; i++) {
    fin >> xCoord[i];
  }

  fin.close();

  fin.open(name2.c_str());
  
  if (!fin) {
    cout << "File \"" << name2 << "\" with y coordinates can not be opened!" << endl;
    exit(20);
  }

  fin >> nb;

  yCoord.resize(nb);
  for (int i=0; i<nb; i++) {
    fin >> yCoord[i];
  }

  fin.close();
}

Grid_xy::Grid_xy(const string name1, const string name2, const int gh, const string& type) {
  vector<double> xCoord, yCoord;

  loadXYFiles(name1, name2, xCoord, yCoord);
    
  //m, n pocet bunek ve smeru i, j
  Mnodes = xCoord.size();
  Nnodes = yCoord.size();
  Mvolumes = Mnodes - 1;
  Nvolumes = Nnodes -1;
  ghost = gh;

  nodes.allocate(-ghost, Mnodes+ghost, -ghost, Nnodes+ghost);
  centers.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  volumes.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  facesI.allocate(-ghost, Mvolumes+ghost, -ghost, Nnodes+ghost);
  facesJ.allocate(-ghost, Mnodes+ghost, -ghost, Nvolumes+ghost);

  //generovani souradnic
  
  for (int i=0; i<Mnodes; i++) {
    for (int j=0; j<Nnodes; j++) {
      nodes[i][j].vertex = Point2d(xCoord[i], yCoord[j]);
    }
  }

  // vytvoreni pomocnych bunek
  for (int j=0; j<Nnodes; j++) {
    for (int k=1; k<=ghost; k++) {
      //leva hranice
      Vector2d s(nodes[1][j].vertex, nodes[0][j].vertex);
       nodes[-k][j].vertex = nodes[0][j].vertex + k * s;

      //prava hranice
      s = Vector2d(nodes[Mnodes-2][j].vertex, nodes[Mnodes-1][j].vertex);
      nodes[Mnodes-1+k][j].vertex = nodes[Mnodes-1][j].vertex + k * s;
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int k=1; k<=ghost; k++) {
      //dolni hranice
      Vector2d s(nodes[i][1].vertex, nodes[i][0].vertex);
      nodes[i][-k].vertex = nodes[i][0].vertex + k * s;

      //horni hranice
      s = Vector2d(nodes[i][Nnodes-2].vertex, nodes[i][Nnodes-1].vertex);
      nodes[i][Nnodes-1+k].vertex = nodes[i][Nnodes-1].vertex + k * s;
    }
  }

  //vypocet objemu bunek a stredu bunek
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      Vector2d u1(nodes[i][j].vertex, nodes[i+1][j+1].vertex);
      Vector2d u2(nodes[i+1][j].vertex, nodes[i][j+1].vertex);

      volumes[i][j] = 0.5 * fabs(cross(u1, u2));
      centers[i][j] = 0.25 * (node(i, j) + node(i+1, j) + node(i+1, j+1) + node(i, j+1));
    }
  }

  map<string, coefficients>::iterator it;

  it = mCoefficients.find(type);

  if (it != mCoefficients.end()) {
    computeCoefficients = it->second;
  }
  else {
    cout << "There is no possibility \"" << type << "\" for evaluating of coefficients in verteces!" << endl;
    cout<< "Possible choise is: ";
    for (it=mCoefficients.begin(); it!=mCoefficients.end(); it++) {
      cout << it->first << ", ";
    }
    cout << endl;
    exit(11);
  }

  computeCoefficients(*this);
  
  update_faces();

  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nnodes+ghost; j++) {
      facesI[i][j].name = "internal";

      if (j == 0) {
	if (i < Mnodes/2) {
	  facesI[i][j].name = "lower";
	}
	else {
	  facesI[i][j].name = "wall";
	}
      }
      
      if (j == (Nnodes-1)) {
	facesI[i][j].name = "upper";
      }
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      facesJ[i][j].name = "internal";

      if (i == 0) facesJ[i][j].name = "inlet";
      if (i == (Mnodes-1)) facesJ[i][j].name = "outlet";
    }
  }
}
