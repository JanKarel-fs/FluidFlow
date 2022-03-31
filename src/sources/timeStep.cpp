#include "timeStep.hpp"

double timeStep(const CellField<Compressible>& w, const Grid& g, const Settings& setting) {
  double dt = std::numeric_limits<double>::max();

#pragma omp parallel default(shared) reduction(min:dt)
  {
#pragma omp for
    for (int i=0; i<w.M(); i++) {
      for (int j=0; j<w.N(); j++) {
	Vector2d u = w[i][j].rhoU / w[i][j].rho;

	Point2d A = g.faceJ(i, j).center;
	Point2d B = g.faceJ(i+1, j).center;
	Point2d C = g.faceI(i, j).center;
	Point2d D = g.faceI(i, j+1).center;
	
	Vector2d s1(A, B);
	Vector2d s2(C, D);
	double dx = s1.length();
	double dy = s2.length();

	s1 /= dx;
	s2 /= dy;

	double lambda = 0.;

	switch (setting.convection) {
	case 0: break;
	case 1: {
	  double uTilde = dot(u, s1);
	  double vTilde = dot(u, s2);
	  double c = w[i][j].a();
	  double lambdaConv = (fabs(uTilde) + c) / dx + (fabs(vTilde) + c) / dy;
	  lambda += lambdaConv;
	}
	  break;
	default:
	  cout << "No such possibility for a computation of a convective part of a time step!" << endl;
	  exit(10);
	}

	switch (setting.diffusion) {
	case 0:break;
	case 1: {
	  double mu = w[i][j].mu();
	  double rho = w[i][j].rho;
	  double lambdaDiff = 2. * mu/rho * (1. / (dx*dx) + 1. / (dy*dy));
	  lambda += lambdaDiff;
	}
	  break;
	default:
	  cout << "No such possibility for a computation of a dissipative part of a time step!" << endl;
	  exit(10);
	}	

	double dt_ij = setting.CFL / lambda;


	if (dt_ij < dt) dt = dt_ij;
      }
    }
  }
  return dt;
}
