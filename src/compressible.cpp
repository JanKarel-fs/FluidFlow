#include "compressible.hpp"

double Compressible::kappa;
double Compressible::cp;
double Compressible::cv;
double Compressible::R;
double Compressible::Pr;

double Compressible::p() const {
  double pressure = (kappa - 1.) * (e - 0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / rho);
  return pressure;
}

double Compressible::a() const {
  double speedOfSound = std::sqrt(kappa * p() / rho);
  return speedOfSound;
}

double Compressible::Ma() const {
  double MachNumber = rhoU.length() / rho / a();
  return MachNumber;
}

Compressible (*Compressible::flux)(const Compressible& wl, const Compressible& wr, const Vector2d& s);


Compressible Compressible::Upwind(const Compressible& wl, const Compressible& wr, const Vector2d& s) {
  // s - normalovy vektor prenasobeny delkou steny
  // numericky tok -> Upwind

  Vector2d n = s / s.length();
  Vector2d u = 0.5 * (wl.rhoU/wl.rho + wr.rhoU/wr.rho);
  double un = dot(u, n);

  Compressible flx;
  if (un >= 0) flx = wl * un;
  else flx = wr * un;

  double p = 0.5 * (wl.p() + wr.p());

  flx += Compressible(0., p*n, p*un);

  return flx * s.length();
}

Compressible Compressible::Rusanov(const Compressible& wl, const Compressible& wr, const Vector2d& s) {
  // Vector2d n = s / s.length();

  // double unL = dot(wl.rhoU/wl.rho, n);
  // double unR = dot(wr.rhoU/wr.rho, n);

  // double lambdaL = std::fabs(unL) + wl.a();
  // double lambdaR = std::fabs(unR) + wr.a();

  // double Slambda = std::max(lambdaL, lambdaR);

  // Compressible Fl(wl.rho*unL, wl.rhoU*unL + wl.p()*n, (wl.e + wl.p()) * unL);
  // Compressible Fr(wr.rho*unR, wr.rhoU*unR + wr.p()*n, (wr.e + wr.p()) * unR);

  // Compressible flx = 0.5 * (Fl + Fr) - 0.5 * Slambda * (wr - wl);

  // return flx * s.length();

  
  // Standardne: 1)rotace do n,t, 2)vypcoteni toku, 3) zpetna rotace do x,y
  Vector2d n = s / s.length();
  Vector2d t(-n.y, n.x);

  Compressible WR = wr;
  Compressible WL = wl;

  double rhouR = dot(wr.rhoU, n);
  double rhovR = dot(wr.rhoU, t);
  WR.rhoU = Vector2d(rhouR, rhovR);

  double rhouL = dot(wl.rhoU, n);
  double rhovL = dot(wl.rhoU, t);
  WL.rhoU = Vector2d(rhouL, rhovL);

  double lambdaL = std::fabs(WL.rhoU.x/WL.rho) + WL.a();
  double lambdaR = std::fabs(WR.rhoU.x/WR.rho) + WR.a();

  double Slambda = std::max(lambdaL, lambdaR);

  Compressible Fl(WL.rhoU.x, Vector2d(pow(WL.rhoU.x,2)/WL.rho + WL.p(), WL.rhoU.x*WL.rhoU.y/WL.rho),
		  (WL.e + WL.p()) * WL.rhoU.x/WL.rho);
  Compressible Fr(WR.rhoU.x, Vector2d(pow(WR.rhoU.x,2)/WR.rho + WR.p(), WR.rhoU.x*WR.rhoU.y/WR.rho),
		  (WR.e + WR.p()) * WR.rhoU.x/WR.rho);

  Compressible flx = 0.5 * (Fl + Fr) - 0.5 * Slambda * (WR - WL);

  Vector2d nInv(n.x, -n.y);
  Vector2d tInv(-t.x, t.y);

  double RHOU = dot(flx.rhoU, nInv);
  double RHOV = dot(flx.rhoU, tInv);

  flx.rhoU = Vector2d(RHOU, RHOV);

  return flx * s.length();
  
}


Compressible Compressible::fabs(const Compressible& w) {
  double RHO = std::fabs(w.rho);
  double RHOU = std::fabs(w.rhoU.x);
  double RHOV = std::fabs(w.rhoU.y);
  double E = std::fabs(w.e);

  return Compressible(RHO, Vector2d(RHOU, RHOV), E);
}

Compressible Compressible::sqrt(const Compressible& w) {
  double RHO = std::sqrt(w.rho);
  double RHOU = std::sqrt(w.rhoU.x);
  double RHOV = std::sqrt(w.rhoU.y);
  double E = std::sqrt(w.e);

  return Compressible(RHO, Vector2d(RHOU, RHOV), E);
}

Compressible Compressible::max(const Compressible& a, const Compressible& b) {
  double RHO = std::max(a.rho, b.rho);
  double RHOU = std::max(a.rhoU.x, b.rhoU.x);
  double RHOV = std::max(a.rhoU.y, b.rhoU.y);
  double E = std::max(a.e, b.e);

  return Compressible(RHO, Vector2d(RHOU, RHOV), E);
}

Compressible Compressible::min(const Compressible& a, const Compressible& b) {
  double RHO = std::min(a.rho, b.rho);
  double RHOU = std::min(a.rhoU.x, b.rhoU.x);
  double RHOV = std::min(a.rhoU.y, b.rhoU.y);
  double E = std::min(a.e, b.e);

  return Compressible(RHO, Vector2d(RHOU, RHOV), E);
}