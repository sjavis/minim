#include <cmath>
#include <vector>
#include "System.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "vec.h"
#include <iostream>


Lbfgs::Lbfgs(System &sys, int m_, int maxIter)
: Minimiser(sys, maxIter), m(m_), g(sys.ndof), rho(m),
  s(m, std::vector<double>(sys.ndof)),
  y(m, std::vector<double>(sys.ndof))
{}


void Lbfgs::iteration() {
  std::vector<double> step(sys.ndof);
  std::vector<double> g0(sys.ndof);

  if (iter == 0) {
    sys.gradient(g0);
  } else {
    g0 = g;
  }
  i_cycle = iter % m;

  getDirection(step);
  linesearch(step);
  sys.gradient(g);

  s[i_cycle] = step;
  y[i_cycle] = vec::diff(g, g0);
  rho[i_cycle] = 1 / vec::dotProduct(step, y[i_cycle]);
}


bool Lbfgs::checkConvergence() {
  double sum = vec::dotProduct(g, g);
  double rms = sqrt(sum/sys.ndof);
  return (rms < sys.convergence_rms);
}


void Lbfgs::getDirection(std::vector<double> &step) {
  double alpha[m] = {0};
  int m_tmp = (m < iter) ? m : iter;

  sys.gradient(g);
  if (iter == 0) {
    step = vec::multiply(-init_hessian, g);
    return;
  } else {
    step = vec::multiply(-1, g);
  }

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (i_cycle - 1 - i1 + m) % m;
    alpha[i] = rho[i] * vec::dotProduct(step, s[i]);
    for (int j=0; j<sys.ndof; j++) {
      step[j] -= alpha[i] * y[i][j];
    }
  }

  int i = (i_cycle - 1 + m) % m;
  double gamma = 1 / (rho[i] * vec::dotProduct(y[i], y[i]));
  step = vec::multiply(gamma, step);

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (i_cycle - m_tmp + i1 + m) % m;
    double beta = rho[i] * vec::dotProduct(step, y[i]);
    step = vec::sum(step, vec::multiply(alpha[i]-beta, s[i]));
  }
}


void Lbfgs::linesearch(std::vector<double> &step) {
  const double c = 0.5; // Armijo control parameter
  const double tau = 0.5; // Shrink factor

  double t = - c * vec::dotProduct(g, step);
  double e0 = sys.energy();

  for (int i=0; i<10; i++) {
    double e = sys.energy(vec::sum(sys.state, step));
    if (e0-e >= t) break;

    step = vec::multiply(step, tau);
    t = t * tau;
  }

  sys.state = vec::sum(sys.state, step);
}
