#include <cmath>
#include <vector>
#include "System.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "linalg.h"


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
  y[i_cycle] = linalg::diff(g, g0);
  rho[i_cycle] = 1 / linalg::dotProduct(step, y[i_cycle]);
}


bool Lbfgs::checkConvergence() {
  double sum = linalg::dotProduct(g, g);
  double rms = sqrt(sum/sys.ndof);
  return (rms < sys.convergence_rms);
}


void Lbfgs::getDirection(std::vector<double> &step) {
  double alpha[m] = {0};
  int m_tmp = (m < iter) ? m : iter;

  sys.gradient(step);
  step = linalg::multiply(-1, step);
  if (iter == 0) return;

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (i_cycle - 1 - i1 + m) % m;
    alpha[i] = rho[i] * linalg::dotProduct(step, s[i]);
    for (int j=0; j<sys.ndof; j++) {
      step[j] -= alpha[i] * y[i][j];
    }
  }

  int i = (i_cycle - 1 + m) % m;
  double gamma = 1 / (rho[i] * linalg::dotProduct(y[i], y[i]));
  step = linalg::multiply(gamma, step);

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (i_cycle - m_tmp + i1 + m) % m;
    double beta = rho[i] * linalg::dotProduct(step, y[i]);
    step = linalg::sum(step, linalg::multiply(alpha[i]-beta, s[i]));
  }
}


void Lbfgs::linesearch(std::vector<double> &step) {
  double e0 = sys.energy();

  for (int i=0; i<5; i++) {
    double e = sys.energy(linalg::sum(sys.state, step));
    if (e < e0) break;
    step = linalg::multiply(0.1, step);
  }

  sys.state = linalg::sum(sys.state, step);
}
