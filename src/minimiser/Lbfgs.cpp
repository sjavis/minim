#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "vec.h"


Lbfgs::Lbfgs(State &state, int m, int maxIter)
  : Minimiser(state, maxIter), _m(m), _g(state.ndof), _rho(m),
    _s(m, std::vector<double>(state.ndof)),
    _y(m, std::vector<double>(state.ndof))
{}


void Lbfgs::iteration() {
  std::vector<double> step(state.ndof);
  std::vector<double> g0(state.ndof);

  if (iter == 0) {
    g0 = state.gradient();
  } else {
    g0 = _g;
  }
  _i_cycle = iter % _m;

  getDirection(step);
  linesearch(step);
  _g = state.gradient();

  _s[_i_cycle] = step;
  _y[_i_cycle] = vec::diff(_g, g0);
  _rho[_i_cycle] = 1 / vec::dotProduct(step, _y[_i_cycle]);
}


void Lbfgs::getDirection(std::vector<double> &step) {
  double alpha[_m] = {0};
  int m_tmp = (_m < iter) ? _m : iter;

  _g = state.gradient();
  if (iter == 0) {
    step = vec::multiply(-_init_hessian, _g);
    return;
  } else {
    step = vec::multiply(-1, _g);
  }

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (_i_cycle - 1 - i1 + _m) % _m;
    alpha[i] = _rho[i] * vec::dotProduct(step, _s[i]);
    for (int j=0; j<state.ndof; j++) {
      step[j] -= alpha[i] * _y[i][j];
    }
  }

  int i = (_i_cycle - 1 + _m) % _m;
  double gamma = 1 / (_rho[i] * vec::dotProduct(_y[i], _y[i]));
  step = vec::multiply(gamma, step);

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (_i_cycle - m_tmp + i1 + _m) % _m;
    double beta = _rho[i] * vec::dotProduct(step, _y[i]);
    step = vec::sum(step, vec::multiply(alpha[i]-beta, _s[i]));
  }
}


void Lbfgs::linesearch(std::vector<double> &step) {
  const double c = 0.5; // Armijo control parameter
  const double tau = 0.5; // Shrink factor

  double t = - c * vec::dotProduct(_g, step);
  double e0 = state.energy();

  for (int i=0; i<10; i++) {
    double e = state.energy(vec::sum(state.coords, step));
    if (e0-e >= t) break;

    step = vec::multiply(step, tau);
    t = t * tau;
  }

  state.coords = vec::sum(state.coords, step);
}
