#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "linesearch.h"
#include "Communicator.h"
#include "vec.h"


Lbfgs::Lbfgs(State &state)
  : Minimiser(state), _m(5), _rho(_m),
    _s(_m, std::vector<double>(state.nblock)),
    _y(_m, std::vector<double>(state.nblock)),
    _g0(state.nblock), _g1(state.nblock), _step(state.nblock)
{}

Lbfgs::Lbfgs(State &state, AdjustFunc adjustModel)
  : Minimiser(state, adjustModel), _rho(_m),
    _s(_m, std::vector<double>(state.nblock)),
    _y(_m, std::vector<double>(state.nblock)),
    _g0(state.nblock), _g1(state.nblock), _step(state.nblock)
{}


Lbfgs& Lbfgs::setM(int m) {
  _m = m;
  _rho.resize(m);
  _s.resize(m, std::vector<double>(state.nblock));
  _y.resize(m, std::vector<double>(state.nblock));
  return *this;
}

Lbfgs& Lbfgs::setMaxIter(int maxIter) {
  Minimiser::setMaxIter(maxIter);
  return *this;
}


void Lbfgs::iteration() {
  if (iter == 0) _g0 = state.gradient();
  _i_cycle = iter % _m;

  // Find and take step
  getDirection();
  backtrackingLinesearch(state, _step, state.comm.dotProduct(_g0, _step));

  // Get new gradient
  _g1 = state.gradient();

  // Store the changes required for LBFGS
  _s[_i_cycle] = _step;
  _y[_i_cycle] = vec::diff(_g1, _g0);
  _rho[_i_cycle] = 1 / state.comm.dotProduct(_step, _y[_i_cycle]);

  _g0 = _g1;
}


void Lbfgs::getDirection() {
  double alpha[_m] = {0};
  int m_tmp = (_m < iter) ? _m : iter;

  if (iter == 0) {
    _step = vec::multiply(-_init_hessian, _g0);
    return;
  } else {
    _step = vec::multiply(-1, _g0);
  }

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (_i_cycle - 1 - i1 + _m) % _m;
    alpha[i] = _rho[i] * state.comm.dotProduct(_step, _s[i]);
    for (int j=0; j<state.ndof; j++) {
      _step[j] -= alpha[i] * _y[i][j];
    }
  }

  int i = (_i_cycle - 1 + _m) % _m;
  double gamma = 1 / (_rho[i] * state.comm.dotProduct(_y[i], _y[i]));
  _step = vec::multiply(gamma, _step);

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (_i_cycle - m_tmp + i1 + _m) % _m;
    double beta = _rho[i] * state.comm.dotProduct(_step, _y[i]);
    _step = vec::sum(_step, vec::multiply(alpha[i]-beta, _s[i]));
  }
}


bool Lbfgs::checkConvergence() {
  double sum = vec::dotProduct(_g0, _g0);
  double rms = sqrt(sum/state.ndof);
  return (rms < state.convergence);
}
