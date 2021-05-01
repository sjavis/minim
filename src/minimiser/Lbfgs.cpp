#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "linesearch.h"
#include "vec.h"


Lbfgs::Lbfgs(State &state, int m, int maxIter)
  : Minimiser(state, maxIter), _m(m), _rho(m),
    _s(m, std::vector<double>(state.ndof)),
    _y(m, std::vector<double>(state.ndof)),
    _g0(state.ndof), _g1(state.ndof), _step(state.ndof)
{}


void Lbfgs::iteration() {
  if (iter == 0) _g0 = state.gradient();
  _i_cycle = iter % _m;

  // Find and take step
  getDirection();
  backtrackingLinesearch(state, _step, vec::dotProduct(_g0, _step));

  // Get new gradient
  _g1 = state.gradient();

  // Store the changes required for LBFGS
  _s[_i_cycle] = _step;
  _y[_i_cycle] = vec::diff(_g1, _g0);
  _rho[_i_cycle] = 1 / vec::dotProduct(_step, _y[_i_cycle]);

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
    alpha[i] = _rho[i] * vec::dotProduct(_step, _s[i]);
    for (int j=0; j<state.ndof; j++) {
      _step[j] -= alpha[i] * _y[i][j];
    }
  }

  int i = (_i_cycle - 1 + _m) % _m;
  double gamma = 1 / (_rho[i] * vec::dotProduct(_y[i], _y[i]));
  _step = vec::multiply(gamma, _step);

  for (int i1=0; i1<m_tmp; i1++) {
    int i = (_i_cycle - m_tmp + i1 + _m) % _m;
    double beta = _rho[i] * vec::dotProduct(_step, _y[i]);
    _step = vec::sum(_step, vec::multiply(alpha[i]-beta, _s[i]));
  }
}
