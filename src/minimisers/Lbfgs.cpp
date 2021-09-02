#include <vector>
#include <math.h>
#include "State.h"
#include "Minimiser.h"
#include "Lbfgs.h"
#include "linesearch.h"
#include "Communicator.h"
#include "vec.h"
#include "minimMpi.h"


Lbfgs::Lbfgs(State &state, AdjustFunc adjustModel)
  : Minimiser(state, adjustModel), _m(5)
{
  if (minim::mpi.rank == 0) {
    _s = std::vector<std::vector<double>>(_m, std::vector<double>(state.ndof));
    _y = std::vector<std::vector<double>>(_m, std::vector<double>(state.ndof));
    _rho = std::vector<double>(_m);
    _g0 = std::vector<double>(state.ndof);
    _g1 = std::vector<double>(state.ndof);
  }
}


Lbfgs& Lbfgs::setM(int m) {
  _m = m;
  if (minim::mpi.rank == 0) {
    _rho.resize(m);
    _s.resize(m, std::vector<double>(state.ndof));
    _y.resize(m, std::vector<double>(state.ndof));
  }
  return *this;
}

Lbfgs& Lbfgs::setMaxIter(int maxIter) {
  Minimiser::setMaxIter(maxIter);
  return *this;
}


void Lbfgs::iteration() {
  if (iter == 0) _g0 = state.comm.gather(state.gradient(), 0);
  _i_cycle = iter % _m;

  // Find minimisation direction
  std::vector<double> step = getDirection();
  std::vector<double> step_block = state.comm.scatter(step, 0);

  // Perform linesearch
  double slope;
  if (minim::mpi.rank==0) slope = vec::dotProduct(_g0, step) / vec::norm(step);
  state.comm.bcast(slope);
  double step_multiplier = backtrackingLinesearch(state, step_block, slope);

  // Get new gradient
  _g1 = state.comm.gather(state.gradient(), 0);

  // Store the changes required for LBFGS
  if (minim::mpi.rank == 0) {
    _s[_i_cycle] = vec::multiply(step_multiplier, step);
    _y[_i_cycle] = vec::diff(_g1, _g0);
    _rho[_i_cycle] = 1 / vec::dotProduct(_s[_i_cycle], _y[_i_cycle]);
  }

  _g0 = _g1;
}


std::vector<double> Lbfgs::getDirection() {
  std::vector<double> step;

  // Compute the step on the main processor
  if (minim::mpi.rank == 0) {
    double alpha[_m] = {0};
    int m_tmp = (_m < iter) ? _m : iter;

    if (iter == 0) {
      step = vec::multiply(-_init_hessian, _g0);
      return step;
    } else {
      step = vec::multiply(-1, _g0);
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

  return step;
}


bool Lbfgs::checkConvergence() {
  double rms;
  if (minim::mpi.rank == 0) rms = sqrt(vec::dotProduct(_g0, _g0) / state.ndof);
  state.comm.bcast(rms);
  return (rms < state.convergence);
}
