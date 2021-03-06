#include "Lbfgs.h"

#include <vector>
#include <math.h>
#include "State.h"
#include "Minimiser.h"
#include "Communicator.h"
#include "linesearch.h"
#include "vec.h"
#include "utils/mpi.h"

namespace minim {

  typedef std::vector<double> Vector;


  Lbfgs::Lbfgs(State &state, AdjustFunc adjustModel)
    : Minimiser(state, adjustModel), _m(5)
  {
    if (minim::mpi.rank == 0) {
      _s = std::vector<Vector>(_m, Vector(state.ndof));
      _y = std::vector<Vector>(_m, Vector(state.ndof));
      _rho = Vector(_m);
      _g0 = Vector(state.ndof);
      _g1 = Vector(state.ndof);
    }
  }


  Lbfgs& Lbfgs::setM(int m) {
    _m = m;
    if (minim::mpi.rank == 0) {
      _rho.resize(m);
      _s.resize(m, Vector(state.ndof));
      _y.resize(m, Vector(state.ndof));
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
    Vector step = getDirection();
    Vector step_block = state.comm.scatter(step, 0);

    // Perform linesearch
    double de0;
    if (minim::mpi.rank==0) de0 = vec::dotProduct(_g0, step);
    state.comm.bcast(de0);
    double step_multiplier = backtrackingLinesearch(state, step_block, de0);

    // Get new gradient
    _g1 = state.comm.gather(state.gradient(), 0);

    // Store the changes required for LBFGS
    if (minim::mpi.rank == 0) {
      _s[_i_cycle] = step_multiplier * step;
      _y[_i_cycle] = _g1 - _g0;
      _rho[_i_cycle] = 1 / vec::dotProduct(_s[_i_cycle], _y[_i_cycle]);
    }

    _g0 = _g1;
  }


  Vector Lbfgs::getDirection() {
    Vector step;

    // Compute the step on the main processor
    if (minim::mpi.rank == 0) {
      double alpha[_m] = {0};
      int m_tmp = (_m < iter) ? _m : iter;

      if (iter == 0) {
        step = -_init_hessian * _g0;
        return step;
      }

      step = -_g0;
      for (int i1=0; i1<m_tmp; i1++) {
        int i = (_i_cycle - 1 - i1 + _m) % _m;
        alpha[i] = _rho[i] * vec::dotProduct(step, _s[i]);
        for (int j=0; j<state.ndof; j++) {
          step[j] -= alpha[i] * _y[i][j];
        }
      }

      int i = (_i_cycle - 1 + _m) % _m;
      double gamma = 1 / (_rho[i] * vec::dotProduct(_y[i], _y[i]));
      step = gamma * step;

      for (int i1=0; i1<m_tmp; i1++) {
        int i = (_i_cycle - m_tmp + i1 + _m) % _m;
        double beta = _rho[i] * vec::dotProduct(step, _y[i]);
        step += (alpha[i]-beta) * _s[i];
      }

      if (vec::dotProduct(step, _g0) > 0) step = -step;
    }

    return step;
  }


  bool Lbfgs::checkConvergence() {
    double rms;
    if (minim::mpi.rank == 0) rms = sqrt(vec::dotProduct(_g0, _g0) / state.ndof);
    state.comm.bcast(rms);
    return (rms < state.convergence);
  }

}
