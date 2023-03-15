#include "Lbfgs.h"

#include <math.h>
#include "State.h"
#include "linesearch.h"
#include "utils/vec.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  Lbfgs& Lbfgs::setM(int m) {
    _m = m;
    return *this;
  }

  Lbfgs& Lbfgs::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }

  Lbfgs& Lbfgs::setMaxStep(double maxStep) {
    _maxStep = maxStep;
    return *this;
  }


  void Lbfgs::init(State& state) {
    if (state.comm.rank() == 0) {
      _root = true;
      _s = vector2d<double>(_m, vector<double>(state.ndof));
      _y = vector2d<double>(_m, vector<double>(state.ndof));
      _rho = vector<double>(_m);
    } else {
      _root = false;
    }
  }


  void Lbfgs::iteration(State& state) {
    if (iter == 0) {
      _g = state.gradient();
      _i = 0;
    } else {
      _i++;
    }

    // Find minimisation direction
    vector<double> step = getDirection();
    // Ensure it is going downhill
    double gs;
    if (_root) {
      gs = vec::dotProduct(_g, step);
      if (gs > 0) {
        gs = -gs;
        step = -step;
      }
    }

    // Perform linesearch
    vector<double> step_block = state.comm.scatter(step, 0);
    if (linesearch == "backtracking") {
      state.comm.bcast(gs);
      double step_multiplier = backtrackingLinesearch(state, step_block, gs);
      if (_root) step *= step_multiplier;
    } else {
      state.blockCoords(state.blockCoords() + step_block);
    }

    // Get new gradient
    _gNew = state.gradient();

    // Store the changes required for LBFGS
    if (_root) {
      auto s = step;
      auto y = _gNew - _g;
      double sy = vec::dotProduct(s, y);
      if (sy != 0) {
        int i_cycle = _i % _m;
        _s[i_cycle] = s;
        _y[i_cycle] = y;
        _rho[i_cycle] = 1 / sy;
      } else {
        _i --;
      }
    }

    _g = _gNew;
  }


  vector<double> Lbfgs::getDirection() {
    vector<double> step;

    // Compute the step on the main processor
    if (_root) {
      double alpha[_m] = {0};
      int m_tmp = std::min(_m, _i);
      int i_cycle = _i % _m;

      if (m_tmp == 0) {
        step = -_init_hessian * _g;
        return step;
      }

      step = -_g;
      int ndof = step.size();
      for (int i1=0; i1<m_tmp; i1++) {
        int i = (i_cycle - 1 - i1 + _m) % _m;
        alpha[i] = _rho[i] * vec::dotProduct(step, _s[i]);
        for (int j=0; j<ndof; j++) {
          step[j] -= alpha[i] * _y[i][j];
        }
      }

      int i = (i_cycle - 1 + _m) % _m;
      double gamma = 1 / (_rho[i] * vec::dotProduct(_y[i], _y[i]));
      step = gamma * step;

      for (int i1=0; i1<m_tmp; i1++) {
        int i = (i_cycle - m_tmp + i1 + _m) % _m;
        double beta = _rho[i] * vec::dotProduct(step, _y[i]);
        step += (alpha[i]-beta) * _s[i];
      }

      if (_maxStep != 0) {
        double stepSize = vec::norm(step);
        if (stepSize > _maxStep) step *= _maxStep / stepSize;
      }
    }

    return step;
  }


  bool Lbfgs::checkConvergence(const State& state) {
    double rms = sqrt(vec::dotProduct(_g, _g) / state.ndof);
    return (rms < state.convergence);
  }

}
