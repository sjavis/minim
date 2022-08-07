#include "Fire.h"

#include <math.h>
#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "vec.h"

namespace minim {

  Fire& Fire::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }

  Fire& Fire::setDtMax(double dt_max) {
    _dt_max = dt_max;
    return *this;
  }


  void Fire::init(State &state) {
    _v = std::vector<double>(state.ndof);
    _g = state.gradient();
    if (_dt_max == 0) _dt_max = 0.1 / sqrt(sqrt(state.comm.dotProduct(_g, _g)));
  }


  void Fire::iteration(State &state) {
    if (iter == 0) {
      _dt = _dt_max;
      _g = state.gradient();
      _gnorm = sqrt(state.comm.dotProduct(_g, _g));
    }
    double p = - state.comm.dotProduct(_v, _g);

    // Update velocity
    if (p > 0) {
      double vnorm = sqrt(state.comm.dotProduct(_v, _v));
      _v = (1-_a)*_v - (_a*vnorm/_gnorm + _dt)*_g;
      _n_steps++;
    } else {
      _v = -_dt * _g;
      _n_steps = 0;
    }

    // Adjust time step
    if ((p > 0) && (_n_steps > _n_min)) {
      _a = _a * _f_a;
      _dt = _dt * _f_inc;
      if (_dt > _dt_max) _dt = _dt_max;
    } else if (p <= 0) {
      _a = _a_start;
      _dt = _dt * _f_dec;
    }

    // Update position
    auto step = _dt * _v;
    state.blockCoords(state.blockCoords() + step);
    _g = state.gradient();
    _gnorm = sqrt(state.comm.dotProduct(_g, _g));
  }


  bool Fire::checkConvergence(const State &state) {
    double rms = _gnorm / sqrt(state.ndof);
    return (rms < state.convergence);
  }

}
