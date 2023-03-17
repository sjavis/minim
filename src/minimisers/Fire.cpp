#include "Fire.h"

#include <math.h>
#include "State.h"
#include "linesearch.h"
#include "utils/vec.h"
#include "utils/print.h"

namespace minim {

  Fire::Fire(double dtMax) : dtMax(dtMax) {
    linesearch = "none";
  }


  Fire& Fire::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }

  Fire& Fire::setDtMax(double dtMax) {
    this->dtMax = dtMax;
    return *this;
  }


  void Fire::init(State& state) {
    _v = std::vector<double>(state.ndof);
    _g = state.gradient();
    _gNorm = vec::norm(_g);
    if (dtMax != 0) return;
    if (_gNorm!=0) {
      dtMax = 0.1 / sqrt(_gNorm);
    } else {
      print("Warning: Unable to estimate dtMax parameter for FIRE.");
      dtMax = 1;
    }
  }


  void Fire::iteration(State& state) {
    if (iter == 0) {
      _dt = dtMax;
      _g = state.gradient();
      _gNorm = vec::norm(_g);
    }
    double p = - vec::dotProduct(_v, _g);

    // Update velocity
    if (p > 0) {
      double vNorm = vec::norm(_v);
      _v = (1-_a)*_v - (_a*vNorm/_gNorm + _dt)*_g;
      _nSteps++;
    } else {
      _v = -_dt * _g;
      _nSteps = 0;
    }

    // Adjust time step
    if ((p > 0) && (_nSteps > _nMin)) {
      _a = _a * _fA;
      _dt = _dt * _fInc;
      if (_dt > dtMax) _dt = dtMax;
    } else if (p <= 0) {
      _a = _aStart;
      _dt = _dt * _fDec;
    }

    // Get step
    auto step = _dt * _v;
    auto stepBlock = state.comm.scatter(step);

    // Perform linesearch (if set)
    if (linesearch == "backtracking") {
      double gs = vec::dotProduct(_g, step);
      backtrackingLinesearch(state, stepBlock, gs);
    } else {
      state.blockCoords(state.blockCoords() + stepBlock);
    }

    // Update gradient
    _g = state.gradient();
    _gNorm = vec::norm(_g);
  }


  bool Fire::checkConvergence(const State& state) {
    double rms = _gNorm / sqrt(state.ndof);
    return (rms < state.convergence);
  }

}
