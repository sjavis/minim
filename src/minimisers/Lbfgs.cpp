#include "minimisers/Lbfgs.h"

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

  Lbfgs& Lbfgs::setInitStep(double initStep) {
    _initStep = initStep;
    return *this;
  }


  void Lbfgs::init(State& state) {
    _s = vector2d<double>(_m);
    _y = vector2d<double>(_m);
    _rho = vector<double>(_m);
  }


  void Lbfgs::iteration(State& state) {
    if (iter == 0) {
      _g = state.procGradient();
      _i = 0;
    } else {
      _i++;
    }

    // Find minimisation direction
    vector<double> step = getDirection(*state.comm);
    state.applyConstraints(step);
    // Ensure it is going downhill
    double gs = state.comm->dotProduct(_g, step);
    if (gs > 0) {
      gs = -gs;
      step = -step;
    }

    // Perform linesearch
    if (linesearch == "backtracking") {
      backtrackingLinesearch(state, step, gs);
    } else {
      state.blockCoords(state.blockCoords() + step);
    }

    // Get new gradient
    vector<double> gNew = state.procGradient();

    // Store the changes required for LBFGS
    double sy = state.comm->dotProduct(step, gNew-_g);
    if (sy != 0) {
      int i_cycle = _i % _m;
      _s[i_cycle] = step;
      _y[i_cycle] = gNew - _g;
      _rho[i_cycle] = 1 / sy;
    } else {
      _i --;
    }

    _g = gNew;
  }


  vector<double> Lbfgs::getDirection(const Communicator& comm) {
    vector<double> alpha(_m);
    int m_tmp = std::min(_m, _i);
    int i_cycle = _i % _m;

    vector<double> step = -_g;

    if (_i == 0) {
      // First iteration: Directly set the step size
      double gnorm = comm.norm(_g);
      if (gnorm > 0) step *= _initStep / gnorm;
      return step;
    }

    int ndof = step.size();
    for (int i1=0; i1<m_tmp; i1++) {
      int i = i_cycle - 1 - i1; // Go backwards from newest
      i = (i + _m) % _m; // Wrap around
      alpha[i] = _rho[i] * comm.dotProduct(step, _s[i]);
      #pragma omp parallel for simd
      for (int j=0; j<ndof; j++) {
        step[j] -= alpha[i] * _y[i][j];
      }
    }

    int i = (i_cycle - 1 + _m) % _m;
    double gamma = 1 / (_rho[i] * comm.dotProduct(_y[i], _y[i]));
    step *= gamma;

    for (int i1=0; i1<m_tmp; i1++) {
      int i = i_cycle - m_tmp + i1; // Go forwards from oldest
      i = (i + _m) % _m; // Wrap around
      double beta = _rho[i] * comm.dotProduct(step, _y[i]);
      step += (alpha[i]-beta) * _s[i];
    }

    // Cap the max step size (if using)
    if (_maxStep != 0) {
      double stepSize = comm.norm(step);
      if (stepSize > _maxStep) step *= _maxStep / stepSize;
    }

    return step;
  }


  bool Lbfgs::checkConvergence(const State& state) {
    double rms = sqrt(state.comm->dotProduct(_g, _g) / state.ndof);
    return (rms < state.convergence);
  }

}
