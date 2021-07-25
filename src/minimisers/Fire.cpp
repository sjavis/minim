#include <math.h>
#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "GradDescent.h"
#include "vec.h"


Fire::Fire(State &state)
: Minimiser(state), _v(state.ndof)
{}


Fire& Fire::setDtMax(double dt_max) {
  _dt_max = dt_max;
  return *this;
}


void Fire::iteration() {
  if (iter == 0) {
    _g = state.gradient();
    _gnorm = vec::norm(_g);
  }
  double p = - vec::dot_product(_v, _g);

  // Update velocity
  if (p > 0) {
    _v = vec::diff( (1-_a)*_v, vec::multiply(_a*vec::norm(_v)/_gnorm + dt, _g) );
    _n_steps++;
  } else {
    _v = vec::multiply(-1*_dt, _g);
    _n_steps = 0;
  }

  // Adjust time step
  if ((p > 0) && (_n_steps > _n_min)) {
    _dt = min(_dt*_f_inc, _dt_max);
    _a = _a * _f_a;
  } else if (p <= 0) {
    _dt = _dt * _f_dec;
    _a = _a_start;
  }

  // Update position
  auto step = vec::multiply(_dt, _v);
  state.blockCoords(vec::sum(state.blockCoords(), step));
  _g = state.gradient();
  _gnorm = vec::norm(_g);
}


bool Fire::checkConvergence() {
  double rms = _gnorm / sqrt(state.ndof);
  return (rms < state.convergence);
}
