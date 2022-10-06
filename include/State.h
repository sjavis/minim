#ifndef STATE_H
#define STATE_H

#include <vector>
#include "Communicator.h"

namespace minim {
  class Potential;

  class State {
    typedef std::vector<double> Vector;

    public:
      int ndof;
      double convergence = 1e-6;
      std::unique_ptr<Potential> pot;
      Communicator comm;

      State(const Potential& pot, const Vector& coords);
      State(const State& state);
      State& operator=(const State& state);
      ~State() {};

      double energy() const;
      double energy(const Vector& coords) const;

      Vector gradient() const;
      Vector gradient(const Vector& coords) const;

      double operator[](int i);

      Vector getCoords() const;
      void setCoords(const Vector& in);

      Vector blockCoords() const;
      void blockCoords(const Vector& in);

      void communicate();

    private:
      int _istart;
      int _iend;
      Vector _coords;
  };

}

#endif
