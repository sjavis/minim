#ifndef STATE_H
#define STATE_H

#include <vector>
#include "Communicator.h"

namespace minim {
  class Potential;
  class Args;

  class State {
    typedef std::vector<double> Vector;

    public:
      int ndof;
      double convergence = 1e-6;
      Args &args;
      Communicator comm;

      State(Potential &pot, const Vector &coords, Args &args);
      ~State() {};

      double energy() const;
      double energy(const Vector &coords) const;

      Vector gradient() const;
      Vector gradient(const Vector &coords) const;

      double operator[](int i);

      Vector getCoords() const;
      void setCoords(const Vector &in);

      Vector blockCoords() const;
      void blockCoords(const Vector &in);

      void communicate();

    private:
      int _istart;
      int _iend;
      Vector _coords;
      Potential &_pot;
  };

}

#endif
