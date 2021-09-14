#ifndef STATE_H
#define STATE_H

#include <vector>
#include "Communicator.h"

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

    double energy();
    double energy(const Vector &coords);

    Vector gradient();
    Vector gradient(const Vector &coords);

    double operator[](int i);

    Vector getCoords();
    void setCoords(const Vector &in);

    Vector blockCoords();
    void blockCoords(const Vector &in);

    void communicate();

  private:
    int _istart;
    int _iend;
    Vector _coords;
    Potential &_pot;
};

#endif
