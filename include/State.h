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
    int nblock;
    double convergence = 1e-6;
    Args &args;
    Communicator comm;

    State(Potential &pot, Vector coords, Args &args);
    ~State() {};

    double energy();
    double energy(Vector coords);

    Vector gradient();
    Vector gradient(Vector coords);

    double operator[](int i);

    Vector getCoords();
    void setCoords(Vector in);

    Vector blockCoords();
    void blockCoords(Vector in);

  private:
    int _istart;
    int _iend;
    Vector _coords;
    Potential &_pot;
};

#endif
