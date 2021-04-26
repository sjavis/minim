#ifndef STATE_H
#define STATE_H

#include <vector>

class Potential;
class Args;

class State {
  typedef std::vector<double> Vector;

  public:
    Vector coords;
    const int ndof;
    double convergence = 1e-6;

    State(Potential &pot, Vector coords, Args &args);
    ~State() {};

    double energy();
    double energy(Vector coords);

    Vector gradient();
    Vector gradient(Vector coords);

    double operator[](int i);

  private:
    Potential &_pot;
    Args &_args;
};

#endif
