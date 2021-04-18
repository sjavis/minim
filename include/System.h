#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>

class System {
  public:
    const int ndof;
    double convergence_rms = 1e-6;
    std::vector<double> state;

    System(int ndof_);
    System(int ndof_, std::vector<double> state_);
    virtual ~System() {};

    virtual double energy(std::vector<double> state_) = 0;
    virtual void gradient(std::vector<double> state_, std::vector<double> &g) = 0;
    double energy() { return energy(state); };
    void gradient(std::vector<double> &g) { gradient(state, g); }

};

#endif
