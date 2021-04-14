#ifndef SYSTEM_H
#define SYSTEM_H

class System {
  public:
    const int ndof;
    double* state;

    System(int ndof_);
    System(int ndof_, double* state_);
    virtual ~System();

    virtual double energy(double* state_) = 0;
    virtual void gradient(double* state_, double* g) = 0;
    double energy() { return energy(state); };
    void gradient(double* g) { return gradient(state, g); }

};

#endif
