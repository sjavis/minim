#ifndef LBFGS_H
#define LBFGS_H


// l-BFGS minimisation
class Lbfgs : public Minimiser {
  public:
    Lbfgs(System &sys, int maxIter = 10000);
    ~Lbfgs() {};

    void iteration();
    bool checkConvergence();
};

#endif
