#ifndef GRADDESCENT_H
#define GRADDESCENT_H


// Gradient desecent minimisation
class GradDescent : public Minimiser {
  public:
    double alpha;

    GradDescent(System &sys, double a = 1e-3, int maxIter = 10000);
    ~GradDescent() {};

    void iteration();
    bool checkConvergence();
};

#endif
