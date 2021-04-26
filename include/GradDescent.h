#ifndef GRADDESCENT_H
#define GRADDESCENT_H

#include <vector>


// Gradient desecent minimisation
class GradDescent : public Minimiser {
  public:
    double alpha;

    GradDescent(State &state, double alpha = 1e-3, int maxIter = 10000);
    ~GradDescent() {};

    void iteration();

  private:
    std::vector<double> _g;
};

#endif
