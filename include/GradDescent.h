#ifndef GRADDESCENT_H
#define GRADDESCENT_H

#include <vector>


// Gradient desecent minimisation
class GradDescent : public Minimiser {
  public:
    GradDescent(State &state);
    ~GradDescent() {};

    GradDescent& setAlpha(double alpha);
    GradDescent& setMaxIter(int maxIter);

    void iteration();

  private:
    double _alpha = 1e-3;
    std::vector<double> _g;
};

#endif
