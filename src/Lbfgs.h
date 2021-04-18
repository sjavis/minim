#ifndef LBFGS_H
#define LBFGS_H


// L-BFGS minimisation
class Lbfgs : public Minimiser {
  public:
    Lbfgs(System &sys, int m = 5, int maxIter = 10000);
    ~Lbfgs() {};

    void iteration();
    bool checkConvergence();

  private:
    int m;
    int i_cycle;
    std::vector<double> g;
    std::vector<double> rho;
    std::vector<std::vector<double>> s;
    std::vector<std::vector<double>> y;

    void getDirection(std::vector<double> &step);
    void linesearch(std::vector<double> &step);
};

#endif
