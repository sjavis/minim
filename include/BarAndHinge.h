#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include "Potential.h"


class BarAndHinge : public Potential {
  public:
    BarAndHinge() {};
    ~BarAndHinge() {};

    double energy(const std::vector<double> &coords, const Args &args) override;

    std::vector<double> gradient(const std::vector<double> &coords, const Args &args) override;

  private:
    void stretching(const std::vector<double> &coords, Args::Element el, double *e, std::vector<double> *g);
    void bending(const std::vector<double> &coords, Args::Element el, double *e, std::vector<double> *g);
};

#endif
