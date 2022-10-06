#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {

  class BarAndHinge : public NewPotential<BarAndHinge> {
    typedef std::vector<double> Vector;

    public:
      BarAndHinge() {};
      ~BarAndHinge() {};

      double energy(const Vector& coords) const override;
      Vector gradient(const Vector& coords) const override;

    private:
      void stretching(const Vector& coords, Potential::Element el, double* e, Vector* g) const;
      void bending(const Vector& coords, Potential::Element el, double* e, Vector* g) const;
  };

}

#endif
