#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {

  class BarAndHinge : public NewPotential<BarAndHinge> {
    typedef std::vector<double> Vector;

    public:
      BarAndHinge() { _blockEnergyGradientDef = true; };
      ~BarAndHinge() {};

      void blockEnergyGradient(const Vector& coords, double* e, Vector* g) const override;

    private:
      void stretching(Potential::Element el, const Vector& coords, double* e, Vector* g) const;
      void bending(Potential::Element el, const Vector& coords, double* e, Vector* g) const;
  };

}

#endif
