#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {

  class BarAndHinge : public NewPotential<BarAndHinge> {
    typedef std::vector<double> Vector;

    public:
      BarAndHinge() { _parallelDef = true; };
      ~BarAndHinge() {};

      void elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const override;

    private:
      void stretching(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void bending(const Vector& coords, const Element& el, double* e, Vector* g) const;
  };

}

#endif
