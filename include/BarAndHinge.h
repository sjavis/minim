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

      BarAndHinge& setFixed(const std::vector<bool>& fixed);

      Vector force = Vector(3);

    private:
      std::vector<bool> fixed;

      void stretching(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void bending(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void forceEnergy(const Vector& coords, const Element& el, double* e, Vector* g) const;
  };

}

#endif
