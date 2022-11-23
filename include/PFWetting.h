#ifndef PFWETTING_H
#define PFWETTING_H

#include <vector>
#include "Potential.h"

namespace minim {

  class PFWetting : public NewPotential<PFWetting> {
    typedef std::vector<double> Vector;

    public:
      int[3] gridSize;
      double epsilon;

      PFWetting() { _blockEnergyGradientDef = true; };
      ~PFWetting() {};

      void blockEnergyGradient(const Vector& coords, double* e, Vector* g) const override;

      State newState(const Vector& coords) override;
    
    private:
      void elementEnergyGradient(const Element el, const Vector& coords, double* e, Vector* g) const;
      Vector getCoord(int i) const;
  };

}

#endif
