#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {

  class BarAndHinge : public Potential {
    typedef std::vector<double> Vector;

    public:
      BarAndHinge() {};
      ~BarAndHinge() {};

      double energy(const Vector &coords, const Args &args) const override;

      Vector gradient(const Vector &coords, const Args &args) const override;

    private:
      void stretching(const Vector &coords, Args::Element el, double *e, Vector *g) const;
      void bending(const Vector &coords, Args::Element el, double *e, Vector *g) const;
  };

}

#endif
