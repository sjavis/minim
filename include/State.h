#ifndef STATE_H
#define STATE_H

#include <vector>
#include "Communicator.h"

namespace minim {
  class Potential;

  class State {
    typedef std::vector<double> Vector;

    public:
      size_t ndof;
      double convergence = 1e-6;
      std::unique_ptr<Potential> pot;
      Communicator comm;

      State(const Potential& pot, const Vector& coords);
      State(const State& state);
      State& operator=(const State& state);
      ~State() {};

      double energy() const;
      double energy(const Vector& coords) const;
      Vector gradient() const;
      Vector gradient(const Vector& coords) const;
      void energyGradient(double* e, Vector* g) const;
      void energyGradient(const Vector& coords, double* e, Vector* g) const;

      double blockEnergy() const;
      double blockEnergy(const Vector& coords) const;
      Vector blockGradient() const;
      Vector blockGradient(const Vector& coords) const;
      void blockEnergyGradient(double* e, Vector* g) const;
      void blockEnergyGradient(const Vector& coords, double* e, Vector* g) const;

      double procEnergy() const;
      double procEnergy(const Vector& coords) const;
      Vector procGradient() const;
      Vector procGradient(const Vector& coords) const;
      void procEnergyGradient(double* e, Vector* g) const;
      void procEnergyGradient(const Vector& coords, double* e, Vector* g) const;

      double operator[](int i);

      Vector getCoords() const;
      void setCoords(const Vector& in);

      Vector blockCoords() const;
      void blockCoords(const Vector& in);

      void communicate();

    private:
      Vector _coords;
  };

}

#endif
