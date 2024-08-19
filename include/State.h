#ifndef STATE_H
#define STATE_H

#include <vector>
#include <memory>
#include <cstddef>
#include "Potential.h"
#include "Communicator.h"

namespace minim {
  class Potential;
  using std::vector;

  class State {
    public:
      size_t ndof;
      double convergence;
      std::unique_ptr<Potential> pot;
      std::unique_ptr<Communicator> comm;
      bool usesThisProc = true;


      State(const Potential& pot, const vector<double>& coords, const vector<int>& ranks={});

      // Copy constructor / assignment
      State(const State& state);
      State& operator=(const State& state);

      // Energy / Gradient
      double energy() const;
      double energy(const vector<double>& coords) const;
      vector<double> gradient() const;
      vector<double> gradient(const vector<double>& coords) const;
      void energyGradient(double* e, vector<double>* g) const;
      void energyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      // Coordinates
      double operator[](int i);
      vector<double> coords() const;
      void coords(const vector<double>& in);

      // Parallel functions
      vector<double> blockCoords() const;
      void blockCoords(const vector<double>& in);

      double blockEnergy() const;
      double blockEnergy(const vector<double>& coords) const;
      vector<double> blockGradient() const;
      vector<double> blockGradient(const vector<double>& coords) const;
      void blockEnergyGradient(double* e, vector<double>* g) const;
      void blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      double procEnergy() const;
      double procEnergy(const vector<double>& coords) const;
      vector<double> procGradient() const;
      vector<double> procGradient(const vector<double>& coords) const;
      void procEnergyGradient(double* e, vector<double>* g) const;
      void procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      double allEnergy() const;
      vector<double> allGradient() const;
      void allEnergyGradient(double* e, vector<double>* g) const;
      vector<double> allCoords() const;

      double componentEnergy(int component) const;

      void communicate();

      bool isFailed = false;
      void failed();

      vector<double> _coords;
  };

}

#endif
