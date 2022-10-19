#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <memory>

namespace minim {
  class State;

  class Potential {
    typedef std::vector<double> Vector;
    typedef double (*EFunc)(const Vector&);
    typedef Vector (*GFunc)(const Vector&);
    typedef void (*EGFunc)(const Vector&, double*, Vector*);
    EFunc _energy;
    GFunc _gradient;
    EGFunc _energyGradient;

    public:
      Potential(EFunc energy, GFunc gradient) : _energy(energy), _gradient(gradient), energyDef(true), gradientDef(true) {};
      Potential(EGFunc energyGradient) : _energyGradient(energyGradient), energyGradientDef(true) {};
      ~Potential() {};
      virtual std::unique_ptr<Potential> clone() const {
        return std::make_unique<Potential>(*this);
      }

      virtual double energy(const Vector& coords) const;
      virtual Vector gradient(const Vector& coords) const;
      virtual void energyGradient(const Vector& coords, double* e, Vector* g) const;
      virtual double blockEnergy(const Vector& coords) const;
      virtual Vector blockGradient(const Vector& coords) const;
      virtual void blockEnergyGradient(const Vector& coords, double* e, Vector* g) const;

      bool energyDef = false;
      bool gradientDef = false;
      bool energyGradientDef = false;
      bool blockEnergyDef = false;
      bool blockGradientDef = false;
      bool blockEnergyGradientDef = false;


      State newState(int ndof);
      virtual State newState(const Vector& coords);


      struct Element {
        int id;
        int type;
        std::vector<int> idof;
        std::vector<double> parameters;
      };
      std::vector<Element> elements;
      std::vector<Element> elements_halo;

      Potential& setElements(std::vector<Element> elements);
      Potential& setElements(std::vector<std::vector<int>> idofs);
      Potential& setElements(std::vector<std::vector<int>> idofs, std::vector<int> types,
                             std::vector<std::vector<double>> parameters);

    protected:
      Potential() {};
  };


  // An intermediate class is used to return the derived type for methods that return a Potential
  template<typename Derived>
  class NewPotential : public Potential {
    public:
      std::unique_ptr<Potential> clone() const override {
        return std::make_unique<Derived>(static_cast<const Derived&>(*this));
      }

      Derived& setElements(std::vector<Element> elements) {
        return static_cast<Derived&>(Potential::setElements(elements));
      }
      Derived& setElements(std::vector<std::vector<int>> idofs) {
        return static_cast<Derived&>(Potential::setElements(idofs));
      }
      Derived& setElements(std::vector<std::vector<int>> idofs, std::vector<int> types,
                           std::vector<std::vector<double>> parameters) {
        return static_cast<Derived&>(Potential::setElements(idofs, types, parameters));
      }
  };

}

#endif
