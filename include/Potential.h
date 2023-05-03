#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <memory>
#include <functional>

namespace minim {
  class State;
  class Communicator;

  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  class Potential {
    typedef double (*EFunc)(const vector<double>&);
    typedef vector<double> (*GFunc)(const vector<double>&);
    typedef void (*EGFunc)(const vector<double>&, double*, vector<double>*);
    EFunc _energy;
    GFunc _gradient;
    EGFunc _energyGradient;

    public:
      // Basic user functions
      Potential(EFunc energy, GFunc gradient);
      Potential(EGFunc energyGradient);

      State newState(int ndof, const vector<int>& ranks={});
      virtual State newState(const vector<double>& coords, const vector<int>& ranks={});

      virtual double energy(const vector<double>& coords) const;
      virtual vector<double> gradient(const vector<double>& coords) const;
      virtual void energyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      // Initialisation
      bool serialDef() const;
      bool parallelDef() const;

      virtual void init(const vector<double>& coords) {};
      virtual void distributeParameters(const Communicator& comm) {}; // Take care using this, if the potential is cloned any distributed parameters will be copied as they are
      bool distributed = false;

      // Energy elements for parallelisation
      struct Element {
        int type;
        vector<int> idof;
        vector<double> parameters;
      };
      vector<Element> elements;
      vector<Element> elements_halo;

      Potential& setElements(vector<Element> elements);
      Potential& setElements(vector2d<int> idofs);
      Potential& setElements(vector2d<int> idofs, vector<int> types, vector2d<double> parameters);

      virtual void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      virtual void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {};

      // Constraints
      struct Constraint {
        vector<int> idof;
        std::function<vector<double>(const vector<double>&)> normal;
        std::function<void(vector<double>&)> correction = nullptr;
      };
      vector<Constraint> constraints;

      Potential& setConstraints(vector<Constraint> constraints);
      Potential& setConstraints(vector<int> iFix);
      Potential& setConstraints(vector2d<int> idofs, vector<double> normal, std::function<void(vector<double>&)> correction=nullptr);
      Potential& setConstraints(vector2d<int> idofs, std::function<vector<double>(const vector<double>&)> normal, std::function<void(vector<double>&)> correction=nullptr);

      // Copy / destruct
      ~Potential() {};
      virtual std::unique_ptr<Potential> clone() const {
        return std::make_unique<Potential>(*this);
      }

    protected:
      Potential() : _energy(nullptr), _gradient(nullptr), _energyGradient(nullptr) {};

      bool _energyDef = false;
      bool _energyGradientDef = false;
      bool _parallelDef = false;
  };


  // An intermediate class is used to return the derived type for methods that return a Potential
  template<typename Derived>
  class NewPotential : public Potential {
    public:
      std::unique_ptr<Potential> clone() const override {
        return std::make_unique<Derived>(static_cast<const Derived&>(*this));
      }

      Derived& setElements(vector<Element> elements) {
        return static_cast<Derived&>(Potential::setElements(elements));
      }
      Derived& setElements(vector2d<int> idofs) {
        return static_cast<Derived&>(Potential::setElements(idofs));
      }
      Derived& setElements(vector2d<int> idofs, vector<int> types, vector2d<double> parameters) {
        return static_cast<Derived&>(Potential::setElements(idofs, types, parameters));
      }

      Derived& setConstraints(vector<Constraint> constraints) {
        return static_cast<Derived&>(Potential::setConstraints(constraints));
      }
      Derived& setConstraints(vector<int> iFix) {
        return static_cast<Derived&>(Potential::setConstraints(iFix));
      }
      Derived& setConstraints(vector2d<int> idofs, vector<double> normal, std::function<void(vector<double>&)> correction=nullptr) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal, correction));
      }
      Derived& setConstraints(vector2d<int> idofs, std::function<vector<double>(const vector<double>&)> normal, std::function<void(vector<double>&)> correction=nullptr) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal, correction));
      }
  };

}

#endif
