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
      double convergence = 1e-6;

      // Basic user functions
      Potential(EFunc energy, GFunc gradient);
      Potential(EGFunc energyGradient);

      virtual double energy(const vector<double>& coords) const;
      virtual vector<double> gradient(const vector<double>& coords) const;
      virtual void energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const;
      // Note: In any derived classes, either energyGradient or (energy and gradient) MUST be overridden

      State newState(int ndof, const vector<int>& ranks={});
      State newState(const vector<double>& coords, const vector<int>& ranks={});

      // Constraints
      struct Constraint {
        using NormalFn = std::function<vector<double>(const vector<int>&, const vector<double>&)>;
        using CorrectionFn = std::function<void(const vector<int>&, vector<double>&)>;
        vector<int> idof;
        NormalFn normalFn;
        CorrectionFn correction = nullptr;
        vector<double> normal(const vector<double>& normalVec) const { return normalFn(idof, normalVec); }
      };
      vector<Constraint> constraints;

      Potential& setConstraints(vector<int> iFix);
      Potential& setConstraints(vector2d<int> idofs, vector<double> normal);
      Potential& setConstraints(vector2d<int> idofs, Constraint::NormalFn normal, Constraint::CorrectionFn correction=nullptr);

      const vector<double>& applyConstraints(const vector<double>& coords, vector<double>& grad) const;
      const vector<double>& correctConstraints(vector<double>& coords) const;

      bool isFixed(int index) const;
      vector<char> isFixed(const vector<int>& indicies) const;

      // Internal
      virtual void init(const vector<double>& coords) {};
      virtual void initLocal(const vector<double>& coords, const Communicator& comm) {}; // Take care using this, if the potential is cloned any distributed parameters will be copied as they are
      bool isSerial() const;

      // Copy / destruct
      virtual ~Potential() = default;
      virtual std::unique_ptr<Potential> clone() const {
        return std::make_unique<Potential>(*this);
      }


      // Members and functions specific to different types of potential
      enum{
        SERIAL = 0,
        UNSTRUCTURED = 1,
        GRID = 2,
      };
      virtual int potentialType() const { return SERIAL; };
      std::unique_ptr<Communicator> newComm() const;

      // UNSTRUCTURED: Energy elements for parallelisation
      bool distributed = false;

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

      virtual void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {};
      virtual void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {};

      // GRID
      vector<int> gridSize;
      int dofPerNode = 1;

    protected:
      Potential() : _energy(nullptr), _gradient(nullptr), _energyGradient(nullptr) {};
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

      Derived& setConstraints(vector<int> iFix) {
        return static_cast<Derived&>(Potential::setConstraints(iFix));
      }
      Derived& setConstraints(vector2d<int> idofs, vector<double> normal) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal));
      }
      Derived& setConstraints(vector2d<int> idofs, Constraint::NormalFn normal, Constraint::CorrectionFn correction=nullptr) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal, correction));
      }
  };

}

#endif
