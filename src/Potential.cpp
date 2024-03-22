#include "Potential.h"

#include "State.h"
#include "utils/vec.h"
#include <stdexcept>

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  Potential::Potential(EFunc energy, GFunc gradient)
    : _energy(energy), _gradient(gradient), _energyGradient(nullptr)
  {
    _energyDef = true;
  }

  Potential::Potential(EGFunc energyGradient)
    : _energy(nullptr), _gradient(nullptr), _energyGradient(energyGradient)
  {
    _energyGradientDef = true;
  }


  double Potential::energy(const vector<double>& coords) const {
    if (_energyDef) {
      if (_energy == nullptr) throw std::logic_error("Energy function marked as defined but not called.");
      return _energy(coords);
    } else {
      if (!_energyGradientDef) throw std::logic_error("Energy function not defined.");
      double e;
      energyGradient(coords, &e, nullptr);
      return e;
    }
  }


  vector<double> Potential::gradient(const vector<double>& coords) const {
    if (_energyDef) {
      if (_gradient == nullptr) throw std::logic_error("Gradient function marked as defined but not called.");
      return _gradient(coords);
    } else {
      if (!_energyGradientDef) throw std::logic_error("Gradient function not defined.");
      vector<double> g(coords.size());
      energyGradient(coords, nullptr, &g);
      return g;
    }
  }


  void Potential::energyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (_energyGradientDef) {
      if (_energyGradient == nullptr) throw std::logic_error("Energy+gradient function marked as defined but not called.");
      return _energyGradient(coords, e, g);
    } else {
      if (!_energyDef) throw std::logic_error("Energy and/or gradient function not defined.");
      if (e != nullptr) *e = energy(coords);
      if (g != nullptr) *g = gradient(coords);
    }
  }


  void Potential::elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    throw std::logic_error("You shouldn't be here. The potential must override elementEnergyGradient if _parallelDef is marked as true.");
  }


  bool Potential::serialDef() const {
    return (_energyGradientDef || _energyDef);
  }


  bool Potential::parallelDef() const {
    return _parallelDef;
  }


  State Potential::newState(const vector<double>& coords, const vector<int>& ranks) {
    return State(*this, coords, ranks);
  }

  State Potential::newState(int ndof, const vector<int>& ranks) {
    vector<double> coords(ndof);
    return newState(coords, ranks);
  }


  Potential& Potential::setElements(vector<Element> elements)
  {
    this->elements = elements;
    return *this;
  }

  Potential& Potential::setElements(vector2d<int> idofs) {
    // Generate energy elements
    elements = {};
    for (const auto& idof: idofs) {
      elements.push_back({0, idof});
    }
    return *this;
  }

  Potential& Potential::setElements(vector2d<int> idofs, vector<int> types, vector2d<double> parameters) {
    // Generate energy elements
    elements = {};
    int nelements = idofs.size();
    for (int i=0; i<nelements; i++) {
      elements.push_back({types[i], idofs[i], parameters[i]});
    }
    return *this;
  }


  Potential& Potential::setConstraints(vector<int> iFix) {
    for (int i: iFix) {
      constraints.push_back({{i}});
    }
    return *this;
  }

  Potential& Potential::setConstraints(vector2d<int> idofs, vector<double> normal) {
    for (const vector<int>& idof: idofs) {
      constraints.push_back({idof, [normal](auto&&, auto&&){return normal;}});
    }
    return *this;
  }

  Potential& Potential::setConstraints(vector2d<int> idofs, Constraint::NormalFn normal, Constraint::CorrectionFn correction) {
    for (const vector<int>& idof: idofs) {
      constraints.push_back({idof, normal, correction});
    }
    return *this;
  }


  vector<double> Potential::applyConstraints(const vector<double>& coords, vector<double>& grad) const {
    for (const auto& constraint: constraints) {
      if (constraint.idof.size() == 1) grad[constraint.idof[0]] = 0;
      else {
        // Remove the component of grad in the direction of the normal
        auto normal = constraint.normal(coords);
        auto gradSlice = vec::slice(grad, constraint.idof);
        double normalSq = vec::dotProduct(normal, normal);
        double gradNormal = vec::dotProduct(gradSlice, normal);
        for (size_t i=0; i<normal.size(); i++) grad[constraint.idof[i]] -= gradNormal * normal[i] / normalSq;
      }
    }
    return grad;
  }

  vector<double> Potential::correctConstraints(vector<double>& coords) const {
    for (const auto& constraint: constraints) {
      if (constraint.correction) constraint.correction(constraint.idof, coords);
    }
    return coords;
  }


  bool Potential::isFixed(int index) const {
    for (const auto& constraint: constraints) {
      if (constraint.idof.size()!=1) continue;
      if (constraint.idof[0]==index) return true;
    }
    return false;
  }

  vector<bool> Potential::isFixed(const vector<int>& indicies) const {
    // Get array of fixed indicies
    vector<int> iFixed;
    iFixed.reserve(constraints.size());
    for (const auto& constraint: constraints) {
      if (constraint.idof.size()!=1) continue;
      iFixed.push_back(constraint.idof[0]);
    }
    // Find if the chosen indicies are included
    vector<bool> fixed(indicies.size(), false);
    for (int i=0; i<(int)indicies.size(); i++) {
      if (vec::isIn(iFixed, indicies[i])) fixed[i] = true;
    }
    return fixed;
  }


}
