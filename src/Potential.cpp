#include "Potential.h"

#include "State.h"
#include "communicators/CommUnstructured.h"
#include "communicators/CommGrid.h"
#include "utils/vec.h"
#include <stdexcept>

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  Potential::Potential(EFunc energy, GFunc gradient)
    : _energy(energy), _gradient(gradient), _energyGradient(nullptr)
  {}

  Potential::Potential(EGFunc energyGradient)
    : _energy(nullptr), _gradient(nullptr), _energyGradient(energyGradient)
  {}


  double Potential::energy(const vector<double>& coords) const {
    if (_energy) return _energy(coords);
    throw std::runtime_error("Potential: No energy function defined. If using parallel you must call energyGradient instead.");
  }


  vector<double> Potential::gradient(const vector<double>& coords) const {
    if (_gradient) return _gradient(coords);
    throw std::runtime_error("Potential: No gradient function defined. If using parallel you must call energyGradient instead.");
  }


  void Potential::energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {
    if (_energyGradient) return _energyGradient(coords, e, g);

    if (e != nullptr) *e = energy(coords);
    if (g != nullptr) *g = gradient(coords);
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

  vector<char> Potential::isFixed(const vector<int>& indicies) const {
    // Get array of fixed indicies
    vector<int> iFixed;
    iFixed.reserve(constraints.size());
    for (const auto& constraint: constraints) {
      if (constraint.idof.size()!=1) continue;
      iFixed.push_back(constraint.idof[0]);
    }
    // Find if the chosen indicies are included
    vector<char> fixed(indicies.size(), false);
    for (int i=0; i<(int)indicies.size(); i++) {
      if (vec::isIn(iFixed, indicies[i])) fixed[i] = true;
    }
    return fixed;
  }


  bool Potential::isSerial() const {
    return potentialType() == SERIAL;
  }


  std::unique_ptr<Communicator> Potential::newComm() const {
    if (potentialType() == UNSTRUCTURED) {
      if (this->distributed) throw std::invalid_argument("Potential: You cannot create a State with a Potential that has already been distributed.");
      return std::make_unique<CommUnstructured>();
    } else if (potentialType() == GRID) {
      return std::make_unique<CommGrid>();
    }
    return std::make_unique<CommUnstructured>();
  }


}
