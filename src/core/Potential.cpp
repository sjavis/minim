#include "Potential.h"

#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;

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


  double Potential::energy(const Vector& coords) const {
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


  Vector Potential::gradient(const Vector& coords) const {
    if (_energyDef) {
      if (_gradient == nullptr) throw std::logic_error("Gradient function marked as defined but not called.");
      return _gradient(coords);
    } else {
      if (!_energyGradientDef) throw std::logic_error("Gradient function not defined.");
      Vector g(coords.size());
      energyGradient(coords, nullptr, &g);
      return g;
    }
  }


  void Potential::energyGradient(const Vector& coords, double* e, Vector* g) const {
    if (_energyGradientDef) {
      if (_energyGradient == nullptr) throw std::logic_error("Energy+gradient function marked as defined but not called.");
      return _energyGradient(coords, e, g);
    } else {
      if (!_energyDef) throw std::logic_error("Energy and/or gradient function not defined.");
      if (e != nullptr) *e = energy(coords);
      if (g != nullptr) *g = gradient(coords);
    }
  }


  double Potential::blockEnergy(const Vector& coords) const {
    if (_blockEnergyDef) throw std::logic_error("Block energy function marked as defined but not called.");
    if (!_blockEnergyGradientDef) throw std::logic_error("Block energy function not defined.");
    double e;
    blockEnergyGradient(coords, &e, nullptr);
    return e;
  }


  Vector Potential::blockGradient(const Vector& coords) const {
    if (_blockEnergyDef) throw std::logic_error("Block gradient function marked as defined but not called.");
    if (!_blockEnergyGradientDef) throw std::logic_error("Block gradient function not defined.");
    Vector g(coords.size());
    blockEnergyGradient(coords, nullptr, &g);
    return g;
  }


  void Potential::blockEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (_blockEnergyGradientDef) throw std::logic_error("Block energy+gradient function marked as defined but not called.");
    if (!_blockEnergyDef) throw std::logic_error("Block energy and/or gradient function not defined.");
    if (e != nullptr) *e = blockEnergy(coords);
    if (g != nullptr) *g = blockGradient(coords);
  }


  bool totalEnergyDef() const {
    return (_energyGradientDef || _energyDef);
  }


  bool blockEnergyDef() const {
    return (_blockEnergyGradientDef || _blockEnergyDef);
  }


  State Potential::newState(const Vector& coords) {
    return State(*this, coords);
  }

  State Potential::newState(int ndof) {
    Vector coords(ndof);
    return newState(coords);
  }


  Potential& Potential::setElements(std::vector<Element> elements)
  {
    this->elements = elements;
    return *this;
  }

  Potential& Potential::setElements(std::vector<std::vector<int>> idofs) {
    // Generate energy elements
    elements = {};
    int id = 0;
    for (const auto& idof: idofs) {
      Element el = {id, 0, idof};
      elements.push_back(el);
      id++;
    }
    return *this;
  }

  Potential& Potential::setElements(std::vector<std::vector<int>> idofs, std::vector<int> types,
                                    std::vector<std::vector<double>> parameters)
  {
    // Generate energy elements
    elements = {};
    int id = 0;
    int nelements = idofs.size();
    for (int i=0; i<nelements; i++) {
      Element el = {id, types[i], idofs[i], parameters[i]};
      elements.push_back(el);
      id++;
    }
    return *this;
  }

}
