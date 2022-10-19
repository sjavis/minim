#include "Potential.h"

#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;


  double Potential::energy(const Vector& coords) const {
    if (energyDef) {
      return _energy(coords);
    } else if (energyGradientDef) {
      double e;
      energyGradient(coords, &e, nullptr);
      return e;
    } else {
      throw std::logic_error("Energy function not defined.");
    }
  }


  Vector Potential::gradient(const Vector& coords) const {
    if (gradientDef) {
      return _gradient(coords);
    } else if (energyGradientDef) {
      Vector g(coords.size());
      energyGradient(coords, nullptr, &g);
      return g;
    } else {
      throw std::logic_error("Gradient function not defined.");
    }
  }


  void Potential::energyGradient(const Vector& coords, double* e, Vector* g) const {
    if (energyGradientDef) {
      _energyGradient(coords, e, g);
    } else if (energyDef && gradientDef) {
      if (e != nullptr) *e = energy(coords);
      if (g != nullptr) *g = gradient(coords);
    } else {
      throw std::logic_error("Energy and/or gradient function not defined.");
    }
  }


  double Potential::blockEnergy(const Vector& coords) const {
    if (blockEnergyGradientDef) {
      double e;
      blockEnergyGradient(coords, &e, nullptr);
      return e;
    } else {
      throw std::logic_error("Energy function not defined.");
    }
  }


  Vector Potential::blockGradient(const Vector& coords) const {
    if (blockEnergyGradientDef) {
      Vector g(coords.size());
      blockEnergyGradient(coords, nullptr, &g);
      return g;
    } else {
      throw std::logic_error("Gradient function not defined.");
    }
  }


  void Potential::blockEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (blockEnergyDef && blockGradientDef) {
      if (e != nullptr) *e = blockEnergy(coords);
      if (g != nullptr) *g = blockGradient(coords);
    } else {
      throw std::logic_error("Energy and/or gradient function not defined.");
    }
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
