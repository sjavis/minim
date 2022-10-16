#include "Potential.h"

#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;


  double Potential::energy(const Vector& coords) const {
    return _energy(coords);
  }


  Vector Potential::gradient(const Vector& coords) const {
    return _gradient(coords);
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
