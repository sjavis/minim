#include "Potential.h"

#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;
  typedef Potential::Args Args;


  double Potential::energy(const Vector& coords, const Args& args) const {
    return _energy(coords, args);
  }


  Vector Potential::gradient(const Vector& coords, const Args& args) const {
    return _gradient(coords, args);
  }


  std::unique_ptr<Args> Potential::newArgs(int ndof) {
    return std::make_unique<Args>(Args(ndof));
  }


  State Potential::newState(Vector coords) {
    std::unique_ptr<Args> args = newArgs(coords.size());
    std::unique_ptr<Potential> pot = this->clone();
    return State(pot, coords, args);
  }

  State Potential::newState(int ndof) {
    Vector coords(ndof);
    return newState(coords);
  }


  // Args
  Args::Args(int ndof, std::vector<std::vector<int>> idofs) : ndof(ndof) {
    // Generate energy elements
    int id = 0;
    for (const auto& idof: idofs) {
      Args::Element el = {id, 0, idof};
      elements.push_back(el);
      id++;
    }
  }

  Args::Args(int ndof, std::vector<std::vector<int>> idofs, std::vector<int> types,
             std::vector<std::vector<double>> parameters)
    : ndof(ndof) {
    // Generate energy elements
    int id = 0;
    int nelements = idofs.size();
    for (int i=0; i<nelements; i++) {
      Args::Element el = {id, types[i], idofs[i], parameters[i]};
      elements.push_back(el);
      id++;
    }
  }
  

  std::unique_ptr<Args> Args::clone() const {
    return std::unique_ptr<Args>(new Args(*this));
  }

}
