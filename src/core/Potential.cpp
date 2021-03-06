#include "Potential.h"

#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;


  Args::Args(int ndof, std::vector<std::vector<int>> idofs) : ndof(ndof) {
    // Generate energy elements
    int id = 0;
    for (int i=0; i<idofs.size(); i++) {
      Args::Element el = {id, 0, idofs[i]};
      elements.push_back(el);
      id++;
    }
  }

  Args::Args(int ndof, std::vector<std::vector<int>> idofs, std::vector<int> types,
             std::vector<std::vector<double>> parameters)
    : ndof(ndof) {
    // Generate energy elements
    int id = 0;
    for (int i=0; i<idofs.size(); i++) {
      Args::Element el = {id, types[i], idofs[i], parameters[i]};
      elements.push_back(el);
      id++;
    }
  }


  double Potential::energy(const Vector &coords, const Args &args) {
    return _energy(coords, args);
  }


  Vector Potential::gradient(const Vector &coords, const Args &args) {
    return _gradient(coords, args);
  }


  Args* Potential::newArgs(int ndof) {
    return new Args(ndof);
  }


  State Potential::newState(Vector coords) {
    Args *args = newArgs(coords.size());
    return State(*this, coords, *args);
  }

  State Potential::newState(int ndof) {
    Vector coords(ndof);
    return newState(coords);
  }

}
