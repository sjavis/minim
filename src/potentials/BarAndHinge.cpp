#include "potentials/BarAndHinge.h"

#include <math.h>
#include <stdexcept>
#include <algorithm>
#include "Communicator.h"
#include "utils/vec.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;

  constexpr double pi() { return 4*atan(1); }


  void BarAndHinge::init(const vector<double>& coords) {
    // Get the information required for the elements
    vector2d<int> bondList = computeBondList();
    vector2d<int> hingeList = computeHingeList(bondList);
    computeRigidities(coords, bondList, hingeList);
    computeLength0(coords, bondList);
    computeTheta0(coords, hingeList);
    vector<char> fixed = isFixed(vec::iota(coords.size()));

    // Check element DOFs are valid TODO: fix this
    // int nNode = coords.size() / 3;
    // vector2d<int> elementList = bondList;
    // elementList.insert(elementList.end(), hingeList.begin(), hingeList.end());
    // for (auto element: elementList) {
    //   for (int iNode: element) {
    //     if (iNode >= nNode) throw std::logic_error("The triangulation expects more coordinates than are given.");
    //   }
    // }

    // Assign elements
    elements = {};
    for (int iB=0; iB<(int)bondList.size(); iB++) {
      vector<int> idofs(6);
      for (int iN=0; iN<2; iN++) {
        auto n = bondList[iB][iN];
        idofs[3*iN+0] = 3 * n;
        idofs[3*iN+1] = 3 * n + 1;
        idofs[3*iN+2] = 3 * n + 2;
      }
      if (vec::all(vec::slice(fixed, idofs))) continue;
      elements.push_back({0, idofs, {kBond[iB], length0[iB]}});
    }
    for (int iH=0; iH<(int)hingeList.size(); iH++) {
      vector<int> idofs(12);
      for (int iN=0; iN<4; iN++) {
        auto n = hingeList[iH][iN];
        idofs[3*iN+0] = 3 * n;
        idofs[3*iN+1] = 3 * n + 1;
        idofs[3*iN+2] = 3 * n + 2;
      }
      if (vec::all(vec::slice(fixed, idofs))) continue;
      elements.push_back({1, idofs, {kHinge[iH], theta0[iH]}});
    }
    for (int iN=0; iN<(int)coords.size()/3; iN++) {
      vector<int> idofs{3*iN, 3*iN+1, 3*iN+2};
      if (vec::all(vec::slice(fixed, idofs))) continue;
      // External force
      elements.push_back({2, idofs, {(double)iN}});
      // Substrate interaction
      if (wallOn) elements.push_back({3, idofs});
    }
  }


  vector2d<int> BarAndHinge::computeBondList() {
    if (!_bondList.empty()) return _bondList;
    if (_triList.empty()) throw std::invalid_argument("BarAndHinge: The triangulation has not been set.");
    vector2d<int> bondList;
    bondList.reserve(3*_triList.size());
    // Get the list of pairs
    for (auto tri: _triList) {
      std::sort(tri.begin(), tri.end());
      bondList.push_back({tri[0], tri[1]});
      bondList.push_back({tri[0], tri[2]});
      bondList.push_back({tri[1], tri[2]});
    }
    bondList = vec::unique(bondList);
    // Add the additional indicies for each triangle
    for (auto& bond: bondList) {
      for (auto tri: _triList) {
        if (!vec::isIn(tri, bond[0]) || !vec::isIn(tri, bond[1])) continue;
        for (int iN: tri) {
          if (!vec::isIn(bond, iN)) bond.push_back(iN);
        }
      }
    }
    return bondList;
  }


  vector2d<int> BarAndHinge::computeHingeList(const vector2d<int>& bondList) {
    if (!_hingeList.empty()) return _hingeList;
    if (_triList.empty()) throw std::invalid_argument("BarAndHinge: The triangulation has not been set.");
    vector2d<int> hingeList;
    hingeList.reserve(3*_triList.size());
    vector<char> shared(6, false);
    for (int it1=0; it1<(int)_triList.size(); it1++) {
      for (int it2=it1+1; it2<(int)_triList.size(); it2++) {
        int nShared = 0;
        for (int in1=0; in1<3; in1++) {
          for (int in2=0; in2<3; in2++) {
            if (_triList[it1][in1] != _triList[it2][in2]) continue;
            nShared++;
            shared[in1] = true;
            shared[3+in2] = true;
          }
        }
        if (nShared == 0) {
          continue;
        } else if (nShared != 2) {
          shared = vector<char>(6, false);
          continue;
        }
        // Keep the cyclic order of the first triangle
        // This allows the user to use theta0 in the full 2pi range
        vector<int> nodes(4);
        for (int in1=0; in1<3; in1++) {
          if (!shared[in1]) {
            nodes[0] = _triList[it1][in1];
            nodes[1] = _triList[it1][(in1+1)%3];
            nodes[2] = _triList[it1][(in1+2)%3];
          }
        }
        for (int in2=0; in2<3; in2++) {
          if (!shared[3+in2]) nodes[3] = _triList[it2][in2];
        }
        shared = vector<char>(6, false);
        hingeList.push_back(nodes);
      }
    }
    // Append the bond number of each hinge
    for (auto& hinge: hingeList) {
      for (int iBond=0; iBond<(int)bondList.size(); iBond++) {
        vector<int> bond = {bondList[iBond][0], bondList[iBond][1]};
        if (vec::isIn(bond, hinge[1]) && vec::isIn(bond, hinge[2])) {
          hinge.push_back(iBond);
          break;
        }
      }
    }
    return hingeList;
  }


  void BarAndHinge::computeRigidities(const vector<double>& coords, const vector2d<int>& bondList, const vector2d<int>& hingeList) {
    // All rigidities the same
    if (kBond.size() == 1) {
      double k = kBond[0];
      kBond = vector<double>(bondList.size(), k);
    }
    if (kHinge.size() == 1) {
      double k = kHinge[0];
      kHinge = vector<double>(hingeList.size(), k);
    }
    if (!kBond.empty() && !kHinge.empty()) return;

    // Ensure thickness is set
    if (thickness.empty()) {
      throw std::invalid_argument("BarAndHinge requires either the rigidities or the thickness to be set");
    } else if (thickness.size() == 1) {
      double t = thickness[0];
      thickness = vector<double>(coords.size()/3, t);
    }

    // Calculate the length and area sum of each bond
    int nBond = bondList.size();
    vector<double> lengthSq(nBond);
    vector<double> areaSum(nBond);
    for (int iBond=0; iBond<nBond; iBond++) {
      vector<int> bond = bondList[iBond];
      if (bond.size()==2) throw std::invalid_argument("The bondList should contain the indicies of the adjacent triangles.");
      vector<double> x1 = {coords[3*bond[0]], coords[3*bond[0]+1], coords[3*bond[0]+2]};
      vector<double> x2 = {coords[3*bond[1]], coords[3*bond[1]+1], coords[3*bond[1]+2]};
      lengthSq[iBond] = vec::sum(vec::pow(x2-x1, 2));
      // Triangle area
      for (int i3: vector<int>(bond.begin()+2, bond.end())) {
        vector<double> x3 = {coords[3*i3], coords[3*i3+1], coords[3*i3+2]};
        areaSum[iBond] += std::abs(vec::norm(vec::crossProduct(x2-x1, x3-x1))) / 2;
      }
    }

    // Stretching rigidities
    if (kBond.empty()) {
      kBond.reserve(nBond);
      for (int iBond=0; iBond<nBond; iBond++) {
        vector<int> bond = bondList[iBond];
        double t = (thickness[bond[0]] + thickness[bond[1]]) / 2;
        double rigidity = modulus * t;
        double factor = areaSum[iBond] / lengthSq[iBond];
        kBond.push_back(rigidity * factor);
      }
    }

    // Bending rigidities
    if (kHinge.empty()) {
      kHinge.reserve(hingeList.size());
      for (auto hinge: hingeList) {
        if (hinge.size()==4) throw std::invalid_argument("The hingeList should contain the index of the bond.");
        double t = (thickness[hinge[0]] + thickness[hinge[1]] + thickness[hinge[2]] + thickness[hinge[3]]) / 4;
        double rigidity = modulus * pow(t,3) / (12 * (1 - pow(poissonRatio,2)));
        double factor = lengthSq[hinge[4]] / areaSum[hinge[4]];
        kHinge.push_back(rigidity * factor);
      }
    }
  }


  void BarAndHinge::computeLength0(const vector<double>& coords, const vector2d<int>& bondList) {
    if (length0.empty()) {
      length0 = vector<double>(bondList.size());
      for (int iB=0; iB<(int)bondList.size(); iB++) {
        auto n = bondList[iB];
        vector<double> x1 = {coords[3*n[0]], coords[3*n[0]+1], coords[3*n[0]+2]};
        vector<double> x2 = {coords[3*n[1]], coords[3*n[1]+1], coords[3*n[1]+2]};
        length0[iB] = vec::norm(x1 - x2);
      }
    } else if (length0.size() == 1) {
      double l0 = length0[0];
      length0 = vector<double>(bondList.size(), l0);
    }
  }


  void BarAndHinge::computeTheta0(const vector<double>& coords, const vector2d<int>& hingeList) {
    if (theta0.empty()) {
      theta0 = vector<double>(hingeList.size());
      for (int iH=0; iH<(int)hingeList.size(); iH++) {
        auto n = hingeList[iH];
        vector<double> x1 = {coords[3*n[0]], coords[3*n[0]+1], coords[3*n[0]+2]};
        vector<double> x2 = {coords[3*n[1]], coords[3*n[1]+1], coords[3*n[1]+2]};
        vector<double> x3 = {coords[3*n[2]], coords[3*n[2]+1], coords[3*n[2]+2]};
        vector<double> x4 = {coords[3*n[3]], coords[3*n[3]+1], coords[3*n[3]+2]};
        auto n1 = vec::crossProduct(x2-x1, x3-x2);
        auto n2 = vec::crossProduct(x3-x2, x4-x3);
        double c = vec::dotProduct(n1, n2) / (vec::norm(n1) * vec::norm(n2));
        double n1b3 = vec::dotProduct(n1, x4-x3);
        theta0[iH] = (n1b3>=0) ? acos(c) : 2*pi()-acos(c);
      }
    } else if (theta0.size() == 1) {
      double t0 = theta0[0];
      theta0 = vector<double>(hingeList.size(), t0);
    }
  }


  void BarAndHinge::elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    switch (el.type) {
      case 0:
        stretching(coords, el, e, g);
        break;
      case 1:
        bending(coords, el, e, g);
        break;
      case 2:
        forceEnergy(coords, el, e, g);
        break;
      case 3:
        substrate(coords, el, e, g);
        break;
    }
  }


  void BarAndHinge::stretching(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    vector<double> x1{coords[el.idof[0]], coords[el.idof[1]], coords[el.idof[2]]};
    vector<double> x2{coords[el.idof[3]], coords[el.idof[4]], coords[el.idof[5]]};

    // Compute distance
    auto dx = x1 - x2;
    auto l = vec::norm(dx);
    auto dl = l - el.parameters[1];

    if (e != nullptr) {
      *e += 0.5 * el.parameters[0] * pow(dl, 2);
    }
    if (g != nullptr) {
      auto g_factor = 2 * el.parameters[0] * dl / l;
      (*g)[el.idof[0]] += g_factor * dx[0];
      (*g)[el.idof[1]] += g_factor * dx[1];
      (*g)[el.idof[2]] += g_factor * dx[2];
      (*g)[el.idof[3]] -= g_factor * dx[0];
      (*g)[el.idof[4]] -= g_factor * dx[1];
      (*g)[el.idof[5]] -= g_factor * dx[2];
    }
  }


  void BarAndHinge::bending(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    vector<double> x1{coords[el.idof[0]], coords[el.idof[1]], coords[el.idof[2]]};
    vector<double> x2{coords[el.idof[3]], coords[el.idof[4]], coords[el.idof[5]]};
    vector<double> x3{coords[el.idof[6]], coords[el.idof[7]], coords[el.idof[8]]};
    vector<double> x4{coords[el.idof[9]], coords[el.idof[10]], coords[el.idof[11]]};

    // Compute bond vectors
    auto b1 = x2 - x1;
    auto b2 = x3 - x2;
    auto b3 = x4 - x3;
    double b2m = vec::norm(b2);

    // Compute normal vectors
    auto n1 = vec::crossProduct(b1, b2);
    auto n2 = vec::crossProduct(b2, b3);
    double n1sq = vec::dotProduct(n1, n1);
    double n2sq = vec::dotProduct(n2, n2);
    double n12m = sqrt(n1sq * n2sq);

    // Compute cos
    double c = vec::dotProduct(n1, n2) / n12m;
    c = std::max(c, -1.0);
    c = std::min(c,  1.0);
    // E = k (1 - cos(\theta - \theta_0))
    // double s = b2m / n12m * vec::dotProduct(n1, b3);
    // double c0 = cos(el.parameters[1]);
    // double s0 = sin(el.parameters[1]);
    // E = k/2 (\theta - \theta_0)^2
    double n1b3 = vec::dotProduct(n1, b3);
    double theta = (n1b3>=0) ? acos(c) : 2*pi()-acos(c);
    double dtheta = theta - el.parameters[1];
    dtheta = fmod(dtheta + pi(), 2*pi()) - pi();

    if (e != nullptr) {
      // E = k (1 - cos(\theta - \theta_0))
      // *e += el.parameters[0] * (1 - c*c0 + s*s0); // Double angle formula for cos(t - t0)
      // E = k/2 (\theta - \theta_0)^2
      *e += 0.5 * el.parameters[0] * pow(dtheta, 2);
    }
    if (g != nullptr) {
      // E = k (1 - cos(\theta - \theta_0))
      // double gfactor = el.parameters[0] * (s*c0 + c*s0);
      // E = k/2 (\theta - \theta_0)^2
      double gfactor = el.parameters[0] * (dtheta);
      // Normal vectors with magnitude 1 / triangle height
      auto n1h = b2m/n1sq * n1;
      auto n2h = b2m/n2sq * n2;
      // Quantify triangular skew, 0.5 if symmetric
      double skew1 = -vec::dotProduct(b1, b2) / (b2m*b2m);
      double skew2 = -vec::dotProduct(b3, b2) / (b2m*b2m);
      (*g)[el.idof[0]] += -gfactor * n1h[0];
      (*g)[el.idof[1]] += -gfactor * n1h[1];
      (*g)[el.idof[2]] += -gfactor * n1h[2];
      (*g)[el.idof[3]] +=  gfactor * ((1-skew1)*n1h[0] - skew2*n2h[0]);
      (*g)[el.idof[4]] +=  gfactor * ((1-skew1)*n1h[1] - skew2*n2h[1]);
      (*g)[el.idof[5]] +=  gfactor * ((1-skew1)*n1h[2] - skew2*n2h[2]);
      (*g)[el.idof[6]] += -gfactor * ((1-skew2)*n2h[0] - skew1*n1h[0]);
      (*g)[el.idof[7]] += -gfactor * ((1-skew2)*n2h[1] - skew1*n1h[1]);
      (*g)[el.idof[8]] += -gfactor * ((1-skew2)*n2h[2] - skew1*n1h[2]);
      (*g)[el.idof[9]]  += gfactor * n2h[0];
      (*g)[el.idof[10]] += gfactor * n2h[1];
      (*g)[el.idof[11]] += gfactor * n2h[2];
    }
  }


  void BarAndHinge::forceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    vector<double> force(3);
    if (this->force.size() == 3) {
      force = this->force;
    } else if (!this->force.empty()) {
      int iN = (int) el.parameters[0];
      force = {this->force[3*iN], this->force[3*iN+1], this->force[3*iN+2]};
    }
    if (e != nullptr) (*e) -= force[2] * coords[el.idof[2]];
    if (g != nullptr) (*g)[el.idof[2]] -= force[2];
  }


  void BarAndHinge::substrate(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    // Height
    double height = coords[el.idof[2]];
    if (!wallAdhesion && height>0) return;
    double lj_r0 = 0.858374218932559*lj_sigma;
    height = std::max(height + lj_r0, lj_r0/2);
    // Lennard-Jones potential
    double lj9 = 2.0/15 * pow(lj_sigma / height, 9);
    double lj3 = pow(lj_sigma / height, 3);
    double sigma3 = pow(lj_sigma, 3);
    double sigma9 = pow(sigma3, 3);
    double lj9_0 = 0.5270462766947298*sigma9;
    double lj3_0 = 1.5811388300841895*sigma3;
    if (e != nullptr) (*e) += lj_epsilon * ((lj9-lj9_0) - (lj3-lj3_0));
    if (g != nullptr) (*g)[el.idof[2]] -= lj_epsilon * (9*lj9 - 3*lj3) / height;
  }


  BarAndHinge& BarAndHinge::setTriangulation(const vector2d<int>& triList) {
    this->_triList = triList;
    return *this;
  }

  BarAndHinge& BarAndHinge::setBondList(const vector2d<int>& bondList) {
    this->_bondList = bondList;
    return *this;
  }

  BarAndHinge& BarAndHinge::setHingeList(const vector2d<int>& hingeList) {
    this->_hingeList = hingeList;
    return *this;
  }

  BarAndHinge& BarAndHinge::setModulus(double modulus) {
    this->modulus = modulus;
    return *this;
  }

  BarAndHinge& BarAndHinge::setThickness(double thickness) {
    this->thickness = {thickness};
    return *this;
  }

  BarAndHinge& BarAndHinge::setThickness(const vector<double>& thickness) {
    this->thickness = thickness;
    return *this;
  }

  BarAndHinge& BarAndHinge::setRigidity(double kBond, double kHinge) {
    this->kBond = {kBond};
    this->kHinge = {kHinge};
    return *this;
  }

  BarAndHinge& BarAndHinge::setRigidity(const vector<double>& kBond, const vector<double>& kHinge) {
    if (!kBond.empty()) this->kBond = kBond;
    if (!kHinge.empty()) this->kHinge = kHinge;
    return *this;
  }

  BarAndHinge& BarAndHinge::setLength0(double length0) {
    this->length0 = {length0};
    return *this;
  }

  BarAndHinge& BarAndHinge::setLength0(const vector<double>& length0) {
    this->length0 = length0;
    return *this;
  }

  BarAndHinge& BarAndHinge::setTheta0(double theta0) {
    this->theta0 = {theta0};
    return *this;
  }

  BarAndHinge& BarAndHinge::setTheta0(const vector<double>& theta0) {
    this->theta0 = theta0;
    return *this;
  }

  BarAndHinge& BarAndHinge::setWall(bool wallOn) {
    this->wallOn = wallOn;
    return *this;
  }

  BarAndHinge& BarAndHinge::setWallAdhesion(bool wallAdhesion) {
    this->wallAdhesion = wallAdhesion;
    return *this;
  }

  BarAndHinge& BarAndHinge::setWallParams(double epsilon, double sigma) {
    this->lj_epsilon = epsilon;
    this->lj_sigma = sigma;
    return *this;
  }

  BarAndHinge& BarAndHinge::setForce(const vector<double>& force) {
    this->force = force;
    return *this;
  }

  BarAndHinge& BarAndHinge::setForce(const vector2d<double>& force) {
    this->force.reserve(3*force.size());
    for (int i=0; i<(int)force.size(); i++) {
      for (int j=0; j<3; j++) {
        this->force.push_back(force[i][j]);
      }
    }
    return *this;
  }

}
