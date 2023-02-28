#include "BarAndHinge.h"

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
    // Check element DOFs are valid
    int nNode = coords.size() / 3;
    vector2d<int> elementList = bondList;
    elementList.insert(elementList.end(), hingeList.begin(), hingeList.end());
    for (auto element: elementList) {
      for (int iNode: element) {
        if (iNode >= nNode) throw std::logic_error("The triangulation expects more coordinates than are given.");
      }
    }
    if (!fixed.empty() && fixed.size()!=coords.size()) {
      throw std::logic_error("The vector for the fixed degrees of freedom has the wrong size");
    }

    // Get rigidities
    if (thickness.size() == 1) {
      double t = thickness[0];
      thickness = vector<double>(coords.size()/3, t);
    }
    if (kBond.size() == 1) {
      double k = kBond[0];
      kBond = vector<double>(bondList.size(), k);
    } else if (kBond.empty() && !thickness.empty()) {
      kBond.reserve(bondList.size());
      for (auto bond: bondList) {
        double t = (thickness[bond[0]] + thickness[bond[1]]) / 2;
        double k = sqrt(3)/4 * modulus * t;
        kBond.push_back(k);
      }
    }
    if (kHinge.size() == 1) {
      double k = kHinge[0];
      kHinge = vector<double>(hingeList.size(), k);
    } else if (kHinge.empty() && !thickness.empty()) {
      kHinge.reserve(hingeList.size());
      for (auto hinge: hingeList) {
        double t = (thickness[hinge[0]] + thickness[hinge[1]] + thickness[hinge[2]] + thickness[hinge[3]]) / 4;
        vector<double> x1 = {coords[3*hinge[0]], coords[3*hinge[0]+1], coords[3*hinge[0]+2]};
        vector<double> x2 = {coords[3*hinge[1]], coords[3*hinge[1]+1], coords[3*hinge[1]+2]};
        vector<double> x3 = {coords[3*hinge[2]], coords[3*hinge[2]+1], coords[3*hinge[2]+2]};
        vector<double> x4 = {coords[3*hinge[3]], coords[3*hinge[3]+1], coords[3*hinge[3]+2]};
        double lengthSq = vec::sum(vec::pow(x3-x2, 2));
        double area1 = std::abs(vec::norm(vec::crossProduct(x2-x1, x3-x1))) / 2;
        double area2 = std::abs(vec::norm(vec::crossProduct(x2-x4, x3-x4))) / 2;
        double areaSum = area1 + area2;
        double k = modulus * pow(t,3) / (12 * (1 - pow(poissonRatio,2)));
        k *= lengthSq / areaSum;
        kHinge.push_back(k);
      }
    }

    // Get equilibrium values
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
      elements.push_back({1, idofs, {kHinge[iH], theta0[iH]}});
    }
    for (int iN=0; iN<(int)coords.size()/3; iN++) {
      if (!fixed.empty() && fixed[3*iN]) continue;
      // External force
      elements.push_back({2, {3*iN, 3*iN+1, 3*iN+2}, {(double)iN}});
      // Substrate interaction
      if (wallOn) elements.push_back({3, {3*iN, 3*iN+1, 3*iN+2}});
    }
  }


  void BarAndHinge::distributeParameters(const Communicator& comm) {
    if (!fixed.empty()) fixed = comm.assignProc(fixed);
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
      *e += el.parameters[0] * pow(dl, 2);
    }
    if (g != nullptr) {
      auto g_factor = 2 * el.parameters[0] * dl / l;
      if (fixed.empty() || !fixed[el.idof[0]]) (*g)[el.idof[0]] += g_factor * dx[0];
      if (fixed.empty() || !fixed[el.idof[1]]) (*g)[el.idof[1]] += g_factor * dx[1];
      if (fixed.empty() || !fixed[el.idof[2]]) (*g)[el.idof[2]] += g_factor * dx[2];
      if (fixed.empty() || !fixed[el.idof[3]]) (*g)[el.idof[3]] -= g_factor * dx[0];
      if (fixed.empty() || !fixed[el.idof[4]]) (*g)[el.idof[4]] -= g_factor * dx[1];
      if (fixed.empty() || !fixed[el.idof[5]]) (*g)[el.idof[5]] -= g_factor * dx[2];
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
      *e += el.parameters[0] * 0.5 * pow(dtheta, 2);
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
      if (fixed.empty() || !fixed[el.idof[0]])  (*g)[el.idof[0]] += -gfactor * n1h[0];
      if (fixed.empty() || !fixed[el.idof[1]])  (*g)[el.idof[1]] += -gfactor * n1h[1];
      if (fixed.empty() || !fixed[el.idof[2]])  (*g)[el.idof[2]] += -gfactor * n1h[2];
      if (fixed.empty() || !fixed[el.idof[3]])  (*g)[el.idof[3]] +=  gfactor * ((1-skew1)*n1h[0] - skew2*n2h[0]);
      if (fixed.empty() || !fixed[el.idof[4]])  (*g)[el.idof[4]] +=  gfactor * ((1-skew1)*n1h[1] - skew2*n2h[1]);
      if (fixed.empty() || !fixed[el.idof[5]])  (*g)[el.idof[5]] +=  gfactor * ((1-skew1)*n1h[2] - skew2*n2h[2]);
      if (fixed.empty() || !fixed[el.idof[6]])  (*g)[el.idof[6]] += -gfactor * ((1-skew2)*n2h[0] - skew1*n1h[0]);
      if (fixed.empty() || !fixed[el.idof[7]])  (*g)[el.idof[7]] += -gfactor * ((1-skew2)*n2h[1] - skew1*n1h[1]);
      if (fixed.empty() || !fixed[el.idof[8]])  (*g)[el.idof[8]] += -gfactor * ((1-skew2)*n2h[2] - skew1*n1h[2]);
      if (fixed.empty() || !fixed[el.idof[9]])  (*g)[el.idof[9]]  += gfactor * n2h[0];
      if (fixed.empty() || !fixed[el.idof[10]]) (*g)[el.idof[10]] += gfactor * n2h[1];
      if (fixed.empty() || !fixed[el.idof[11]]) (*g)[el.idof[11]] += gfactor * n2h[2];
    }
  }


  void BarAndHinge::forceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    if (!fixed.empty() && fixed[el.idof[2]]) return;
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
    if (!fixed.empty() && fixed[el.idof[2]]) return;
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
    // BondList
    bondList.reserve(3*triList.size());
    for (auto tri: triList) {
      std::sort(tri.begin(), tri.end());
      bondList.push_back({tri[0], tri[1]});
      bondList.push_back({tri[0], tri[2]});
      bondList.push_back({tri[1], tri[2]});
    }
    std::sort(bondList.begin(), bondList.end());
    bondList.erase(std::unique(bondList.begin(), bondList.end()), bondList.end());
    // HingeList
    hingeList.reserve(3*triList.size());
    vector<bool> shared(6, false);
    for (int it1=0; it1<(int)triList.size(); it1++) {
      for (int it2=it1+1; it2<(int)triList.size(); it2++) {
        int nShared = 0;
        for (int in1=0; in1<3; in1++) {
          for (int in2=0; in2<3; in2++) {
            if (triList[it1][in1] != triList[it2][in2]) continue;
            nShared++;
            shared[in1] = true;
            shared[3+in2] = true;
          }
        }
        if (nShared == 0) {
          continue;
        } else if (nShared != 2) {
          shared = vector<bool>(6, false);
          continue;
        }
        // Keep the cyclic order of the first triangle
        // This allows the user to use theta0 in the full 2pi range
        vector<int> nodes(4);
        for (int in1=0; in1<3; in1++) {
          if (!shared[in1]) {
            nodes[0] = triList[it1][in1];
            nodes[1] = triList[it1][(in1+1)%3];
            nodes[2] = triList[it1][(in1+2)%3];
          }
        }
        for (int in2=0; in2<3; in2++) {
          if (!shared[3+in2]) nodes[3] = triList[it2][in2];
        }
        shared = vector<bool>(6, false);
        hingeList.push_back(nodes);
      }
    }
    return *this;
  }

  BarAndHinge& BarAndHinge::setBondList(const vector2d<int>& bondList) {
    this->bondList = bondList;
    return *this;
  }

  BarAndHinge& BarAndHinge::setHingeList(const vector2d<int>& hingeList) {
    this->hingeList = hingeList;
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
    this->kBond = kBond;
    this->kHinge = kHinge;
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

  BarAndHinge& BarAndHinge::setFixed(const vector<bool>& fixed) {
    this->fixed = fixed;
    return *this;
  }

}
