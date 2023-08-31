#include "PFWetting.h"

#include <math.h>
#include <stdexcept>
#include <functional>
#include "State.h"
#include "utils/vec.h"


namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  std::array<int,3> getCoord(int i, std::array<int,3> gridSize) {
    int z = i % gridSize[2];
    int y = (i - z) / gridSize[2] % gridSize[1];
    int x = i / (gridSize[1] * gridSize[2]);
    return {x, y, z};
  }

  // Used to calculate the indicies of the neighbouring nodes
  class Neighbours {
    public:
      std::array<int,26> di;
      vector2d<int> dx = {{
        // Adjacent faces
        {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, { 0,  0,  1}, { 0,  1,  0}, { 1,  0,  0},
        // Adjacent edges
        {-1, -1,  0}, {-1,  0, -1}, {-1,  0,  1}, {-1,  1,  0},
        { 0, -1, -1}, { 0, -1,  1}, { 0,  1, -1}, { 0,  1,  1},
        { 1, -1,  0}, { 1,  0, -1}, { 1,  0,  1}, { 1,  1,  0},
        // Adjacent corners
        {-1, -1, -1}, {-1, -1,  1}, {-1,  1, -1}, {-1,  1,  1},
        { 1, -1, -1}, { 1, -1,  1}, { 1,  1, -1}, { 1,  1,  1}
      }};

      int operator[](int i) { return di[i]; };

      Neighbours(std::array<int,3> gridSize, int i0) {
        auto x0 = getCoord(i0, gridSize);
        for (int i=0; i<26; i++) {
          int x = (x0[0] + dx[i][0] + gridSize[0]) % gridSize[0];
          int y = (x0[1] + dx[i][1] + gridSize[1]) % gridSize[1];
          int z = (x0[2] + dx[i][2] + gridSize[2]) % gridSize[2];
          di[i] = (x*gridSize[1] + y)*gridSize[2] + z;
        }
      }
  };


  void PFWetting::assignKappa() {
    if ((int)interfaceSize.size()==1 && nFluid>1) {
      interfaceSize = vector<double>(nFluid, interfaceSize[0]);
    } else if ((int)interfaceSize.size() != nFluid) {
      throw std::invalid_argument("Invalid size of interfaceSize array.");
    }
    if ((int)surfaceTension.size()==1 && nFluid>1) {
      surfaceTension = vector<double>(nFluid, surfaceTension[0]);
    } else if ((int)surfaceTension.size() != nFluid) {
      throw std::invalid_argument("Invalid size of surfaceTension array.");
    }
    if (nFluid <= 2) {
      kappa = vector<double>(nFluid, 3*surfaceTension[0]/interfaceSize[0]);
      kappaP = vector<double>(nFluid, pow(interfaceSize[0], 2)*kappa[0]);
    } else {
      kappa = vector<double>(nFluid);
      kappaP = vector<double>(nFluid);
      auto kappaSums = 6 * surfaceTension / interfaceSize;
      auto kappaPSums = 6 * surfaceTension * interfaceSize;
      kappa[0] = 0.5 * ( kappaSums[0] + kappaSums[1] - kappaSums[nFluid-1]);
      kappa[1] = 0.5 * ( kappaSums[0] - kappaSums[1] + kappaSums[nFluid-1]);
      kappa[2] = 0.5 * (-kappaSums[0] + kappaSums[1] + kappaSums[nFluid-1]);
      kappaP[0] = 0.5 * ( kappaPSums[0] + kappaPSums[1] - kappaPSums[nFluid-1]);
      kappaP[1] = 0.5 * ( kappaPSums[0] - kappaPSums[1] + kappaPSums[nFluid-1]);
      kappaP[2] = 0.5 * (-kappaPSums[0] + kappaPSums[1] + kappaPSums[nFluid-1]);
      for (int i=3; i<nFluid; i++) {
        kappa[i] = 0.5 * (-kappaSums[0] - kappaSums[1] + kappaSums[nFluid-1]) + kappaSums[i-1];
        kappaP[i] = 0.5 * (-kappaPSums[0] - kappaPSums[1] + kappaPSums[nFluid-1]) + kappaPSums[i-1];
      }
    }
  }


  void PFWetting::init(const vector<double>& coords) {
    if ((int)coords.size() != nGrid*nFluid) {
      throw std::invalid_argument("Size of coordinates array does not match the grid size.");
    }

    setDefaults();
    checkArraySizes();
    assignKappa();

    // Forces
    vector<double> fMag(nFluid);
    vector2d<double> fNorm(nFluid);
    if (!force.empty()) {
      for (int i=0; i<nFluid; i++) {
        fMag[i] = vec::norm(force[i]);
        fNorm[i] = force[i] / fMag[i];
      }
    }

    // Get fluid volume and solid surface area for each node
    nodeVol = vector<double>(nGrid, 0);
    auto surfaceArea = vector<double>(nGrid, 0);
    for (int iGrid=0; iGrid<nGrid; iGrid++) {
      int type = getType(iGrid);
      if (type == 0) {
        nodeVol[iGrid] = 1;
      } else if (type > 0) {
        nodeVol[iGrid] = type / 8.0;
        if (type == 2 || type == 4 || type == 6) surfaceArea[iGrid] = 1;
        if (type == 1 || type == 7) surfaceArea[iGrid] = 0.75;
        if (type == 3 || type == 5) surfaceArea[iGrid] = 1.25;
      }
      nodeVol[iGrid] = nodeVol[iGrid] * pow(resolution, 3);
      surfaceArea[iGrid] = surfaceArea[iGrid] * pow(resolution, 2);
    }

    // Set initial volumes
    if (volume.empty() && volumeFixed) {
      volume = vector<double>(nFluid, 0);
      for (int iGrid=0; iGrid<nGrid; iGrid++) {
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          volume[iFluid] += coords[iFluid+nFluid*iGrid] * nodeVol[iGrid];
        }
      }
    }

    // Assign elements
    elements = {};
    for (int i=0; i<nGrid; i++) {
      if (solid[i]) continue;

      // Set bulk fluid elements
      Neighbours di(gridSize, i);
      vector<int> idofs = {i, di[0], di[1], di[2], di[3], di[4], di[5]};
      for (auto &idof: idofs) {
        if (solid[idof]) idof = i;
      }
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        elements.push_back({0, idofs*nFluid+iFluid, {nodeVol[i], (double)iFluid}});
      }

      // Set soft density constraint elements
      if (nFluid>1 && densityConstraint==2) {
        idofs = vector<int>(nFluid);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          idofs[iFluid] = i * nFluid + iFluid;
        }
        elements.push_back({1, idofs});
      }

      // Set surface fluid elements
      if (surfaceArea[i] > 0 && !contactAngle.empty() && contactAngle[i]!=90) {
        double wettingParam = 1/sqrt(2.0) * cos(contactAngle[i] * 3.1415926536/180);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          int idof = i * nFluid + iFluid;
          elements.push_back({2, {idof}, {surfaceArea[i], wettingParam}});
        }
      }

      // Set external force elements
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (fMag[iFluid]>0 && !solid[i]) {
          std::array<int,3> coordI = getCoord(i);
          vector<double> coord{coordI[0]-(gridSize[0]-1)/2.0, coordI[1]-(gridSize[1]-1)/2.0, coordI[2]-(gridSize[2]-1)/2.0};
          double h = - vec::dotProduct(coord, fNorm[iFluid]) * resolution;
          vector<double> params{nodeVol[i], fMag[iFluid], h};
          int idof = i * nFluid + iFluid;
          elements.push_back({3, {idof}, params});
        }
      }

      // Set confining potential elements for the frozen fluid method
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (confinementStrength[iFluid] == 0) continue;
        int idof = i * nFluid + iFluid;
        elements.push_back({4, {idof}, {confinementStrength[iFluid], coords[idof]}});
      }
    }

    // Set constraints
    constraints = {};
    // Fixed fluid
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (!fixFluid[iFluid]) continue;
      vector<int> idofs = iFluid + nFluid*vec::iota(nGrid);
      setConstraints(idofs);
    }
    // Hard density constraint
    if (nFluid>1 && densityConstraint==1) {
      vector<int> iVariableFluid;
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (!fixFluid[iFluid]) iVariableFluid.push_back(iFluid);
      }
      vector2d<int> idofs(nGrid);
      for (int iGrid=0; iGrid<nGrid; iGrid++) idofs[iGrid] = iGrid*nFluid + iVariableFluid;
      setConstraints(idofs, vector<double>(iVariableFluid.size(), 1));
    }
  }


  void PFWetting::distributeParameters(const Communicator& comm) {
    if (nGrid % comm.size() != 0) {
      throw std::invalid_argument("The total grid size must be a multiple of the number of processors.");
    }
    // Store only the relevant node volumes and fluid numbers
    // These are required by pressure / volume constraints so cannot be stored in element parameters
    vector<double> nodeVolTmp(nGrid * nFluid);
    vector<double> fluidTypeTmp(nGrid * nFluid);
    for (int iNode=0; iNode<nGrid; iNode++) {
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        nodeVolTmp[nFluid*iNode+iFluid] = nodeVol[iNode];
        fluidTypeTmp[nFluid*iNode+iFluid] = iFluid;
      }
    }
    nodeVol = comm.assignProc(nodeVolTmp);
    fluidTypeTmp = comm.assignProc(fluidTypeTmp);
    fluidType = vector<int>(fluidTypeTmp.begin(), fluidTypeTmp.end());
  }


  void PFWetting::blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {
    // Constant volume / pressure constraints rely upon the whole system
    if (!volumeFixed && !vec::any(pressure)) return;
    vector<double> volFluid(nFluid, 0);
    for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
      if (nFluid == 1) {
        volFluid[0] += 0.5*(coords[iDof]+1) * nodeVol[iDof];
      } else {
        volFluid[fluidType[iDof]] += coords[iDof] * nodeVol[iDof];
      }
    }
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (!volumeFixed && pressure[iFluid]==0) continue;
      volFluid[iFluid] = comm.sum(volFluid[iFluid]);
      if (e) {
        if (volumeFixed) *e += volConst * pow(volFluid[iFluid] - volume[iFluid], 2) / comm.size();
        if (pressure[iFluid]!=0) *e -= pressure[iFluid] * volFluid[iFluid] / comm.size();
      }
    }

    if (!g) return;
    if (volumeFixed) {
      vector<double> vFactor = volConst * (volFluid - volume);
      for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
        int f = fluidType[iDof];
        (*g)[iDof] += vFactor[f] * nodeVol[iDof];
      }
    }
    if (vec::any(pressure)) {
      vector<double> pFactor = (nFluid==1) ? -0.5*pressure : -pressure;
      for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
        int f = fluidType[iDof];
        if (pressure[f]!=0) (*g)[iDof] += pFactor[f] * nodeVol[iDof];
      }
    }
  }


  void PFWetting::elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    switch (el.type) {
      case 0: {
        double vol = el.parameters[0];
        int iFluid = el.parameters[1];
        double c = coords[el.idof[0]];
        // Bulk energy
        if (nFluid == 1) {
          double factor = kappa[0] / 16 * vol;
          if (e) *e += factor * pow(c+1, 2) * pow(c-1, 2);
          if (g) (*g)[el.idof[0]] += factor * 4 * c * (c*c - 1);
        } else {
          double factor = 0.5 * kappa[iFluid] * vol;
          if (e) *e += factor * pow(c, 2) * pow(c-1, 2);
          if (g) (*g)[el.idof[0]] += factor * 2 * c * (c-1) * (2*c-1);
        }
        // Gradient energy
        double factor = (nFluid==1) ? 0.25*kappaP[0]*vol : 0.5*kappaP[iFluid]*vol;
        double res2 = pow(resolution, 2);
        double grad2 = 0;

        if ((el.idof[1]!=el.idof[0]) && (el.idof[6]!=el.idof[0])) { // No solid in x direction
          double diffxm = c - coords[el.idof[1]];
          double diffxp = c - coords[el.idof[6]];
          if (e) grad2 += (pow(diffxm,2) + pow(diffxp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffxm + diffxp) / res2;
            (*g)[el.idof[1]] -= factor * diffxm / res2;
            (*g)[el.idof[6]] -= factor * diffxp / res2;
          }
        } else if (el.idof[6] != el.idof[0]) { // Solid on negative side in x direction
          double diffxp = c - coords[el.idof[6]];
          if (e) grad2 += pow(diffxp,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffxp / res2;
            (*g)[el.idof[6]] -= factor * 2*diffxp / res2;
          }
        } else if (el.idof[1] != el.idof[0]) { // Solid on positive side in x direction
          double diffxm = c - coords[el.idof[1]];
          if (e) grad2 += pow(diffxm,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffxm / res2;
            (*g)[el.idof[1]] -= factor * 2*diffxm / res2;
          }
        }

        if ((el.idof[2]!=el.idof[0]) && (el.idof[5]!=el.idof[0])) { // No solid in y direction
          double diffym = c - coords[el.idof[2]];
          double diffyp = c - coords[el.idof[5]];
          if (e) grad2 += (pow(diffym,2) + pow(diffyp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffym + diffyp) / res2;
            (*g)[el.idof[2]] -= factor * diffym / res2;
            (*g)[el.idof[5]] -= factor * diffyp / res2;
          }
        } else if (el.idof[5] != el.idof[0]) { // Solid on negative side in y direction
          double diffyp = c - coords[el.idof[5]];
          if (e) grad2 += pow(diffyp,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffyp / res2;
            (*g)[el.idof[5]] -= factor * 2*diffyp / res2;
          }
        } else if (el.idof[2] != el.idof[0]) { // Solid on positive side in y direction
          double diffym = c - coords[el.idof[2]];
          if (e) grad2 += pow(diffym,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffym / res2;
            (*g)[el.idof[2]] -= factor * 2*diffym / res2;
          }
        }

        if ((el.idof[3]!=el.idof[0]) && (el.idof[4]!=el.idof[0])) { // No solid in z direction
          double diffzm = c - coords[el.idof[3]];
          double diffzp = c - coords[el.idof[4]];
          if (e) grad2 += (pow(diffzm,2) + pow(diffzp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffzm + diffzp) / res2;
            (*g)[el.idof[3]] -= factor * diffzm / res2;
            (*g)[el.idof[4]] -= factor * diffzp / res2;
          }
        } else if (el.idof[4] != el.idof[0]) { // Solid on negative side in z direction
          double diffzp = c - coords[el.idof[4]];
          if (e) grad2 += pow(diffzp,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffzp / res2;
            (*g)[el.idof[4]] -= factor * 2*diffzp / res2;
          }
        } else if (el.idof[3] != el.idof[0]) { // Solid on positive side in z direction
          double diffzm = c - coords[el.idof[3]];
          if (e) grad2 += pow(diffzm,2) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*diffzm / res2;
            (*g)[el.idof[3]] -= factor * 2*diffzm / res2;
          }
        }
        if (e) *e += factor * grad2;
      } break;

      case 1: {
        // Total density soft constraint
        if (nFluid == 1) break;
        double coef = volConst * pow(resolution, 3);
        double rhoDiff = coords[el.idof[0]] + coords[el.idof[1]] + coords[el.idof[2]] - 1;
        if (e) *e += coef * pow(rhoDiff, 2);
        if (g) {
          (*g)[el.idof[0]] += 2 * coef * rhoDiff;
          (*g)[el.idof[1]] += 2 * coef * rhoDiff;
          (*g)[el.idof[2]] += 2 * coef * rhoDiff;
        }
      } break;

      case 2: {
        // Surface energy
        // parameter[0]: Surface area
        // parameter[1]: Wetting parameter 1/sqrt(2) cos(theta)
        if (nFluid == 1) {
          double phi = coords[el.idof[0]];
          if (e) *e += el.parameters[1] * (pow(phi,3)/3 - phi - 2.0/3) * el.parameters[0];
          if (g) (*g)[el.idof[0]] += el.parameters[1] * (pow(phi,2) - 1) * el.parameters[0];
        }
      } break;

      case 3: {
        // External force
        // Parameters:
        //   0: Volume
        //   1: Magnitude of force
        //   2: Height
        double c = (nFluid==1) ? 0.5*(1+coords[el.idof[0]]) : coords[el.idof[0]];
        if (c < 0.01) return;
        double vol = el.parameters[0];
        double f = el.parameters[1];
        double h = el.parameters[2];
        if (e) *e += c * vol * f * h;
        if (g) {
          if (nFluid == 1) {
            (*g)[el.idof[0]] += 0.5 * vol * f * h;
          } else {
            (*g)[el.idof[0]] += vol * f * h;
          }
        }
      } break;

      case 4: {
        // Confining potential for the frozen fluid method
        // Parameters:
        //   0: Confinement strength
        //   1: Initial concentration
        double c = coords[el.idof[0]];
        double strength = el.parameters[0];
        double c0 = el.parameters[1];
        if ((c0-0.5)*(c-0.5) < 0) {
          if (e) *e += strength * pow(c-0.5, 2);
          if (g) (*g)[el.idof[0]] += strength * 2 * (c-0.5);
        }
      } break;

      default:
        throw std::invalid_argument("Unknown energy element type.");
    }
  }


  PFWetting& PFWetting::setGridSize(std::array<int,3> gridSize) {
    this->gridSize = gridSize;
    this->nGrid = gridSize[0] * gridSize[1] * gridSize[2];
    return *this;
  }

  PFWetting& PFWetting::setNFluid(int nFluid) {
    this->nFluid = nFluid;
    return *this;
  }

  PFWetting& PFWetting::setResolution(double resolution) {
    this->resolution = resolution;
    return *this;
  }

  PFWetting& PFWetting::setInterfaceSize(double interfaceSize) {
    this->interfaceSize = vector<double>(nFluid, interfaceSize);
    return *this;
  }

  PFWetting& PFWetting::setInterfaceSize(vector<double> interfaceSize) {
    this->interfaceSize = interfaceSize;
    return *this;
  }

  PFWetting& PFWetting::setSurfaceTension(double surfaceTension) {
    this->surfaceTension = vector<double>(nFluid, surfaceTension);
    return *this;
  }

  PFWetting& PFWetting::setSurfaceTension(vector<double> surfaceTension) {
    this->surfaceTension = surfaceTension;
    return *this;
  }

  PFWetting& PFWetting::setDensityConstraint(std::string method) {
    if (method == "none") {
      densityConstraint = 0;
    } else if (vec::isIn({"gradient","hard"}, method)) {
      densityConstraint = 1;
    } else if (vec::isIn({"energy penalty","soft"}, method)) {
      densityConstraint = 2;
    } else {
      throw std::invalid_argument("Invalid density constraint. Allowed methods are: none, gradient / hard, or energy penalty / soft");
    }
    return *this;
  }

  PFWetting& PFWetting::setPressure(vector<double> pressure) {
    this->pressure = pressure;
    return *this;
  }

  PFWetting& PFWetting::setVolume(vector<double> volume, double volConst) {
    this->volumeFixed = true;
    this->volConst = volConst;
    this->volume = volume;
    return *this;
  }

  PFWetting& PFWetting::setVolumeFixed(bool volumeFixed, double volConst) {
    this->volumeFixed = volumeFixed;
    this->volConst = volConst;
    return *this;
  }

  PFWetting& PFWetting::setSolid(vector<bool> solid) {
    this->solid = solid;
    return *this;
  }

  PFWetting& PFWetting::setSolid(std::function<bool(int,int,int)> solidFn) {
    solid = vector<bool>(nGrid);
    int itot = 0;
    for (int i=0; i<gridSize[0]; i++) {
      for (int j=0; j<gridSize[1]; j++) {
        for (int k=0; k<gridSize[2]; k++) {
          solid[itot] = solidFn(i, j, k);
          itot ++;
        }
      }
    }
    return *this;
  }

  PFWetting& PFWetting::setContactAngle(vector<double> contactAngle) {
    this->contactAngle = contactAngle;
    return *this;
  }

  PFWetting& PFWetting::setContactAngle(std::function<double(int,int,int)> contactAngleFn) {
    contactAngle = vector<double>(nGrid);
    int itot = 0;
    for (int i=0; i<gridSize[0]; i++) {
      for (int j=0; j<gridSize[1]; j++) {
        for (int k=0; k<gridSize[2]; k++) {
          contactAngle[itot] = contactAngleFn(i, j, k);
          itot ++;
        }
      }
    }
    return *this;
  }

  PFWetting& PFWetting::setForce(vector<double> force, vector<int> iFluid) {
    if ((int)force.size() != 3) throw std::invalid_argument("Invalid size of force array.");
    if (nFluid==1 || iFluid.empty()) {
      this->force = vector2d<double>(nFluid, force);
    } else {
      this->force = vector2d<double>(nFluid, {0,0,0});
      for (int i: iFluid) {
        this->force[i] = force;
      }
    }
    return *this;
  }

  PFWetting& PFWetting::setFixFluid(int iFluid, bool fix) {
    if ((int)fixFluid.size()!=nFluid) fixFluid = vector<bool>(nFluid, false);
    fixFluid[iFluid] = fix;
    return *this;
  }

  PFWetting& PFWetting::setConfinement(vector<double> strength) {
    this->confinementStrength = strength;
    return *this;
  }


  std::array<int,3> PFWetting::getCoord(int i) const {
    return minim::getCoord(i, gridSize);
  }


  int PFWetting::getType(int i) const {
    // Types:
    // -1: Solid (▪▪)   0: Bulk fluid (  )
    //           (▪▪)                 (  )
    // 1: ▫    2: ▫▫   3: ▫▫   4: ▫▫   5: ▪▫   6: ▪▪   7: ▪▪
    //                    ▫       ▫▫      ▫▫      ▫▫      ▪▫
    if (solid[i]) return -1;

    // Get which neighbouring nodes are solid
    Neighbours di(gridSize, i);
    std::array<bool,26> neiSolid;
    for (int iNei=0; iNei<26; iNei++) {
      neiSolid[iNei] = solid[di[iNei]];
    }
    int nSolidF = neiSolid[0] + neiSolid[1] + neiSolid[2] + neiSolid[3] + neiSolid[4] + neiSolid[5];
    int nSolidE = neiSolid[6]  + neiSolid[7]  + neiSolid[8]  + neiSolid[9]  +
                  neiSolid[10] + neiSolid[11] + neiSolid[12] + neiSolid[13] +
                  neiSolid[14] + neiSolid[15] + neiSolid[16] + neiSolid[17];
    int nSolidC = neiSolid[18] + neiSolid[19] + neiSolid[20] + neiSolid[21] +
                  neiSolid[22] + neiSolid[23] + neiSolid[24] + neiSolid[25];
    int nSolid = nSolidF + nSolidE + nSolidC;

    // Bulk node
    if (nSolid == 0) return 0;

    // Determine surface type
    if (nSolidF == 0) {
      if (nSolidE == 0 && nSolidC == 1) {
        return 1;
      } else if (nSolidE == 1 && nSolidC <= 2) {
        return 2;
      } else if (nSolidE == 2 && nSolidC <= 3) {
        return 3;
      }

    } else if (nSolidF == 1) {
      bool oneSide =
        ((neiSolid[0] || neiSolid[5]) && !(neiSolid[10] || neiSolid[11] || neiSolid[12] || neiSolid[13])) ||
        ((neiSolid[1] || neiSolid[4]) && !(neiSolid[7]  || neiSolid[8]  || neiSolid[15] || neiSolid[16])) ||
        ((neiSolid[2] || neiSolid[3]) && !(neiSolid[6]  || neiSolid[9]  || neiSolid[14] || neiSolid[17]));
      if (oneSide) {
        return 4;
      } else {
        return 5;
      }

    } else if (nSolidF == 2) {
      return 6;

    } else if (nSolidF == 3) {
      return 7;
    }
    throw std::runtime_error("Undefined surface type");
  }


  void PFWetting::setDefaults() {
    if (interfaceSize.empty()) interfaceSize = vector<double>(nFluid, resolution);
    if (solid.empty()) solid = vector<bool>(nGrid, false);
    if (pressure.empty()) pressure = vector<double>(nFluid, 0);
    if (fixFluid.empty()) fixFluid = vector<bool>(nFluid, false);
    if (confinementStrength.empty()) confinementStrength = vector<double>(nFluid, 0);
  }


  void PFWetting::checkArraySizes() {
    if ((int)interfaceSize.size() != nFluid) {
      throw std::invalid_argument("Invalid size interfaceSize array.");
    }
    if ((int)solid.size() != nGrid) {
      throw std::invalid_argument("Invalid size of solid array.");
    }
    if ((int)contactAngle.size() != nGrid*nFluid && !contactAngle.empty()) {
      throw std::invalid_argument("Invalid size of contactAngle array.");
    }
    if ((int)force.size() != nFluid && !force.empty()) {
      throw std::invalid_argument("Invalid size of force array.");
    }
    if ((int)volume.size() != nFluid && !volume.empty()) {
      throw std::invalid_argument("Invalid size of volume array.");
    }
    if ((int)pressure.size() != nFluid) {
      throw std::invalid_argument("Invalid size of pressure array.");
    }
    if ((int)fixFluid.size() != nFluid) {
      throw std::invalid_argument("Invalid size of fixed fluid array.");
    }
    if ((int)confinementStrength.size() != nFluid) {
      throw std::invalid_argument("Invalid size of confinement strength array.");
    }
  }

}
