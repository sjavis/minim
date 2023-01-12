#include "PFWetting.h"

#include <array>
#include <math.h>
#include <stdexcept>
#include <functional>
#include "State.h"
#include "utils/vec.h"


namespace minim {
  typedef std::vector<double> Vector;


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
      std::array<std::array<int,3>,26> dx {{
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
      interfaceSize = Vector(nFluid, interfaceSize[0]);
    } else if ((int)interfaceSize.size() != nFluid) {
      throw std::invalid_argument("Invalid size of interfaceSize array.");
    }
    if ((int)surfaceTension.size()==1 && nFluid>1) {
      surfaceTension = Vector(nFluid, surfaceTension[0]);
    } else if ((int)surfaceTension.size() != nFluid) {
      throw std::invalid_argument("Invalid size of surfaceTension array.");
    }
    // TODO: Calculate different values of kappa and kappaP
    kappa = Vector(nFluid, 3*surfaceTension[0]/interfaceSize[0]);
    kappaP = Vector(nFluid, pow(interfaceSize[0], 2)*kappa[0]);
  }


  void PFWetting::init() {
    int nGrid = gridSize[0] * gridSize[1] * gridSize[2];

    // Check the arrays
    if (solid.empty()) {
      solid = std::vector<bool>(nGrid, false);
    } else if ((int)solid.size() != nGrid) {
      throw std::invalid_argument("Invalid size of solid array.");
    }
    if (!contactAngle.empty() && (int)contactAngle.size() != nGrid*nFluid) {
      throw std::invalid_argument("Invalid size of contactAngle array.");
    }
    if (!force.empty() && (int)force.size()!=nFluid) {
      throw std::invalid_argument("Invalid size of force array.");
    }
    if (volume.empty()) {
      volume = Vector(nFluid, 0);
    } else if ((int)volume.size() != nFluid) {
      throw std::invalid_argument("Invalid size of volume array.");
    }
    if (pressure.empty()) {
      pressure = Vector(nFluid, 0);
    } else if ((int)pressure.size() != nFluid) {
      throw std::invalid_argument("Invalid size of pressure array.");
    }

    // Set values
    assignKappa();
    nodeVol = Vector(nGrid);
    Vector fMag(nFluid);
    std::vector<Vector> fNorm(nFluid);
    if (!force.empty()) {
      for (int i=0; i<nFluid; i++) {
        fMag[i] = vec::norm(force[i]);
        fNorm[i] = force[i] / fMag[i];
      }
    }

    // Assign elements
    elements = {};
    for (int i=0; i<nGrid; i++) {
      if (solid[i]) continue;

      // Get fluid volume and solid surface area for each node
      nodeVol[i] = 1;
      double surfaceArea = 0;
      int type = getType(i);
      if (type > 0) {
        nodeVol[i] = type / 8.0;
        if (type == 2 || type == 4 || type == 6) surfaceArea = 1;
        if (type == 1 || type == 7) surfaceArea = 0.75;
        if (type == 3 || type == 5) surfaceArea = 1.25;
      }
      nodeVol[i] = nodeVol[i] * pow(resolution, 3);
      surfaceArea = surfaceArea * pow(resolution, 2);

      // Set bulk fluid elements
      Neighbours di(gridSize, i);
      std::vector<int> idofs = {i, di[0], di[1], di[2], di[3], di[4], di[5]};
      for (auto &idof: idofs) {
        if (solid[idof]) idof = i;
      }
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        elements.push_back({0, idofs*nFluid+iFluid, {nodeVol[i], (double)iFluid}});
      }

      // Set density constraint elements
      if (nFluid > 1) {
        idofs = std::vector<int>(nFluid);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          idofs[iFluid] = i * nFluid + iFluid;
        }
        elements.push_back({1, idofs});
      }

      // Set surface fluid elements
      if (type > 0 && !contactAngle.empty() && contactAngle[i]!=90) {
        double wettingParam = sqrt(2.0) * cos(contactAngle[i] * 3.1415926536/180);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          int idof = i * nFluid + iFluid;
          elements.push_back({2, {idof}, {surfaceArea, wettingParam}});
        }
      }

      // Set external force elements
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (fMag[iFluid]>0 && !solid[i]) {
          Vector params{nodeVol[i], fMag[iFluid], fNorm[iFluid][0], fNorm[iFluid][1], fNorm[iFluid][2]};
          for (int iFluid=0; iFluid<nFluid; iFluid++) {
            int idof = i * nFluid + iFluid;
            elements.push_back({3, {idof}, params});
          }
        }
      }
    }
  }


  void PFWetting::distributeParameters(const Communicator& comm) {
    // Store only the relevant node volumes and fluid numbers
    // These are required by pressure / volume constraints so cannot be stored in element parameters
    int nGrid = gridSize[0] * gridSize[1] * gridSize[2];
    Vector nodeVolTmp(nGrid * nFluid);
    Vector fluidTypeTmp(nGrid * nFluid);
    for (int iNode=0; iNode<nGrid; iNode++) {
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        nodeVolTmp[nFluid*iNode+iFluid] = nodeVol[iNode];
        fluidTypeTmp[nFluid*iNode+iFluid] = iFluid;
      }
    }
    nodeVol = comm.assignProc(nodeVolTmp);
    fluidTypeTmp = comm.assignProc(fluidTypeTmp);
    fluidType = std::vector<int>(fluidTypeTmp.begin(), fluidTypeTmp.end());
  }


  void PFWetting::blockEnergyGradient(const Vector& coords, const Communicator& comm, double* e, Vector* g) const {
    // Constant volume / pressure constraints rely upon the whole system
    if (!vec::any(volume) && !vec::any(pressure)) return;
    Vector volFluid(nFluid, 0);
    for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
      if (nFluid == 1) {
        volFluid[0] += 0.5*(coords[iDof]+1) * nodeVol[iDof];
      } else {
        volFluid[fluidType[iDof]] += coords[iDof] * nodeVol[iDof];
      }
    }
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (volume[iFluid]==0 && pressure[iFluid]==0) continue;
      volFluid[iFluid] = comm.sum(volFluid[iFluid]);
      if (e) {
        if (volume[iFluid]!=0) *e += volConst * pow(volFluid[iFluid] - volume[iFluid], 2) / comm.size();
        if (pressure[iFluid]!=0) *e -= pressure[iFluid] * volFluid[iFluid] / comm.size();
      }
    }

    if (!g) return;
    Vector vFactor = volConst * (volFluid - volume);
    Vector pFactor = (nFluid==1) ? -0.5*pressure : -pressure;
    for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
      int f = fluidType[iDof];
      if (volume[f]!=0) (*g)[iDof] += vFactor[f] * nodeVol[iDof];
      if (pressure[f]!=0) (*g)[iDof] += pFactor[f] * nodeVol[iDof];
    }
  }


  void PFWetting::elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const {
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
        double diffxm = (el.idof[1]!=el.idof[0]) ? c - coords[el.idof[1]] : 0;
        double diffym = (el.idof[2]!=el.idof[0]) ? c - coords[el.idof[2]] : 0;
        double diffzm = (el.idof[3]!=el.idof[0]) ? c - coords[el.idof[3]] : 0;
        double diffzp = (el.idof[4]!=el.idof[0]) ? c - coords[el.idof[4]] : 0;
        double diffyp = (el.idof[5]!=el.idof[0]) ? c - coords[el.idof[5]] : 0;
        double diffxp = (el.idof[6]!=el.idof[0]) ? c - coords[el.idof[6]] : 0;
        double grad2 = 0;
        double res2 = pow(resolution, 2);
        if (diffxm != 0 && diffxp != 0) { // No solid in x direction
          if (e) grad2 += (pow(diffxm,2) + pow(diffxp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffxm + diffxp) / res2;
            (*g)[el.idof[1]] -= factor * diffxm / res2;
            (*g)[el.idof[6]] -= factor * diffxp / res2;
          }
        } else { // Solid on one side in x direction
          if (e) grad2 += (pow(diffxm,2) + pow(diffxp,2)) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*(diffxm + diffxp) / res2;
            (*g)[el.idof[1]] -= factor * 2*diffxm / res2;
            (*g)[el.idof[6]] -= factor * 2*diffxp / res2;
          }
        }
        if (diffym != 0 && diffyp != 0) { // No solid in y direction
          if (e) grad2 += (pow(diffym,2) + pow(diffyp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffym + diffyp) / res2;
            (*g)[el.idof[2]] -= factor * diffym / res2;
            (*g)[el.idof[5]] -= factor * diffyp / res2;
          }
        } else { // Solid on one side in y direction
          if (e) grad2 += (pow(diffym,2) + pow(diffyp,2)) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*(diffym + diffyp) / res2;
            (*g)[el.idof[2]] -= factor * 2*diffym / res2;
            (*g)[el.idof[5]] -= factor * 2*diffyp / res2;
          }
        }
        if (diffzm != 0 && diffzp != 0) { // No solid in z direction
          if (e) grad2 += (pow(diffzm,2) + pow(diffzp,2)) / (2 * res2);
          if (g) {
            (*g)[el.idof[0]] += factor * (diffzm + diffzp) / res2;
            (*g)[el.idof[3]] -= factor * diffzm / res2;
            (*g)[el.idof[4]] -= factor * diffzp / res2;
          }
        } else { // Solid on one side in z direction
          if (e) grad2 += (pow(diffzm,2) + pow(diffzp,2)) / res2;
          if (g) {
            (*g)[el.idof[0]] += factor * 2*(diffzm + diffzp) / res2;
            (*g)[el.idof[3]] -= factor * 2*diffzm / res2;
            (*g)[el.idof[4]] -= factor * 2*diffzp / res2;
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
        // parameter[1]: Wetting parameter sqrt(2)cos(theta)
        if (nFluid == 1) {
          double phi = coords[el.idof[0]];
          if (e) *e += el.parameters[1] * (pow(phi,3)/3 - pow(phi,2)/2) * el.parameters[0];
          if (g) (*g)[el.idof[0]] += el.parameters[1] * (pow(phi,2) - phi) * el.parameters[0];
        }
      } break;

      case 3: {
        // External force
        // Parameters:
        //   0: Volume
        //   1: Magnitude of force on component 1
        //   2-4: Direction force on component 1
        double vol = el.parameters[0];
        double f = el.parameters[1];
        Vector fNorm = {el.parameters[2], el.parameters[3], el.parameters[4]};
        std::array<int,3> coordI = getCoord(el.idof[0]);
        Vector coord{coordI[0]-(gridSize[0]-1)/2.0, coordI[1]-(gridSize[1]-1)/2.0, coordI[2]-(gridSize[2]-1)/2.0};
        double h = - vec::dotProduct(coord, fNorm) * resolution;
        if (nFluid == 1) {
          if (e) *e += 0.5*(1+coords[el.idof[0]]) * vol * f * h;
          if (g) (*g)[el.idof[0]] += 0.5 * vol * f * h;
        } else {
          if (e) *e += coords[el.idof[0]] * vol * f * h;
          if (g) (*g)[el.idof[0]] += vol * f * h;
        }
      } break;

      default:
        throw std::invalid_argument("Unknown energy element type.");
    }
  }


  PFWetting& PFWetting::setGridSize(std::array<int,3> gridSize) {
    this->gridSize = gridSize;
    return *this;
  }

  PFWetting& PFWetting::setNFluid(int nFluid) {
    this->nFluid = nFluid;
    return *this;
  }

  PFWetting& PFWetting::setInterfaceSize(double interfaceSize) {
    this->interfaceSize = Vector(nFluid, interfaceSize);
    return *this;
  }

  PFWetting& PFWetting::setSurfaceTension(double surfaceTension) {
    this->surfaceTension = Vector(nFluid, surfaceTension);
    return *this;
  }

  PFWetting& PFWetting::setResolution(double resolution) {
    this->resolution = resolution;
    return *this;
  }

  PFWetting& PFWetting::setPressure(Vector pressure) {
    this->pressure = pressure;
    return *this;
  }

  PFWetting& PFWetting::setVolume(Vector volume, double volConst) {
    this->volume = volume;
    this->volConst = volConst;
    return *this;
  }

  PFWetting& PFWetting::setSolid(std::vector<bool> solid) {
    this->solid = solid;
    return *this;
  }

  PFWetting& PFWetting::setSolid(std::function<bool(int,int,int)> solidFn) {
    solid = std::vector<bool>(gridSize[0]*gridSize[1]*gridSize[2]);
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

  PFWetting& PFWetting::setContactAngle(Vector contactAngle) {
    this->contactAngle = contactAngle;
    return *this;
  }

  PFWetting& PFWetting::setContactAngle(std::function<double(int,int,int)> contactAngleFn) {
    contactAngle = Vector(gridSize[0]*gridSize[1]*gridSize[2]);
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

  PFWetting& PFWetting::setForce(Vector force, std::vector<int> iFluid) {
    if ((int)force.size() != 3) throw std::invalid_argument("Invalid size of force array.");
    if (nFluid==1 || iFluid.empty()) {
      this->force = std::vector<Vector>(nFluid, force);
    } else {
      this->force = std::vector<Vector>(nFluid, {0,0,0});
      for (int i: iFluid) {
        this->force[i] = force;
      }
    }
    return *this;
  }


  State PFWetting::newState(const Vector& coords, const std::vector<int>& ranks) {
    return State(*this, coords, ranks);
  }


  int PFWetting::nGrid() const {
    return gridSize[0] * gridSize[1] * gridSize[2];
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

}
