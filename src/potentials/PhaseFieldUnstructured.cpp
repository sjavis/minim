#include "potentials/PhaseFieldUnstructured.h"

#include <math.h>
#include <stdexcept>
#include <functional>
#include "State.h"
#include "utils/vec.h"
#include "minimisers/Lbfgs.h"


namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;

  enum{
    DENSITY_HARD = 0, // Apply a constraint to the gradient
    DENSITY_SOFT = 1, // Use an energy penalty
    DENSITY_NONE = 2, // No density constraint, independent fluids
    DENSITY_FIXED = 3, // Use N-1 fluids, calculate Nth fluid concentration
  };

  enum{
    FLUID_ENERGY = 0,
    FLUID_ENERGY_ALL = 1,
    FLUID_PAIR_ENERGY = 2,
    PRESSURE_ENERGY = 3,
    DENSITY_CONSTRAINT_ENERGY = 4,
    SURFACE_ENERGY = 5,
    FORCE_ENERGY = 6,
    FF_CONFINEMENT_ENERGY = 7,
  };


  vector<int> getCoord(int i, vector<int> gridSize) {
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

      Neighbours(vector<int> gridSize, int i0) {
        auto x0 = getCoord(i0, gridSize);
        for (int i=0; i<26; i++) {
          int x = (x0[0] + dx[i][0] + gridSize[0]) % gridSize[0];
          int y = (x0[1] + dx[i][1] + gridSize[1]) % gridSize[1];
          int z = (x0[2] + dx[i][2] + gridSize[2]) % gridSize[2];
          di[i] = (x*gridSize[1] + y)*gridSize[2] + z;
        }
      }
  };


  void PhaseFieldUnstructured::assignFluidCoefficients() {
    int nParams = (model==MODEL_BASIC) ? nFluid : 0.5*nFluid*(nFluid-1);

    // Ensure surface tension and interface widths are the correct size
    if ((int)interfaceSize.size() == 1) {
      interfaceSize = vector<double>(nParams, interfaceSize[0]);
    } else if ((int)interfaceSize.size() != nParams) {
      throw std::invalid_argument("Invalid size of interfaceSize array.");
    }
    if ((int)surfaceTension.size() == 1) {
      surfaceTension = vector<double>(nParams, surfaceTension[0]);
    } else if ((int)surfaceTension.size() != nParams) {
      throw std::invalid_argument("Invalid size of surfaceTension array.");
    }
    surfaceTensionMean = vec::sum(surfaceTension) / nParams; // Used to scale energies

    // Set the parameters
    if (nFluid<=2 || model==MODEL_NCOMP) {
      kappa = 3 * surfaceTension / interfaceSize;
      kappaP = 3 * surfaceTension * interfaceSize;
    } else { // Compute kappa from the subset of surface tensions
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


  void PhaseFieldUnstructured::init(const vector<double>& coords) {
    if ((int)coords.size() != nGrid*nFluid) {
      throw std::invalid_argument("Size of coordinates array does not match the grid size.");
    }

    setDefaults();
    checkArraySizes();
    assignFluidCoefficients();
    this->convergence = 1e-8 * surfaceTensionMean * pow(resolution, 2);
    if (densityConstraint == DENSITY_FIXED) fixFluid[nFluid-1] = true;

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
    for (int iGrid=0; iGrid<nGrid; iGrid++) {
      if (solid[iGrid]) continue;

      // Set bulk fluid elements
      Neighbours di(gridSize, iGrid);
      vector<int> iNodes = {iGrid, di[0], di[1], di[2], di[3], di[4], di[5]};
      for (auto &iNode: iNodes) {
        if (solid[iNode]) iNode = iGrid;
      }
      if (model == MODEL_BASIC) {
        if (densityConstraint == DENSITY_FIXED) {
          elements.push_back({FLUID_ENERGY_ALL, iNodes*nFluid, {nodeVol[iGrid]}});
        } else {
          for (int iFluid=0; iFluid<nFluid; iFluid++) {
            vector<int> idofs = iNodes*nFluid+iFluid;
            elements.push_back({FLUID_ENERGY, idofs, {nodeVol[iGrid], (double)iFluid}});
          }
        }
      } else if (model == MODEL_NCOMP) {
        int iPair = 0;
        for (int iFluid1=0; iFluid1<nFluid; iFluid1++) {
          for (int iFluid2=iFluid1+1; iFluid2<nFluid; iFluid2++) {
            vector<int> idofs;
            idofs.reserve(14);
            for (auto iNode: iNodes) {
              idofs.push_back(iNode*nFluid+iFluid1);
              idofs.push_back(iNode*nFluid+iFluid2);
            }
            elements.push_back({FLUID_PAIR_ENERGY, idofs, {nodeVol[iGrid], (double)iPair}});
            iPair++;
          }
        }
      }

      // Pressure elements
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (pressure[iFluid]==0) continue;
        int iDof = iGrid*nFluid + iFluid;
        elements.push_back({PRESSURE_ENERGY, {iDof}});
      }

      // Set soft density constraint elements
      if (nFluid>1 && densityConstraint==DENSITY_SOFT) {
        vector<int> idofs = vector<int>(nFluid);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          idofs[iFluid] = iGrid * nFluid + iFluid;
        }
        elements.push_back({DENSITY_CONSTRAINT_ENERGY, idofs});
      }

      // Set surface fluid elements
      if (surfaceArea[iGrid] > 0 && !contactAngle.empty() && contactAngle[iGrid]!=90) {
        double wettingParam = 1/sqrt(2.0) * cos(contactAngle[iGrid] * 3.1415926536/180);
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          int idof = iGrid * nFluid + iFluid;
          elements.push_back({SURFACE_ENERGY, {idof}, {surfaceArea[iGrid], wettingParam}});
        }
      }

      // Set external force elements
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (fMag[iFluid]>0 && !solid[iGrid]) {
          vector<int> coordI = getCoord(iGrid);
          vector<double> coord{coordI[0]-(gridSize[0]-1)/2.0, coordI[1]-(gridSize[1]-1)/2.0, coordI[2]-(gridSize[2]-1)/2.0};
          double h = - vec::dotProduct(coord, fNorm[iFluid]) * resolution;
          vector<double> params{nodeVol[iGrid], fMag[iFluid], h};
          int idof = iGrid * nFluid + iFluid;
          elements.push_back({FORCE_ENERGY, {idof}, params});
        }
      }

      // Set confining potential elements for the frozen fluid method
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (confinementStrength[iFluid] == 0) continue;
        int idof = iGrid * nFluid + iFluid;
        elements.push_back({FF_CONFINEMENT_ENERGY, {idof}, {confinementStrength[iFluid], coords[idof]}});
      }
    }

    // Set constraints
    constraints = {};
    // Hard density constraint
    if (nFluid>1 && densityConstraint==DENSITY_HARD) {
      vector<int> iVariableFluid;
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        if (!fixFluid[iFluid]) iVariableFluid.push_back(iFluid);
      }
      vector2d<int> idofs(nGrid);
      for (int iGrid=0; iGrid<nGrid; iGrid++) idofs[iGrid] = iGrid*nFluid + iVariableFluid;
      setConstraints(idofs, vector<double>(iVariableFluid.size(), 1));
    }
    // Fixed fluid
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (!fixFluid[iFluid]) continue;
      vector<int> idofs = iFluid + nFluid*vec::iota(nGrid);
      setConstraints(idofs);
    }
  }


  void PhaseFieldUnstructured::distributeParameters(const Communicator& comm) {
    if (nGrid % comm.size() != 0) {
      throw std::invalid_argument("The total grid size must be a multiple of the number of processors.");
    }
    // Store only the relevant node volumes and fluid numbers
    // These are required by volume constraint so cannot be stored in element parameters
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


  void PhaseFieldUnstructured::blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {
    // Constant volume constraint relies upon the whole system
    if (!volumeFixed) return;

    // Get the fluid volumes
    vector<double> volFluid(nFluid, 0);
    for (int iDof=0; iDof<(int)comm.nblock; iDof++) {
      if (nFluid == 1) {
        volFluid[0] += 0.5*(coords[iDof]+1) * nodeVol[iDof];
      } else {
        volFluid[fluidType[iDof]] += coords[iDof] * nodeVol[iDof];
      }
    }
    for (double &vol : volFluid) vol = comm.sum(vol);

    // Compute the energy and gradient
    double volCoef = volConst * surfaceTensionMean / pow(resolution, 4);
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      double volDiff = volFluid[iFluid] - volume[iFluid];
      if (e) *e += volCoef * pow(volDiff, 2) / comm.size();
      if (!g) continue;
      int nGrid = comm.nblock / nFluid;
      for (int iGrid=0; iGrid<nGrid; iGrid++) {
        int iDof = iGrid*nFluid + iFluid;
        double interfaceWeight = std::max(0.0, 4*coords[iDof]*(1-coords[iDof])); // Only apply the force to the interface nodes
        (*g)[iDof] += volCoef * volDiff * nodeVol[iDof] * interfaceWeight;
      }
    }
  }


  void phaseGradient(const vector<double>& coords, const vector<int> idof, double factor, double* e, vector<double>* g) {
    double c1 = coords[idof[0]];
    double grad2 = 0;

    // Single phase
    if (idof.size() == 7) {

      // X-gradient
      if ((idof[1]!=idof[0]) && (idof[6]!=idof[0])) { // No solid in x direction
        double gradm = c1 - coords[idof[1]];
        double gradp = c1 - coords[idof[6]];
        if (e) grad2 += 0.5 * (pow(gradm,2) + pow(gradp,2));
        if (g) {
          (*g)[idof[0]] += factor * (gradm + gradp);
          (*g)[idof[1]] -= factor * gradm;
          (*g)[idof[6]] -= factor * gradp;
        }
      } else if (idof[6] != idof[0]) { // Solid on negative side in x direction
        double gradp = c1 - coords[idof[6]];
        if (e) grad2 += pow(gradp,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradp;
          (*g)[idof[6]] -= factor * 2*gradp;
        }
      } else if (idof[1] != idof[0]) { // Solid on positive side in x direction
        double gradm = c1 - coords[idof[1]];
        if (e) grad2 += pow(gradm,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm;
          (*g)[idof[1]] -= factor * 2*gradm;
        }
      }

      // Y-gradient
      if ((idof[2]!=idof[0]) && (idof[5]!=idof[0])) { // No solid in y direction
        double gradm = c1 - coords[idof[2]];
        double gradp = c1 - coords[idof[5]];
        if (e) grad2 += 0.5 * (pow(gradm,2) + pow(gradp,2));
        if (g) {
          (*g)[idof[0]] += factor * (gradm + gradp);
          (*g)[idof[2]] -= factor * gradm;
          (*g)[idof[5]] -= factor * gradp;
        }
      } else if (idof[5] != idof[0]) { // Solid on negative side in y direction
        double gradp = c1 - coords[idof[5]];
        if (e) grad2 += pow(gradp,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradp;
          (*g)[idof[5]] -= factor * 2*gradp;
        }
      } else if (idof[2] != idof[0]) { // Solid on positive side in y direction
        double gradm = c1 - coords[idof[2]];
        if (e) grad2 += pow(gradm,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm;
          (*g)[idof[2]] -= factor * 2*gradm;
        }
      }

      // Z-gradient
      if ((idof[3]!=idof[0]) && (idof[4]!=idof[0])) { // No solid in z direction
        double gradm = c1 - coords[idof[3]];
        double gradp = c1 - coords[idof[4]];
        if (e) grad2 += 0.5 * (pow(gradm,2) + pow(gradp,2));
        if (g) {
          (*g)[idof[0]] += factor * (gradm + gradp);
          (*g)[idof[3]] -= factor * gradm;
          (*g)[idof[4]] -= factor * gradp;
        }
      } else if (idof[4] != idof[0]) { // Solid on negative side in z direction
        double gradp = c1 - coords[idof[4]];
        if (e) grad2 += pow(gradp,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradp;
          (*g)[idof[4]] -= factor * 2*gradp;
        }
      } else if (idof[3] != idof[0]) { // Solid on positive side in z direction
        double gradm = c1 - coords[idof[3]];
        if (e) grad2 += pow(gradm,2);
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm;
          (*g)[idof[3]] -= factor * 2*gradm;
        }
      }

      if (e) *e += factor * grad2;

    // Pair of phases
    } else if (idof.size() == 14) {
      double c2 = coords[idof[1]];

      // X-gradient
      if ((idof[2]!=idof[0]) && (idof[12]!=idof[0])) { // No solid in x direction
        double gradm1 = c1 - coords[idof[2]];
        double gradm2 = c2 - coords[idof[3]];
        double gradp1 = c1 - coords[idof[12]];
        double gradp2 = c2 - coords[idof[13]];
        if (e) grad2 += 0.5 * (gradm1*gradm2 + gradp1*gradp2);
        if (g) {
          (*g)[idof[0]] += factor * (gradm2 + gradp2);
          (*g)[idof[1]] += factor * (gradm1 + gradp1);
          (*g)[idof[2]] -= factor * gradm2;
          (*g)[idof[3]] -= factor * gradm1;
          (*g)[idof[12]] -= factor * gradp2;
          (*g)[idof[13]] -= factor * gradp1;
        }
      } else if (idof[12] != idof[0]) { // Solid on negative side in x direction
        double gradp1 = c1 - coords[idof[12]];
        double gradp2 = c2 - coords[idof[13]];
        if (e) grad2 += gradp1*gradp2;
        if (g) {
          (*g)[idof[0]] += factor * 2*gradp2;
          (*g)[idof[1]] += factor * 2*gradp1;
          (*g)[idof[12]] -= factor * 2*gradp2;
          (*g)[idof[13]] -= factor * 2*gradp1;
        }
      } else if (idof[2] != idof[0]) { // Solid on positive side in x direction
        double gradm1 = c1 - coords[idof[2]];
        double gradm2 = c2 - coords[idof[3]];
        if (e) grad2 += gradm1*gradm2;
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm2;
          (*g)[idof[1]] += factor * 2*gradm1;
          (*g)[idof[2]] -= factor * 2*gradm2;
          (*g)[idof[3]] -= factor * 2*gradm1;
        }
      }

      // Y-gradient
      if ((idof[4]!=idof[0]) && (idof[10]!=idof[0])) { // No solid in y direction
        double gradm1 = c1 - coords[idof[4]];
        double gradm2 = c2 - coords[idof[5]];
        double gradp1 = c1 - coords[idof[10]];
        double gradp2 = c2 - coords[idof[11]];
        if (e) grad2 += 0.5 * (gradm1*gradm2 + gradp1*gradp2);
        if (g) {
          (*g)[idof[0]] += factor * (gradm2 + gradp2);
          (*g)[idof[1]] += factor * (gradm1 + gradp1);
          (*g)[idof[4]] -= factor * gradm2;
          (*g)[idof[5]] -= factor * gradm1;
          (*g)[idof[10]] -= factor * gradp2;
          (*g)[idof[11]] -= factor * gradp1;
        }
      } else if (idof[10] != idof[0]) { // Solid on negative side in y direction
        double gradp1 = c1 - coords[idof[10]];
        double gradp2 = c2 - coords[idof[11]];
        if (e) grad2 += gradp1*gradp2;
        if (g) {
          (*g)[idof[0]] += factor * gradp2;
          (*g)[idof[1]] += factor * gradp1;
          (*g)[idof[10]] -= factor * gradp2;
          (*g)[idof[11]] -= factor * gradp1;
        }
      } else if (idof[4] != idof[0]) { // Solid on positive side in y direction
        double gradm1 = c1 - coords[idof[4]];
        double gradm2 = c2 - coords[idof[5]];
        if (e) grad2 += gradm1*gradm2;
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm2;
          (*g)[idof[1]] += factor * 2*gradm1;
          (*g)[idof[4]] -= factor * 2*gradm2;
          (*g)[idof[5]] -= factor * 2*gradm1;
        }
      }

      // Z-gradient
      if ((idof[6]!=idof[0]) && (idof[8]!=idof[0])) { // No solid in z direction
        double gradm1 = c1 - coords[idof[6]];
        double gradm2 = c2 - coords[idof[7]];
        double gradp1 = c1 - coords[idof[8]];
        double gradp2 = c2 - coords[idof[9]];
        if (e) grad2 += 0.5 * (gradm1*gradm2 + gradp1*gradp2);
        if (g) {
          (*g)[idof[0]] += factor * (gradm2 + gradp2);
          (*g)[idof[1]] += factor * (gradm1 + gradp1);
          (*g)[idof[6]] -= factor * gradm2;
          (*g)[idof[7]] -= factor * gradm1;
          (*g)[idof[8]] -= factor * gradp2;
          (*g)[idof[9]] -= factor * gradp1;
        }
      } else if (idof[8] != idof[0]) { // Solid on negative side in z direction
        double gradp1 = c1 - coords[idof[8]];
        double gradp2 = c2 - coords[idof[9]];
        if (e) grad2 += gradp1*gradp2;
        if (g) {
          (*g)[idof[0]] += factor * gradp2;
          (*g)[idof[1]] += factor * gradp1;
          (*g)[idof[8]] -= factor * gradp2;
          (*g)[idof[9]] -= factor * gradp1;
        }
      } else if (idof[6] != idof[0]) { // Solid on positive side in z direction
        double gradm1 = c1 - coords[idof[6]];
        double gradm2 = c2 - coords[idof[7]];
        if (e) grad2 += gradm1*gradm2;
        if (g) {
          (*g)[idof[0]] += factor * 2*gradm2;
          (*g)[idof[1]] += factor * 2*gradm1;
          (*g)[idof[6]] -= factor * 2*gradm2;
          (*g)[idof[7]] -= factor * 2*gradm1;
        }
      }

      if (e) *e += factor * grad2;
    }
  }


  void PhaseFieldUnstructured::fluidEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
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
    factor = factor / pow(resolution, 2);
    phaseGradient(coords, el.idof, factor, e, g);
  }


  void PhaseFieldUnstructured::fluidEnergyAll(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    // Compute the energy for all N-1 phases
    double vol = el.parameters[0];

    // Get energy and gradient of the Nth fluid
    int nNodes = 7;
    vector<double> cN(nNodes, 1);
    vector<double> gN(nNodes);
    for (int iFluid=0; iFluid<nFluid-1; iFluid++) {
      for (int i=0; i<nNodes; i++) {
        cN[i] -= coords[el.idof[i] + iFluid];
      }
    }
    double bfactorN = 0.5 * kappa[nFluid-1] * vol;
    double gfactorN = 0.5 * kappaP[nFluid-1] * vol;
    if (e) *e += bfactorN * pow(cN[0], 2) * pow(cN[0]-1, 2);
    if (g) {
      gN[0] += bfactorN * 2 * cN[0] * (cN[0]-1) * (2*cN[0]-1);
      phaseGradient(cN, vec::iota(7), gfactorN, e, &gN);
    } else {
      phaseGradient(cN, vec::iota(7), gfactorN, e, nullptr);
    }

    // All other fluids
    for (int iFluid=0; iFluid<nFluid-1; iFluid++) {
      vector<int> idof = el.idof + iFluid;

      // Bulk energy
      double c = coords[idof[0]];
      double bfactor = 0.5 * kappa[iFluid] * vol;
      if (e) *e += bfactor * pow(c, 2) * pow(c-1, 2);
      if (g) (*g)[idof[0]] += bfactor * 2 * c * (c-1) * (2*c-1);

      // Gradient energy
      double gfactor = 0.5 * kappaP[iFluid] * vol;
      gfactor = gfactor / pow(resolution, 2);
      phaseGradient(coords, idof, gfactor, e, g);

      // Contribution to Nth fluid
      if (g) {
        for (int i=0; i<nNodes; i++) {
          (*g)[idof[i]] -= gN[i];
        }
      }
    }
  }


  void PhaseFieldUnstructured::fluidPairEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    double vol = el.parameters[0];
    int iPair = el.parameters[1];
    double c1 = coords[el.idof[0]];
    double c2 = coords[el.idof[1]];

    // Bulk energy
    double factor = 2 * kappa[iPair] * vol; // kappa = beta
    if (e) {
      auto quartic = [](double c){ return pow(c,2)*pow(c-1,2); };
      *e += factor * (quartic(c1) + quartic(c2) + quartic(c1+c2));
    }
    if (g) {
      auto gQuartic = [](double c){ return 2*c*(c-1)*(2*c-1); };
      (*g)[el.idof[0]] += factor * (gQuartic(c1) + gQuartic(c2) + gQuartic(c1+c2));
    }

    // Gradient energy
    double res2 = pow(resolution, 2);
    factor = -0.25 * kappaP[iPair] * vol / res2; // kappaP = -4 lambda
    phaseGradient(coords, el.idof, factor, e, g);
  }


  void PhaseFieldUnstructured::pressureEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    // Constant volume / pressure constraints rely upon the whole system
    int iDof = el.idof[0];
    int iFluid = fluidType[iDof];

    if (nFluid == 1) {
      double volume = 0.5*(coords[iDof]+1) * nodeVol[iDof];
      if (e) *e -= pressure[iFluid] * volume;
      if (g) (*g)[iDof] -= 0.5 * pressure[iFluid] * nodeVol[iDof];

    } else {
      double volume = coords[iDof] * nodeVol[iDof];
      if (e) *e -= pressure[iFluid] * volume;
      if (g) (*g)[iDof] -= pressure[iFluid] * nodeVol[iDof];
    }
  }


  void PhaseFieldUnstructured::densityConstraintEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
     // Total density soft constraint
     if (nFluid == 1) return;
     double coef = densityConst * surfaceTensionMean * pow(resolution, 2);
     double rhoDiff = coords[el.idof[0]] + coords[el.idof[1]] + coords[el.idof[2]] - 1;
     if (e) *e += coef * pow(rhoDiff, 2);
     if (g) {
       (*g)[el.idof[0]] += 2 * coef * rhoDiff;
       (*g)[el.idof[1]] += 2 * coef * rhoDiff;
       (*g)[el.idof[2]] += 2 * coef * rhoDiff;
     }
  }


  void PhaseFieldUnstructured::surfaceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    // Solid surface energy
    // parameter[0]: Surface area
    // parameter[1]: Wetting parameter 1/sqrt(2) cos(theta)
    if (nFluid == 1) {
      double phi = coords[el.idof[0]];
      if (e) *e += el.parameters[1] * (pow(phi,3)/3 - phi - 2.0/3) * el.parameters[0];
      if (g) (*g)[el.idof[0]] += el.parameters[1] * (pow(phi,2) - 1) * el.parameters[0];
    }
  }


  void PhaseFieldUnstructured::forceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
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
  }


  void PhaseFieldUnstructured::ffConfinementEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    // Confining potential for the frozen fluid method
    // Parameters:
    //   0: Confinement strength
    //   1: Initial concentration
    double c = coords[el.idof[0]];
    double strength = el.parameters[0];
    double c0 = el.parameters[1];
    double coef = strength * surfaceTensionMean * pow(resolution, 2);
    if ((c0-0.5)*(c-0.5) < 0) {
      if (e) *e += coef * pow(c-0.5, 2);
      if (g) (*g)[el.idof[0]] += coef * 2 * (c-0.5);
    }
  }


  void PhaseFieldUnstructured::elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    switch (el.type) {
      case FLUID_ENERGY:
        fluidEnergy(coords, el, e, g);
        break;
      case FLUID_ENERGY_ALL:
        fluidEnergyAll(coords, el, e, g);
        break;
      case FLUID_PAIR_ENERGY:
        fluidPairEnergy(coords, el, e, g);
        break;
      case PRESSURE_ENERGY:
        pressureEnergy(coords, el, e, g);
        break;
      case DENSITY_CONSTRAINT_ENERGY:
        densityConstraintEnergy(coords, el, e, g);
        break;
      case SURFACE_ENERGY:
        surfaceEnergy(coords, el, e, g);
        break;
      case FORCE_ENERGY:
        forceEnergy(coords, el, e, g);
        break;
      case FF_CONFINEMENT_ENERGY:
        ffConfinementEnergy(coords, el, e, g);
        break;
      default:
        throw std::invalid_argument("Unknown energy element type.");
    }
  }


  PhaseFieldUnstructured& PhaseFieldUnstructured::setGridSize(vector<int> gridSize) {
    this->gridSize = gridSize;
    this->nGrid = gridSize[0] * gridSize[1] * gridSize[2];
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setNFluid(int nFluid) {
    this->nFluid = nFluid;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setResolution(double resolution) {
    this->resolution = resolution;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setInterfaceSize(double interfaceSize) {
    this->interfaceSize = vector<double>(nFluid, interfaceSize);
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setInterfaceSize(vector<double> interfaceSize) {
    this->interfaceSize = interfaceSize;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setSurfaceTension(double surfaceTension) {
    this->surfaceTension = vector<double>(nFluid, surfaceTension);
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setSurfaceTension(vector<double> surfaceTension) {
    this->surfaceTension = surfaceTension;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setDensityConstraint(std::string method) {
    if (vec::isIn({"gradient","hard"}, method)) {
      densityConstraint = DENSITY_HARD;
    } else if (vec::isIn({"energy penalty","soft"}, method)) {
      densityConstraint = DENSITY_SOFT;
    } else if (method == "none") {
      densityConstraint = DENSITY_NONE;
    } else if (method == "fixed") {
      densityConstraint = DENSITY_FIXED;
    } else {
      throw std::invalid_argument("Invalid density constraint. Allowed methods are: fixed, none, gradient / hard, or energy penalty / soft");
    }
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setPressure(vector<double> pressure) {
    this->pressure = pressure;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setVolume(vector<double> volume, double volConst) {
    this->volumeFixed = true;
    this->volConst = volConst;
    this->volume = volume;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setVolumeFixed(bool volumeFixed, double volConst) {
    this->volumeFixed = volumeFixed;
    this->volConst = volConst;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setSolid(vector<bool> solid) {
    this->solid = solid;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setSolid(std::function<bool(int,int,int)> solidFn) {
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

  PhaseFieldUnstructured& PhaseFieldUnstructured::setContactAngle(vector<double> contactAngle) {
    this->contactAngle = contactAngle;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setContactAngle(std::function<double(int,int,int)> contactAngleFn) {
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

  PhaseFieldUnstructured& PhaseFieldUnstructured::setForce(vector<double> force, vector<int> iFluid) {
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

  PhaseFieldUnstructured& PhaseFieldUnstructured::setFixFluid(int iFluid, bool fix) {
    if ((int)fixFluid.size()!=nFluid) fixFluid = vector<bool>(nFluid, false);
    fixFluid[iFluid] = fix;
    return *this;
  }

  PhaseFieldUnstructured& PhaseFieldUnstructured::setConfinement(vector<double> strength) {
    this->confinementStrength = strength;
    return *this;
  }


  vector<int> PhaseFieldUnstructured::getCoord(int i) const {
    return minim::getCoord(i, gridSize);
  }


  int PhaseFieldUnstructured::getType(int i) const {
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


  void PhaseFieldUnstructured::setDefaults() {
    if (interfaceSize.empty()) interfaceSize = vector<double>(nFluid, resolution);
    if (solid.empty()) solid = vector<bool>(nGrid, false);
    if (pressure.empty()) pressure = vector<double>(nFluid, 0);
    if (fixFluid.empty()) fixFluid = vector<bool>(nFluid, false);
    if (confinementStrength.empty()) confinementStrength = vector<double>(nFluid, 0);
  }


  void PhaseFieldUnstructured::checkArraySizes() {
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


  vector<double> PhaseFieldUnstructured::diffuseSolid(vector<bool> solid, int iFluid, bool twoStep) {
    return diffuseSolid(solid, *this, iFluid, twoStep);
  }

  vector<double> PhaseFieldUnstructured::diffuseSolid(vector<bool> solid, vector<int> gridSize, int nFluid, int iFluid, bool twoStep) {
    PhaseFieldUnstructured potential;
    potential.setNFluid(nFluid);
    potential.setGridSize(gridSize);
    return diffuseSolid(solid, potential, iFluid, twoStep);
  }

  vector<double> PhaseFieldUnstructured::diffuseSolid(vector<bool> solid, PhaseFieldUnstructured potential, int iFluid, bool twoStep) {
    vector<double> confinement(potential.nFluid);
    confinement[iFluid] = 100;
    potential.setConfinement(confinement);
    potential.setDensityConstraint("none");

    double nGrid = potential.gridSize[0] * potential.gridSize[1] * potential.gridSize[2];
    vector<double> init(potential.nFluid * nGrid);
    for (int iGrid=0; iGrid<nGrid; iGrid++) {
      int i = iGrid * potential.nFluid + iFluid;
      init[i] = solid[iGrid];
    }
    State state(potential, init);
    state.convergence = 1e-8 * potential.surfaceTension[iFluid] * pow(potential.resolution, 2);

    Lbfgs min;
    min.setLinesearch("none");
    min.setMaxIter(100);
    auto minimum = min.minimise(state);

    if (twoStep) {
      state = State(potential, minimum);
      minimum = min.minimise(state);
    }

    return minimum;
  }

}
