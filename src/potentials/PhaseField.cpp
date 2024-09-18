#include "potentials/PhaseField.h"

#include <array>
#include <math.h>
#include <stdexcept>
#include <functional>
#include "State.h"
#include "utils/vec.h"
#include "utils/range.h"
#include "communicators/CommGrid.h"
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


  static vector<int> getCoord(int i, vector<int> gridSize) {
    int z = i % gridSize[2];
    int y = (i - z) / gridSize[2] % gridSize[1];
    int x = i / (gridSize[1] * gridSize[2]);
    return {x, y, z};
  }

  static int getIdx(const vector<int>& coord, const vector<int>& gridSize) {
    int x = (coord[0] + gridSize[0]) % gridSize[0];
    int y = (coord[1] + gridSize[1]) % gridSize[1];
    int z = (coord[2] + gridSize[2]) % gridSize[2];
    return (x*gridSize[1] + y)*gridSize[2] + z;
  }

  static vector2d<int> getNeighbours(const vector<int>& xGrid, const vector<int>& gridSizes) {
    vector2d<int> neighbours(3);
    vector<int> x = xGrid;
    for (int iDim=0; iDim<3; ++iDim) {
      neighbours[iDim] = vector<int>(2);
      int x0 = x[iDim];
      x[iDim] = x0 - 1;
      neighbours[iDim][0] = getIdx(x, gridSizes);
      x[iDim] = x0 + 1;
      neighbours[iDim][1] = getIdx(x, gridSizes);
      x[iDim] = x0;
    }
    return neighbours;
  }

  static vector2d<int> dx = {{
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


  void PhaseField::assignFluidCoefficients() {
    int nParams = (model==MODEL_BASIC) ? nFluid : 0.5*nFluid*(nFluid-1);

    // Ensure surface tension and interface widths are the correct size
    if ((int)interfaceSize.size() == 1) {
      interfaceSize = vector<double>(nParams, interfaceSize[0]);
    } else if ((int)interfaceSize.size() != nParams) {
      throw std::invalid_argument("PhaseField: Invalid size of interfaceSize array.");
    }
    if ((int)surfaceTension.size() == 1) {
      surfaceTension = vector<double>(nParams, surfaceTension[0]);
    } else if ((int)surfaceTension.size() != nParams) {
      throw std::invalid_argument("PhaseField: Invalid size of surfaceTension array.");
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


  void PhaseField::init(const vector<double>& coords) {
    if ((int)coords.size() != nGrid*nFluid) {
      throw std::invalid_argument("PhaseField: Size of coordinates array does not match the grid size.");
    }

    setDefaults();
    checkArraySizes();
    assignFluidCoefficients();
    this->convergence = 1e-8 * surfaceTensionMean * pow(resolution, 2);
    if (densityConstraint == DENSITY_FIXED) fixFluid[nFluid-1] = true;

    // Forces
    fMag = vector<double>(nFluid);
    fNorm = vector2d<double>(nFluid);
    if (!force.empty()) {
      for (int i=0; i<nFluid; i++) {
        fMag[i] = vec::norm(force[i]);
        fNorm[i] = force[i] / fMag[i];
      }
    }
  }


  void PhaseField::initLocal(const vector<double>& coords, const Communicator& comm) {
    // Update local grid size
    auto commGrid = static_cast<const CommGrid&>(comm);
    procSizes = vector<int>(3);
    procStart = vector<int>(3);
    haloWidths = vector<int>(3);
    for (int iDim=0; iDim<3; iDim++) {
      // Do not copy the last dimension (fluid no.), store only the grid size
      procSizes[iDim] = commGrid.procSizes[iDim];
      procStart[iDim] = commGrid.procStart[iDim];
      haloWidths[iDim] = commGrid.haloWidths[iDim];
    }
    nGrid = vec::product(procSizes);
    if (!comm.usesThisProc) return;

    // Distribute global parameters
    vector<double> coordsLocal = comm.assignProc(coords);
    if (!contactAngle.empty()) contactAngle = comm.assignProc(contactAngle);
    if ((int)solid.size() != nGrid) {
      vector<bool> solidGlobal = solid;
      solid = vector<bool>(nGrid);
      int iLocal = 0;
      for (int x=procStart[0]; x<procStart[0]+procSizes[0]; x++) {
        for (int y=procStart[1]; y<procStart[1]+procSizes[1]; y++) {
          for (int z=procStart[2]; z<procStart[2]+procSizes[2]; z++) {
            int iGlobal = getIdx({x,y,z}, gridSize);
            solid[iLocal++] = solidGlobal[iGlobal];
          }
        }
      }
    }

    // Get fluid volume and solid surface area for each node (not in halo)
    nodeVol = vector<double>(nGrid, 0);
    surfaceArea = vector<double>(nGrid, 0);
    for (int iGrid : RangeI(procSizes, haloWidths)) {
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

    neighbours.resize(nGrid);
    for (int iGrid : RangeI(procSizes, haloWidths)) {
      vector<int> xGrid = getCoord(iGrid, procSizes);
      neighbours[iGrid] = getNeighbours(xGrid, procSizes);
    }

    // Set initial volumes for constant volume constraint
    if (volume.empty() && volumeFixed) {
      volume = vector<double>(nFluid, 0);
      for (int iGrid : RangeI(procSizes, haloWidths)) {
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          volume[iFluid] += coordsLocal[iFluid+nFluid*iGrid] * nodeVol[iGrid];
        }
      }
      for (int iFluid=0; iFluid<nFluid; iFluid++) {
        volume[iFluid] = comm.sum(volume[iFluid]);
      }
      totalVolume = 0;
      for (int iGrid : RangeI(procSizes, haloWidths)) {
        totalVolume += nodeVol[iGrid];
      }
      totalVolume = comm.sum(totalVolume);
    }

    // Set initial values for confinement potential
    if (vec::any(confinementStrength)) ffInit = coordsLocal;
  }


  void PhaseField::phaseGradient(const vector<double>& coords, int iGrid, int iFluid, const vector<int>& xGrid,
                                 const vector2d<int>& neighbours, double factor, double* e, vector<double>* g) const {
    int i0 = iGrid * nFluid + iFluid;
    double c1 = coords[i0];

    double grad2 = 0;
    for (int iDir=0; iDir<3; iDir++) {
      if (procSizes[iDir] == 1) continue;

      int imGrid = neighbours[iDir][0];
      int ipGrid = neighbours[iDir][1];
      int im = imGrid * nFluid + iFluid;
      int ip = ipGrid * nFluid + iFluid;

      if (!solid[imGrid] && !solid[ipGrid]) { // No solid
        double gradm = c1 - coords[im];
        double gradp = c1 - coords[ip];
        if (e) grad2 += 0.5 * (pow(gradm,2) + pow(gradp,2));
        if (g) {
          (*g)[i0] += factor * (gradm + gradp);
          (*g)[im] -= factor * gradm;
          (*g)[ip] -= factor * gradp;
        }

      } else if (!solid[ipGrid]) { // Solid on negative side
        double gradp = c1 - coords[ip];
        if (e) grad2 += pow(gradp,2);
        if (g) {
          (*g)[i0] += factor * 2*gradp;
          (*g)[ip] -= factor * 2*gradp;
        }

      } else if (!solid[imGrid]) { // Solid on positive side
        double gradm = c1 - coords[im];
        if (e) grad2 += pow(gradm,2);
        if (g) {
          (*g)[i0] += factor * 2*gradm;
          (*g)[im] -= factor * 2*gradm;
        }
      }
    }

    if (e) *e += factor * grad2;
  }


  inline void PhaseField::phasePairGradient(const vector<double>& coords, int iGrid, int iFluid1, int iFluid2, const vector<int>& xGrid,
                                            const vector2d<int>& neighbours, double factor, double* e, vector<double>* g) const {
    int i01 = iGrid * nFluid + iFluid1;
    int i02 = iGrid * nFluid + iFluid2;
    double c1 = coords[i01];
    double c2 = coords[i02];

    double grad2 = 0;
    for (int iDir=0; iDir<3; iDir++) {
      if (procSizes[iDir] == 1) continue;

      int imGrid = neighbours[iDir][0];
      int ipGrid = neighbours[iDir][1];
      int im1 = imGrid * nFluid + iFluid1;
      int ip1 = ipGrid * nFluid + iFluid1;
      int im2 = imGrid * nFluid + iFluid2;
      int ip2 = ipGrid * nFluid + iFluid2;

      if (!solid[imGrid] && !solid[ipGrid]) { // No solid
        double gradm1 = c1 - coords[im1];
        double gradm2 = c2 - coords[im2];
        double gradp1 = c1 - coords[ip1];
        double gradp2 = c2 - coords[ip2];
        if (e) grad2 += 0.5 * (gradm1*gradm2 + gradp1*gradp2);
        if (g) {
          (*g)[i01] += factor * (gradm2 + gradp2);
          (*g)[i02] += factor * (gradm1 + gradp1);
          (*g)[im1] -= factor * gradm2;
          (*g)[im2] -= factor * gradm1;
          (*g)[ip1] -= factor * gradp2;
          (*g)[im2] -= factor * gradm1;
        }

      } else if (!solid[ipGrid]) { // Solid on negative side
        double gradp1 = c1 - coords[ip1];
        double gradp2 = c2 - coords[ip2];
        if (e) grad2 += gradp1 * gradp2;
        if (g) {
          (*g)[i01] += factor * 2*gradp2;
          (*g)[i02] += factor * 2*gradp1;
          (*g)[ip1] -= factor * 2*gradp2;
          (*g)[ip2] -= factor * 2*gradp1;
        }

      } else if (!solid[imGrid]) { // Solid on positive side
        double gradm1 = c1 - coords[im1];
        double gradm2 = c2 - coords[im2];
        if (e) grad2 += gradm1 * gradm2;
        if (g) {
          (*g)[i01] += factor * 2*gradm2;
          (*g)[i02] += factor * 2*gradm1;
          (*g)[im1] -= factor * 2*gradm2;
          (*g)[im2] -= factor * 2*gradm1;
        }
      }
    }

    if (e) *e += factor * grad2;
  }

  void PhaseField::fluidEnergy(const vector<double>& coords, int iGrid, const vector<int>& xGrid, double* e, vector<double>* g) const {
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      int iDof = iGrid * nFluid + iFluid;
      double c = coords[iDof];

      // Bulk energy
      if (nFluid == 1) {
        double factor = kappa[0] / 16 * nodeVol[iGrid];
        if (e) *e += factor * pow(c+1, 2) * pow(c-1, 2);
        if (g) (*g)[iDof] += factor * 4 * c * (c*c - 1);
      } else {
        double factor = 0.5 * kappa[iFluid] * nodeVol[iGrid];
        if (e) *e += factor * pow(c, 2) * pow(c-1, 2);
        if (g) (*g)[iDof] += factor * 2 * c * (c-1) * (2*c-1);
      }

      // Gradient energy
      double factor = (nFluid==1) ? 0.25*kappaP[0]*nodeVol[iGrid] : 0.5*kappaP[iFluid]*nodeVol[iGrid];
      factor = factor / pow(resolution, 2);
      phaseGradient(coords, iGrid, iFluid, xGrid, neighbours[iGrid], factor, e, g);
    }
  }


  void PhaseField::fluidPairEnergy(const vector<double>& coords, int iGrid, const vector<int>& xGrid, double* e, vector<double>* g) const {
    int iPair = 0;
    for (int iFluid1=0; iFluid1<nFluid; iFluid1++) {
      for (int iFluid2=iFluid1+1; iFluid2<nFluid; iFluid2++) {
        int iDof1 = iGrid * nFluid + iFluid1;
        int iDof2 = iGrid * nFluid + iFluid2;
        double c1 = coords[iDof1];
        double c2 = coords[iDof2];

        // Bulk energy
        double factor = 2 * kappa[iPair] * nodeVol[iGrid]; // kappa = beta
        if (e) {
          auto quartic = [](double c){ return pow(c,2)*pow(c-1,2); };
          *e += factor * (quartic(c1) + quartic(c2) + quartic(c1+c2));
        }
        if (g) {
          auto gQuartic = [](double c){ return 2*c*(c-1)*(2*c-1); };
          double gQ12 = gQuartic(c1 + c2);
          (*g)[iDof1] += factor * (gQuartic(c1) + gQ12);
          (*g)[iDof2] += factor * (gQuartic(c2) + gQ12);
        }

        // Gradient energy
        double res2 = pow(resolution, 2);
        factor = -0.25 * kappaP[iPair] * nodeVol[iGrid] / res2; // kappaP = -4 lambda
        phasePairGradient(coords, iGrid, iFluid1, iFluid2, xGrid, neighbours[iGrid], factor, e, g);

        iPair++;
      }
    }
  }


  void PhaseField::pressureEnergy(const vector<double>& coords, int iGrid, double* e, vector<double>* g) const {
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (pressure[iFluid]==0 || fixFluid[iFluid]) continue;

      if (nFluid == 1) {
        double volume = 0.5*(coords[iGrid]+1) * nodeVol[iGrid];
        if (e) *e -= pressure[iFluid] * volume;
        if (g) (*g)[iGrid] -= 0.5 * pressure[iFluid] * nodeVol[iGrid];

      } else {
        int iDof = iGrid * nFluid + iFluid;
        double volume = coords[iDof] * nodeVol[iGrid];
        if (e) *e -= pressure[iFluid] * volume;
        if (g) (*g)[iDof] -= pressure[iFluid] * nodeVol[iGrid];
      }
    }
  }


  void PhaseField::densityConstraintEnergy(const vector<double>& coords, int iGrid, double* e, vector<double>* g) const {
    if (nFluid==1 || densityConstraint!=DENSITY_SOFT) return;

    // Total density soft constraint
    double coef = densityConst * surfaceTensionMean * pow(resolution, 2);
    double rhoDiff = -1;
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      rhoDiff += coords[iGrid*nFluid+iFluid];
    }

    if (e) *e += coef * pow(rhoDiff, 2);

    if (!g) return;
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      (*g)[iGrid*nFluid+iFluid] += 2 * coef * rhoDiff;
    }
  }


  void PhaseField::surfaceEnergy(const vector<double>& coords, int iGrid, double* e, vector<double>* g) const {
    if (surfaceArea[iGrid]==0 || contactAngle.empty() || contactAngle[iGrid]==90) return;
    if (nFluid > 1) return; // Currently only implemented for binary fluids

    double wettingParam = 1/sqrt(2.0) * cos(contactAngle[iGrid] * 3.1415926536/180);
    double phi = coords[iGrid];
    if (e) *e += wettingParam * (pow(phi,3)/3 - phi - 2.0/3) * surfaceArea[iGrid];
    if (g) (*g)[iGrid] += wettingParam * (pow(phi,2) - 1) * surfaceArea[iGrid];
  }


  void PhaseField::forceEnergy(const vector<double>& coords, int iGrid, const vector<int>& xGrid, double* e, vector<double>* g) const {
    if (solid[iGrid] || !vec::any(fMag)) return;

    // Get the global Cartesian coordinates
    vector<int> coordI = xGrid + procStart;
    vector<double> coord{coordI[0]-(gridSize[0]-1)/2.0, coordI[1]-(gridSize[1]-1)/2.0, coordI[2]-(gridSize[2]-1)/2.0};

    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (fMag[iFluid]==0) continue;

      // Get liquid concentration
      int iDof = iGrid * nFluid + iFluid;
      double c = (nFluid==1) ? 0.5*(1+coords[iDof]) : coords[iDof];
      if (c < 0.01) continue; // So gradient is zero in bulk gas phase

      // Get the height
      double h = - vec::dotProduct(coord, fNorm[iFluid]) * resolution;
      double ePerConc = nodeVol[iGrid] * fMag[iFluid] * h;

      // Energy and gradient
      if (e) *e += c * ePerConc;
      if (!g) return;
      if (nFluid == 1) {
        (*g)[iDof] += 0.5 * ePerConc;
      } else {
        (*g)[iDof] += ePerConc;
      }
    }
  }


  void PhaseField::ffConfinementEnergy(const vector<double>& coords, int iGrid, double* e, vector<double>* g) const {
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (confinementStrength[iFluid] == 0) continue;

      int iDof = iGrid*nFluid+iFluid;
      double coef = confinementStrength[iFluid] * surfaceTensionMean * pow(resolution, 2);
      double c = coords[iDof];
      double c0 = ffInit[iDof];
      if ((c0-0.5)*(c-0.5) < 0) {
        if (e) *e += coef * pow(c-0.5, 2);
        if (g) (*g)[iDof] += coef * 2 * (c-0.5);
      }
    }
  }


  void PhaseField::volumeConstraintEnergy(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {
    // Get the (local) fluid volumes
    vector<double> volFluid(nFluid, 0);
    for (int iGrid : RangeI(procSizes, haloWidths)) {
      if (nFluid == 1) {
        volFluid[0] += 0.5*(coords[iGrid]+1) * nodeVol[iGrid];
      } else {
        for (int iFluid=0; iFluid<nFluid; iFluid++) {
          int iDof = iGrid*nFluid+iFluid;
          volFluid[iFluid] += coords[iDof] * nodeVol[iGrid];
        }
      }
    }

    // Compute the energy and gradient
    double volCoef = volConst * surfaceTensionMean / pow(resolution, 4);
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (fixFluid[iFluid]) continue;
      double volDiff = comm.sum(volFluid[iFluid]) - volume[iFluid];
      if (e) *e += volCoef * pow(volDiff, 2) / comm.size();
      if (!g) continue;
      for (int iGrid : RangeI(procSizes, haloWidths)) {
        int iDof = iGrid*nFluid + iFluid;
        double interfaceWeight = std::max(0.0, 4*coords[iDof]*(1-coords[iDof])); // Only apply the force to the interface nodes
        (*g)[iDof] += volCoef * volDiff * nodeVol[iGrid] * interfaceWeight;
      }
    }
  }


  void PhaseField::applyConstraints(const vector<double>& coords, const Communicator& comm, vector<double>* g) const {
    // Fixed fluid
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (!fixFluid[iFluid]) continue;
      for (int iGrid=0; iGrid<nGrid; iGrid++) {
        (*g)[iGrid*nFluid+iFluid] = 0;
      }
    }

    // Get a list of the non-fixed fluids
    vector<int> iVariableFluid;
    for (int iFluid=0; iFluid<nFluid; iFluid++) {
      if (!fixFluid[iFluid]) iVariableFluid.push_back(iFluid);
    }

    // Hard density constraint
    if (nFluid>1 && densityConstraint==DENSITY_HARD) {
      double normFactor = 1.0 / iVariableFluid.size();

      for (int iGrid=0; iGrid<nGrid; iGrid++) {
        // Dot product to get component of increasing density
        double component = 0;
        for (int iFluid : iVariableFluid) {
          component += (*g)[iGrid*nFluid+iFluid];
        }
        // Normalise to get correct corrections
        component *= normFactor;
        // Remove the component
        for (int iFluid : iVariableFluid) {
          (*g)[iGrid*nFluid+iFluid] -= component;
        }
      }
    }

    // Fixed volume
    if (volumeFixed) {
      // Dot product for components of gradient in direction of increasing volumes
      vector<double> component(nFluid, 0);
      for (int iGrid : RangeI(procSizes, haloWidths)) {
        for (int iFluid : iVariableFluid) {
          component[iFluid] += (*g)[iGrid*nFluid+iFluid] * sqrt(nodeVol[iGrid]);
        }
      }
      // Finish dot product + normalise to get correct corrections
      for (int iFluid : iVariableFluid) {
        component[iFluid] = comm.sum(component[iFluid]) / totalVolume;
      }
      // Remove the component
      for (int iGrid : RangeI(procSizes, haloWidths)) {
        for (int iFluid : iVariableFluid) {
          (*g)[iGrid*nFluid+iFluid] -= component[iFluid] * sqrt(nodeVol[iGrid]);
        }
      }
    }
  }


  void PhaseField::energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {
    if (g) *g = vector<double>(coords.size());

    for (int x=haloWidths[0]; x<procSizes[0]-haloWidths[0]; x++) {
      for (int y=haloWidths[1]; y<procSizes[1]-haloWidths[1]; y++) {
        for (int z=haloWidths[2]; z<procSizes[2]-haloWidths[2]; z++) {
          vector<int> xGrid{x, y, z};
          int iGrid = getIdx(xGrid, procSizes);

          if (model == MODEL_BASIC) {
            fluidEnergy(coords, iGrid, xGrid, e, g);
          } else if (model == MODEL_NCOMP) {
            fluidPairEnergy(coords, iGrid, xGrid, e, g);
          }

          surfaceEnergy(coords, iGrid, e, g);
          pressureEnergy(coords, iGrid, e, g);
          densityConstraintEnergy(coords, iGrid, e, g);
          forceEnergy(coords, iGrid, xGrid, e, g);
          ffConfinementEnergy(coords, iGrid, e, g);
        }
      }
    }

    // if (volumeFixed) volumeConstraintEnergy(coords, comm, e, g);

    if (g) applyConstraints(coords, comm, g);
  }


  std::map<std::string,vector<double>> PhaseField::energyComponents(const vector<double>& coords, const Communicator& comm) const {
    vector<std::string> components = {"fluid", "surface", "pressure", "density constraint", "force", "confinement", "volume constraint"};
    std::map<std::string,double> e;
    std::map<std::string,vector<double>> g;
    for(const auto &component : components) {
      e[component] = 0;
      g[component] = vector<double>(coords.size());
    }


    for (int x=haloWidths[0]; x<procSizes[0]-haloWidths[0]; x++) {
      for (int y=haloWidths[1]; y<procSizes[1]-haloWidths[1]; y++) {
        for (int z=haloWidths[2]; z<procSizes[2]-haloWidths[2]; z++) {
          vector<int> xGrid{x, y, z};
          int iGrid = getIdx(xGrid, procSizes);

          if (model == MODEL_BASIC) {
            fluidEnergy(coords, iGrid, xGrid, &e["fluid"], &g["fluid"]);
          } else if (model == MODEL_NCOMP) {
            fluidPairEnergy(coords, iGrid, xGrid, &e["fluid"], &g["fluid"]);
          }

          surfaceEnergy(coords, iGrid, &e["surface"], &g["surface"]);
          pressureEnergy(coords, iGrid, &e["pressure"], &g["pressure"]);
          densityConstraintEnergy(coords, iGrid, &e["density constraint"], &g["density constraint"]);
          forceEnergy(coords, iGrid, xGrid, &e["force"], &g["force"]);
          ffConfinementEnergy(coords, iGrid, &e["confinement"], &g["confinement"]);
        }
      }
    }

    // if (volumeFixed) volumeConstraintEnergy(coords, comm, &e["volume constraint"], &g["volume constraint"]);

    std::map<std::string,vector<double>> eg;
    for(const auto &component : components) {
      eg[component] = {comm.sum(e[component]), vec::norm(comm.gather(g[component]))};
    }
    return eg;
  }


  PhaseField& PhaseField::setGridSize(vector<int> gridSize) {
    this->gridSize = gridSize;
    this->nGrid = gridSize[0] * gridSize[1] * gridSize[2];
    return *this;
  }

  PhaseField& PhaseField::setNFluid(int nFluid) {
    this->nFluid = nFluid;
    this->dofPerNode = nFluid;
    return *this;
  }

  PhaseField& PhaseField::setResolution(double resolution) {
    this->resolution = resolution;
    return *this;
  }

  PhaseField& PhaseField::setInterfaceSize(double interfaceSize) {
    this->interfaceSize = vector<double>(nFluid, interfaceSize);
    return *this;
  }

  PhaseField& PhaseField::setInterfaceSize(vector<double> interfaceSize) {
    this->interfaceSize = interfaceSize;
    return *this;
  }

  PhaseField& PhaseField::setSurfaceTension(double surfaceTension) {
    this->surfaceTension = vector<double>(nFluid, surfaceTension);
    return *this;
  }

  PhaseField& PhaseField::setSurfaceTension(vector<double> surfaceTension) {
    this->surfaceTension = surfaceTension;
    return *this;
  }

  PhaseField& PhaseField::setDensityConstraint(std::string method) {
    if (vec::isIn({"gradient","hard"}, method)) {
      densityConstraint = DENSITY_HARD;
    } else if (vec::isIn({"energy penalty","soft"}, method)) {
      densityConstraint = DENSITY_SOFT;
    } else if (method == "none") {
      densityConstraint = DENSITY_NONE;
    } else if (method == "fixed") {
      densityConstraint = DENSITY_FIXED;
    } else {
      throw std::invalid_argument("PhaseField: Invalid density constraint. Allowed methods are: fixed, none, gradient / hard, or energy penalty / soft");
    }
    return *this;
  }

  PhaseField& PhaseField::setPressure(vector<double> pressure) {
    this->pressure = pressure;
    return *this;
  }

  PhaseField& PhaseField::setVolume(vector<double> volume, double volConst) {
    this->volumeFixed = true;
    this->volConst = volConst;
    this->volume = volume;
    return *this;
  }

  PhaseField& PhaseField::setVolumeFixed(bool volumeFixed, double volConst) {
    this->volumeFixed = volumeFixed;
    this->volConst = volConst;
    return *this;
  }

  PhaseField& PhaseField::setSolid(vector<bool> solid) {
    this->solid = solid;
    return *this;
  }

  PhaseField& PhaseField::setSolid(std::function<bool(int,int,int)> solidFn) {
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

  PhaseField& PhaseField::setContactAngle(double contactAngle) {
    this->contactAngle = vector<double>(nGrid, contactAngle);
    return *this;
  }

  PhaseField& PhaseField::setContactAngle(vector<double> contactAngle) {
    this->contactAngle = contactAngle;
    return *this;
  }

  PhaseField& PhaseField::setContactAngle(std::function<double(int,int,int)> contactAngleFn) {
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

  PhaseField& PhaseField::setForce(vector<double> force, vector<int> iFluid) {
    if ((int)force.size() != 3) throw std::invalid_argument("PhaseField: Invalid size of force array.");
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

  PhaseField& PhaseField::setFixFluid(int iFluid, bool fix) {
    if ((int)fixFluid.size()!=nFluid) fixFluid = vector<bool>(nFluid, false);
    fixFluid[iFluid] = fix;
    return *this;
  }

  PhaseField& PhaseField::setConfinement(vector<double> strength) {
    this->confinementStrength = strength;
    return *this;
  }


  int PhaseField::getType(int i) const {
    // Types:
    // -1: Solid (▪▪)   0: Bulk fluid (  )
    //           (▪▪)                 (  )
    // 1: ▫    2: ▫▫   3: ▫▫   4: ▫▫   5: ▪▫   6: ▪▪   7: ▪▪
    //                    ▫       ▫▫      ▫▫      ▫▫      ▪▫
    if (solid[i]) return -1;

    // Get which neighbouring nodes are solid
    std::array<bool,26> neiSolid;
    auto x0 = getCoord(i, procSizes);
    for (int iDir=0; iDir<26; iDir++) {
      auto xNei = x0 + dx[iDir];
      int iNei = getIdx(xNei, procSizes);
      neiSolid[iDir] = solid[iNei];
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
    throw std::runtime_error("PhaseField: Undefined surface type");
  }


  void PhaseField::setDefaults() {
    if (interfaceSize.empty()) interfaceSize = vector<double>(nFluid, resolution);
    if (solid.empty()) solid = vector<bool>(nGrid, false);
    if (pressure.empty()) pressure = vector<double>(nFluid, 0);
    if (fixFluid.empty()) fixFluid = vector<bool>(nFluid, false);
    if (confinementStrength.empty()) confinementStrength = vector<double>(nFluid, 0);
  }


  void PhaseField::checkArraySizes() {
    if ((int)interfaceSize.size() != nFluid) {
      throw std::invalid_argument("PhaseField: Invalid size interfaceSize array.");
    }
    if ((int)solid.size() != nGrid) {
      throw std::invalid_argument("PhaseField: Invalid size of solid array.");
    }
    if ((int)contactAngle.size() != nGrid*nFluid && !contactAngle.empty()) {
      throw std::invalid_argument("PhaseField: Invalid size of contactAngle array.");
    }
    if ((int)force.size() != nFluid && !force.empty()) {
      throw std::invalid_argument("PhaseField: Invalid size of force array.");
    }
    if ((int)volume.size() != nFluid && !volume.empty()) {
      throw std::invalid_argument("PhaseField: Invalid size of volume array.");
    }
    if ((int)pressure.size() != nFluid) {
      throw std::invalid_argument("PhaseField: Invalid size of pressure array.");
    }
    if ((int)fixFluid.size() != nFluid) {
      throw std::invalid_argument("PhaseField: Invalid size of fixed fluid array.");
    }
    if ((int)confinementStrength.size() != nFluid) {
      throw std::invalid_argument("PhaseField: Invalid size of confinement strength array.");
    }
  }


  vector<double> PhaseField::diffuseSolid(vector<bool> solid, int iFluid, bool twoStep) {
    return diffuseSolid(solid, *this, iFluid, twoStep);
  }

  vector<double> PhaseField::diffuseSolid(vector<bool> solid, vector<int> gridSize, int nFluid, int iFluid, bool twoStep) {
    PhaseField potential;
    potential.setNFluid(nFluid);
    potential.setGridSize(gridSize);
    return diffuseSolid(solid, potential, iFluid, twoStep);
  }

  vector<double> PhaseField::diffuseSolid(vector<bool> solid, PhaseField potential, int iFluid, bool twoStep) {
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
