#include "communicators/CommGrid.h"

#include "Potential.h"
#include "utils/vec.h"
#include "utils/mpi.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  inline int make1dIndex(const vector<int>& indices, const vector<int>& sizes) {
    int index = (indices[0] + sizes[0]) % sizes[0];
    for (int i=1; i<(int)indices.size(); i++) {
      int dimIdx = (indices[i] + sizes[i]) % sizes[i];
      index = index*sizes[i] + dimIdx;
    }
    return index;
  }

  inline vector<int> makeNdIndices(int index, const vector<int>& sizes) {
    int nDim = sizes.size();
    vector<int> indices(nDim);
    for (int iDim=nDim-1; iDim>0; iDim--) {
      indices[iDim] = index % sizes[iDim];
      index = index / sizes[iDim];
    }
    indices[0] = index;
    return indices;
  }

  // Recursive function to populate an N-dimensional sub-grid (out)
  // dim : current dimension iterating over
  // indices : current index of 'in' for all dimensions (passed by reference)
  // iOut : current index of sub-grid 'out' (passed by reference)
  // supSizes : the dimension sizes of the super-grid input
  // subSizes : the dimension sizes of the sub-grid output
  // subStart : the index on which the sub-grid starts on the super-grid for each dimension
  // in : The super-grid vector to copy from
  // out : The sub-grid vector to copy to
template <typename T>
  inline void assignNDimSubgrid(int dim, vector<int>& indices, int& iOut, const vector<int>& supSizes, const vector<int>& subSizes, const vector<int>& subStart, const vector<T>& in, vector<T>& out) {
    // Assign values once iterated through all dimensions
    if (dim == (int)subSizes.size()) {
      int iIn = make1dIndex(indices, supSizes);
      out[iOut++] = in[iIn];
      return;
    }

    // Iterate over this dimension and recursively go to next dimension
    int iEnd = subStart[dim] + subSizes[dim];
    for (indices[dim]=subStart[dim]; indices[dim]<iEnd; indices[dim]++) {
      assignNDimSubgrid(dim+1, indices, iOut, supSizes, subSizes, subStart, in, out);
    }
  }


  // Recursive function to populate part of an N-dimensional super-grid (out)
  // dim : current dimension iterating over
  // indices : current index of 'out' for all dimensions (passed by reference)
  // iIn : current index of sub-grid 'in' (passed by reference)
  // supSizes : the dimension sizes of the super-grid output
  // subSizes : the dimension sizes of the sub-grid input
  // subStart : the index on which the sub-grid starts on the super-grid for each dimension
  // in : The sub-grid vector to copy from
  // out : The super-grid vector to copy to
  template <typename T>
  inline void assignNDimSupergrid(int dim, vector<int>& indices, int& iIn, const vector<int>& supSizes, const vector<int>& subSizes, const vector<int>& subStart, const vector<T>& in, vector<T>& out) {
    // Assign values once iterated through all dimensions
    if (dim == (int)supSizes.size()) {
      int iOut = make1dIndex(indices, supSizes);
      out[iOut] = in[iIn++];
      return;
    }

    // Iterate over this dimension and recursively go to next dimension
    int iEnd = subStart[dim] + subSizes[dim];
    for (indices[dim]=subStart[dim]; indices[dim]<iEnd; indices[dim]++) {
      assignNDimSupergrid(dim+1, indices, iIn, supSizes, subSizes, subStart, in, out);
    }
  }


  //===== Assign data =====//

  vector<int> CommGrid::assignBlock(const vector<int>& in) const {
    return assignBlockImpl(in);
  }

  vector<bool> CommGrid::assignBlock(const vector<bool>& in) const {
    return assignBlockImpl(in);
  }

  vector<double> CommGrid::assignBlock(const vector<double>& in) const {
    return assignBlockImpl(in);
  }

  template <typename T>
  vector<T> CommGrid::assignBlockImpl(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (commSize == 1) return in;

    if (in.size() == ndof) return assignProcImpl<T>(in);

    vector<T> out(nproc);
    vector<int> start = blockStart - procStart;
    vector<int> indices(nDim);
    int iIn = 0;
    assignNDimSupergrid(0, indices, iIn, procSizes, blockSizes, start, in, out);
    return out;
  }
  template vector<int> CommGrid::assignBlockImpl(const vector<int>&) const;
  template vector<bool> CommGrid::assignBlockImpl(const vector<bool>&) const;
  template vector<double> CommGrid::assignBlockImpl(const vector<double>&) const;


  vector<int> CommGrid::assignProc(const vector<int>& in) const {
    return assignProcImpl(in);
  }

  vector<bool> CommGrid::assignProc(const vector<bool>& in) const {
    return assignProcImpl(in);
  }

  vector<double> CommGrid::assignProc(const vector<double>& in) const {
    return assignProcImpl(in);
  }

  template <typename T>
  vector<T> CommGrid::assignProcImpl(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (in.size() != ndof) throw std::invalid_argument("Input data has incorrect size. All degrees of freedom required.");
    if (commSize == 1) return in;

    vector<T> out(nproc);
    vector<int> indices(nDim);
    int iOut = 0;
    assignNDimSubgrid(0, indices, iOut, globalSizes, procSizes, procStart, in, out);
    return out;
  }
  template vector<int> CommGrid::assignProcImpl(const vector<int>&) const;
  template vector<bool> CommGrid::assignProcImpl(const vector<bool>&) const;
  template vector<double> CommGrid::assignProcImpl(const vector<double>&) const;


  //===== Access data =====//
  vector<int> CommGrid::getCoords(int loc) const {
    return makeNdIndices(loc, globalSizes);
  }

  int CommGrid::getBlock(int loc) const {
    vector<int> coords = makeNdIndices(loc, globalSizes);
    vector<int> blockIndices = coords / blockSizes;
    return make1dIndex(blockIndices, commArray);
  }


  int CommGrid::getLocalIdx(int loc, int block) const {
    vector<int> blockCoords(nDim);
    // Check if in this block, otherwise get local coordinates
    // if-else to optimise when block is provided
    if (block == -1) {
      vector<int> coords = makeNdIndices(loc, globalSizes);
      vector<int> blockIndices = coords / blockSizes;
      block = make1dIndex(blockIndices, commArray);
      if (commRank != block) return -1;
      blockCoords = coords - blockIndices * blockSizes;

    } else {
      if (commRank != block) return -1;
      vector<int> coords = makeNdIndices(loc, globalSizes);
      for (int i=0; i<nDim; i++) {
        blockCoords[i] = coords[i] % blockSizes[i];
      }
    }

    vector<int> procCoords = blockCoords + blockStart - procStart;
    return make1dIndex(procCoords, procSizes);
  }


  //===== Internal functions =====//


  //===== Constructors, etc =====//
  CommGrid::CommGrid(const CommGrid& other)
    : Communicator(other),
    nDim(other.nDim),
    haloWidth(other.haloWidth),
    commIndices(other.commIndices),
    globalSizes(other.globalSizes),
    blockSizes(other.blockSizes),
    blockStart(other.blockStart),
    procSizes(other.procSizes),
    procStart(other.procStart)
  {
    if (usesThisProc) makeMPITypes();
  }


  CommGrid& CommGrid::operator=(const CommGrid& other) {
    Communicator::operator=(other);
    nDim = other.nDim;
    haloWidth = other.haloWidth;
    globalSizes = other.globalSizes;
    blockSizes = other.blockSizes;
    blockStart = other.blockStart;
    procSizes = other.procSizes;
    procStart = other.procStart;
    commIndices = other.commIndices;
    if (usesThisProc) makeMPITypes();
    return *this;
  }

  std::unique_ptr<Communicator> CommGrid::clone() const {
    return std::make_unique<CommGrid>(static_cast<const CommGrid&>(*this));
  }


  //===== Setup =====//

  vector<int> primeFactorisation(int n) {
    // Return the prime factors in decreasing order
    vector<int> factors;
    int testFactor = (int)sqrt(n);
    while (testFactor > 1) {
      if (n % testFactor == 0) {
        factors.push_back(testFactor);
        n /= testFactor;
      } else {
        testFactor--;
      }
    }
    if (n > 1) factors.push_back(n);
    return factors;
  }


  vector<int> assignCommArray(int commSize, vector<int> globalSizes) {
    int nDim = globalSizes.size();
    vector<int> commArray(nDim, 1);
    // Compute the ideal number of procs for each dimension if the blocks where cubic (not necessarily integers)
    int nGridBlock = vec::product(globalSizes) / commSize;
    vector<double> idealSize = globalSizes / pow(nGridBlock, 1/nDim);
    // Compute the prime factors of the number of processors
    vector<int> factors = primeFactorisation(commSize);
    // Use a greedy allocation to choose the dimension that is brought closest to the ideal value
    for (int factor : factors) {
      int closestDim = 0;
      double closestRatio = 0;
      for (int iDim=0; iDim<nDim; iDim++) {
        // Measure the closeness by the 'absolute' ratio (<=1)
        double ratio = (commArray[iDim] * factor) / idealSize[iDim];
        if (ratio > 1) ratio = 1 / ratio;
        if (ratio > closestRatio) {
          closestRatio = ratio;
          closestDim = iDim;
        }
      }
      commArray[closestDim] *= factor;
    }
    // Could look into MPI_Dims_create?
    return commArray;
  }


  void CommGrid::setup(Potential& pot, size_t ndof, vector<int> ranks) {
    defaultSetup(ndof, ranks);
    if (mpi.size==1 || ranks.size()==1 || !usesThisProc) return; // If running serial
    setComm(ranks);

    // Set the globalSizes, comm array sizes, blockSizes, etc
    globalSizes = vector<int>(pot.gridSize.begin(), pot.gridSize.end());
    nDim = globalSizes.size();
    if (commArray.empty()) commArray = assignCommArray(commSize, globalSizes);
    commIndices = makeNdIndices(commRank, commArray);

    blockSizes = globalSizes / commArray;
    blockStart = blockSizes * commIndices;
    procSizes = blockSizes;
    procStart = blockStart;
    for (int i=0; i<nDim; i++) {
      if (commArray[i] > 1) {
        procSizes[i] += 2 * haloWidth;
        procStart[i] -= haloWidth;
      }
    }

    // Assign the base class values
    nblock = vec::product(blockSizes);
    nproc = vec::product(procSizes);

    // Assign the send receive and gather MPI types
    makeMPITypes();
  }


  // Create a vector of directions to neighbouring processors
  vector<vector<int>> getNeighbourDirections(const vector<int>& commArray) {
    int nDim = commArray.size();
    int nDir = pow(3, nDim);

    vector<int> ignoreDims;
    for (int iDim=0; iDim<nDim; iDim++) {
      if (commArray[iDim]==1) ignoreDims.push_back(iDim);
    }

    vector<vector<int>> directions;
    vector<int> dirArraySizes(nDim, 3);
    for (int iDir=0; iDir<nDir; iDir++) {
      if (iDir == (nDir-1)/2) continue; // Skip all zero direction
      // Convert 'iDir' to a ternary direction vector
      vector<int> direction = makeNdIndices(iDir, dirArraySizes) - 1;

      bool keep = true;
      for (int ignoreDim : ignoreDims) {
        if (direction[ignoreDim] != 0) keep = false;
      }
      if (keep) directions.push_back(direction);
    }
    return directions;
  }


  #ifdef PARALLEL
  // Create send and receive subarrays for each communication direction
  MPI_Datatype createSubarray(int type, vector<int> direction, vector<int> procSizes, vector<int> blockSizes, vector<int> commArray, int haloWidth) {
    MPI_Datatype subarray;
    int nDim = commArray.size();
    int sizes[nDim];
    int start[nDim];
    for (int iDim=0; iDim<nDim; iDim++) {
      if (direction[iDim] == -1) {
        sizes[iDim] = haloWidth;
        if (type == 0) start[iDim] = haloWidth; // send
        if (type == 1) start[iDim] = 0; // recv
      } else if (direction[iDim] == 1) {
        sizes[iDim] = haloWidth;
        if (type == 0) start[iDim] = procSizes[iDim] - 2*haloWidth; // send
        if (type == 1) start[iDim] = procSizes[iDim] - haloWidth; // recv
      } else if (direction[iDim] == 0) {
        sizes[iDim] = blockSizes[iDim];
        if (commArray[iDim] == 1) {
          start[iDim] = 0;
        } else {
          start[iDim] = haloWidth;
        }
      }
    }
    MPI_Type_create_subarray(nDim, &procSizes[0], sizes, start, MPI_ORDER_C, MPI_DOUBLE, &subarray);
    return subarray;
  }
  #endif


  // Make the MPI datatypes to send and receive from each proc
  void CommGrid::makeMPITypes() {
    if (commSize <= 1) return;

    #ifdef PARALLEL
    // Send and receive types
    auto directions = getNeighbourDirections(commArray);
    for (const auto& direction : directions) {
      // Create the MPI datatypes
      MPI_Datatype sendSubarray = createSubarray(0, direction, procSizes, blockSizes, commArray, haloWidth);
      MPI_Datatype recvSubarray = createSubarray(1, direction, procSizes, blockSizes, commArray, haloWidth);
      MPI_Type_commit(&sendSubarray);
      MPI_Type_commit(&recvSubarray);
      // Create the communication objects
      int iNeighbour = make1dIndex(commIndices + direction, commArray);
      int sendTag = make1dIndex(direction+1, vector<int>(nDim, 3));
      int recvTag = make1dIndex(-direction+1, vector<int>(nDim, 3));
      sendTypes.push_back({iNeighbour, sendTag, sendSubarray});
      recvTypes.push_back({iNeighbour, recvTag, recvSubarray});
    }

    // Types for gather
    nGather = vector<int>(commSize, 1);
    iGather = vector<int>(commSize);
    for (int iComm=0; iComm<commSize; iComm++) {
      vector<int> commIndices = makeNdIndices(iComm, commArray);
      vector<int> blockStart = blockSizes * commIndices;
      iGather[iComm] = make1dIndex(blockStart, globalSizes);
    }
    // Local block type for sending
    vector<int> start = blockStart - procStart;
    MPI_Type_create_subarray(nDim, &procSizes[0], &blockSizes[0], &start[0], MPI_ORDER_C, MPI_DOUBLE, &blockType);
    MPI_Type_commit(&blockType);
    // Gather type for receiving
    // Make subarray datatype with 0 displacement (displacement is controlled by iGather in Gatherv)
    start = vector<int>(nDim, 0);
    MPI_Datatype gatherBlockType;
    MPI_Type_create_subarray(nDim, &globalSizes[0], &blockSizes[0], &start[0], MPI_ORDER_C, MPI_DOUBLE, &gatherBlockType);
    // Change extent to one double so displacements are in correct units and so gather works properly with overlapping data (I think?)
    MPI_Type_create_resized(gatherBlockType, 0, 1*sizeof(double), &gatherType);
    MPI_Type_commit(&gatherType);
    #endif
  }

}
