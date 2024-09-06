#include "communicators/CommUnstructured.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <set>
//#include <limits>
#include <stdexcept>
#include "Potential.h"
#include "utils/vec.h"
#include "utils/mpi.h"
#include "utils/print.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  //===== Assign data =====//
  vector<int> CommUnstructured::assignBlock(const vector<int>& in) const {
    return assignBlockImpl(in);
  }

  vector<bool> CommUnstructured::assignBlock(const vector<bool>& in) const {
    return assignBlockImpl(in);
  }

  vector<double> CommUnstructured::assignBlock(const vector<double>& in) const {
    return assignBlockImpl(in);
  }

  template <typename T>
  vector<T> CommUnstructured::assignBlockImpl(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (commSize == 1) return in;

    vector<T> out(nblock);
    int i0 = (in.size() == ndof) ? iblock : 0; // Potential issue here if in.size() == nproc == ndof
    for (size_t i=0; i<nblock; i++) {
      out[i] = in[i0+i];
    }
    return out;
  }
  template vector<int> CommUnstructured::assignBlockImpl(const vector<int>&) const;
  template vector<bool> CommUnstructured::assignBlockImpl(const vector<bool>&) const;
  template vector<double> CommUnstructured::assignBlockImpl(const vector<double>&) const;


  vector<int> CommUnstructured::assignProc(const vector<int>& in) const {
    return assignProcImpl(in);
  }

  vector<bool> CommUnstructured::assignProc(const vector<bool>& in) const {
    return assignProcImpl(in);
  }

  vector<double> CommUnstructured::assignProc(const vector<double>& in) const {
    return assignProcImpl(in);
  }

  template <typename T>
  vector<T> CommUnstructured::assignProcImpl(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (in.size() != ndof) throw std::invalid_argument("CommUnstructured: Input data has incorrect size. All degrees of freedom required.");
    if (commSize == 1) return in;

    vector<T> out(nproc);
    // Assign the main blocks
    for (size_t i=0; i<nblock; i++) {
      out[i] = in[iblock+i];
    }
    // Assign the halo regions
    for (int i=0; i<commSize; i++) {
      for (int j=0; j<nrecv[i]; j++) {
        out[irecv[i]+j] = in[recv_lists[i][j]];
      }
    }
    return out;
  }
  template vector<int> CommUnstructured::assignProcImpl(const vector<int>&) const;
  template vector<bool> CommUnstructured::assignProcImpl(const vector<bool>&) const;
  template vector<double> CommUnstructured::assignProcImpl(const vector<double>&) const;


  //===== Access data =====//
  int CommUnstructured::getBlock(int loc) const {
    for (int i=commSize-1; i>=0; i--) {
      if (loc >= iblocks[i]) return i;
    }
    throw std::invalid_argument("CommUnstructured: Invalid location");
  }


  int CommUnstructured::getLocalIdx(int loc, int block) const {
    if (block == -1) block = getBlock(loc);
    if (commRank == block) {
      return loc - iblocks[block];
    } else {
      return -1;
    }
  }


  //===== Internal functions =====//


  //===== Constructors, etc =====//
  CommUnstructured::CommUnstructured()
    : iblock(0), nblocks({}), iblocks({}), nrecv({}), irecv({}), recv_lists({}), send_lists({})
  {}


  CommUnstructured::CommUnstructured(const CommUnstructured& other)
    : Communicator(other), iblock(other.iblock), nblocks(other.nblocks), iblocks(other.iblocks), nrecv(other.nrecv), irecv(other.irecv), recv_lists(other.recv_lists), send_lists(other.send_lists)
  {
    if (usesThisProc) makeMPITypes();
  }


  CommUnstructured& CommUnstructured::operator=(const CommUnstructured& other) {
    Communicator::operator=(other);
    iblock = other.iblock;
    nblocks = other.nblocks;
    iblocks = other.iblocks;
    nrecv = other.nrecv;
    irecv = other.irecv;
    recv_lists = other.recv_lists;
    send_lists = other.send_lists;
    if (usesThisProc) makeMPITypes();
    return *this;
  }

  std::unique_ptr<Communicator> CommUnstructured::clone() const {
    return std::make_unique<CommUnstructured>(static_cast<const CommUnstructured&>(*this));
  }


  //===== Setup =====//
  void CommUnstructured::setup(Potential& pot, size_t ndof, vector<int> ranks) {
    defaultSetup(pot, ndof, ranks);
    if (commSize <= 1) return;

    if (pot.distributed) {
      print("Warning: Attempting to parallelise a potential that has already been distributed.");
      return;
    }

    // Define the sizes of the blocks for each processor
    iblocks = vector<int>(commSize);
    nblocks = vector<int>(commSize, ndof/commSize);
    for (int iProc=0; iProc<commSize-1; iProc++) {
      if (iProc < (int)ndof % commSize) nblocks[iProc]++;
      iblocks[iProc+1] = iblocks[iProc] + nblocks[iProc];
    }

    // Distribute the elements across the processors
    // For each element's dofs get the block that contains it and if it is this block
    vector2d<int> blocks;
    vector2d<bool> in_block;
    getElementBlocks(pot, blocks, in_block);
    setCommLists(pot, blocks, in_block); // Assign send_lists and recv_lists
    setRecvSizes(); // Assign nrecv and irecv
    if (!checkWellDistributed(ndof)) {
      // Move all onto proc 0 to avoid issues arising from nproc == ndof
      nblocks = vector<int>(commSize, 0);
      nblocks[0] = ndof;
      iblocks = vector<int>(commSize, ndof);
      iblocks[0] = 0;
      getElementBlocks(pot, blocks, in_block);
      setCommLists(pot, blocks, in_block);
      setRecvSizes();
    }
    // Update the elements and constraints
    distributeElements(pot, blocks, in_block);

    this->nblock = nblocks[commRank];
    this->nproc = nblock + vec::sum(nrecv);
    this->iblock = iblocks[commRank];

    // Assign the send receive and gather MPI types
    makeMPITypes();
  }


  void CommUnstructured::getElementBlocks(const Potential& pot, vector2d<int>& blocks, vector2d<bool>& in_block) {
    int nElements = pot.elements.size();
    int nConstraints = pot.constraints.size();
    int nTot = nElements + nConstraints;
    blocks = vector2d<int>(nTot);
    in_block = vector2d<bool>(nTot);
    for (int ie=0; ie<nTot; ie++) {
      auto e_idof = (ie<nElements) ? pot.elements[ie].idof : pot.constraints[ie-nElements].idof;
      int e_ndof = e_idof.size();
      blocks[ie] = vector<int>(e_ndof);
      in_block[ie] = vector<bool>(e_ndof);
      for (int i=0; i<e_ndof; i++) {
        blocks[ie][i] = getBlock(e_idof[i]);
        in_block[ie][i] = (blocks[ie][i] == commRank);
      }
    }
  }


  // Populate lists of indicies being sent to and received from each proc
  void CommUnstructured::setCommLists(const Potential& pot, const vector2d<int>& blocks, const vector2d<bool>& in_block) {
    int nElements = pot.elements.size();
    int nConstraints = pot.constraints.size();
    int nTot = nElements + nConstraints;
    recv_lists = vector2d<int>(commSize);
    send_lists = vector2d<int>(commSize);
    for (int ie=0; ie<nTot; ie++) {
      if (vec::all(in_block[ie]) || !vec::any(in_block[ie])) continue;
      auto e_idof = (ie<nElements) ? pot.elements[ie].idof : pot.constraints[ie-nElements].idof;
      int e_ndof = e_idof.size();
      for (int i=0; i<e_ndof; i++) {
        if (in_block[ie][i]) continue;
        vec::insert_unique(recv_lists[blocks[ie][i]], e_idof[i]);
        for (int j=0; j<e_ndof; j++) { // TODO: Improve this part, it will attempt to add the same values multiple times
          if (in_block[ie][j]) {
            vec::insert_unique(send_lists[blocks[ie][i]], e_idof[j] - iblocks[commRank]);
          }
        }
      }
    }
  }


  // Assign the number of coordinates to receive from each proc and the starting index
  void CommUnstructured::setRecvSizes() {
    nrecv = vector<int>(commSize);
    irecv = vector<int>(commSize);
    irecv[0] = nblocks[commRank];
    for (int i=0; i<commSize; i++) {
      nrecv[i] = recv_lists[i].size();
      if (i < commSize-1) irecv[i+1] = irecv[i] + nrecv[i];
    }
  }


  // Get the processor number for each element
  vector<int> CommUnstructured::assignElements(int nElements, const vector2d<int>& blocks) {
    vector<int> el_proc(nElements);
    // Give to the proc with the most DoF contained in the element
    for (int ie=0; ie<nElements; ie++) {
      std::set<int> uniqueBlocks(blocks[ie].begin(), blocks[ie].end());
      int mostDof = 0;
      for (int block : uniqueBlocks) {
        int nDof = 0;
        for (int ib : blocks[ie]) {
          if (ib == block) nDof++;
        }
        if (nDof > mostDof) {
          mostDof = nDof;
          el_proc[ie] = block;
        }
      }
    }
    // Give to whichever proc has the fewest elements
    // vector<int> proc_ne(commSize);
    // for (size_t ie=0; ie<nElements; ie++) {
    //   int fewest_elements = std::numeric_limits<int>::max();
    //   std::set<int> uniqueBlocks(blocks[ie].begin(), blocks[ie].end());
    //   for (int i : uniqueBlocks) {
    //     if (proc_ne[i] < fewest_elements) {
    //       el_proc[ie] = i;
    //       fewest_elements = proc_ne[i];
    //     }
    //   }
    //   proc_ne[el_proc[ie]] ++;
    // }
    return el_proc;
  }


  void CommUnstructured::distributeElements(Potential& pot, const vector2d<int>& blocks, const vector2d<bool>& in_block) {
    vector<int> el_proc = assignElements(pot.elements.size(), blocks);

    int nElements = pot.elements.size();
    int nConstraints = pot.constraints.size();
    int nTot = nElements + nConstraints;
    vector<Potential::Element> elements_tmp;
    vector<Potential::Constraint> constraints_tmp;

    for (int ie=0; ie<nTot; ie++) {
      if (! vec::any(in_block[ie])) continue;
      // Update element.idof with local index
      auto &idof = (ie<nElements) ? pot.elements[ie].idof : pot.constraints[ie-nElements].idof;
      for (int i=0; i<(int)idof.size(); i++) {
        if (in_block[ie][i]) {
          idof[i] = idof[i] - iblocks[commRank];
        } else {
          int block = blocks[ie][i];
          for (int j=0; j<nrecv[block]; j++) {
            if (idof[i] == recv_lists[block][j]) {
              idof[i] = irecv[block] + j;
              break;
            }
          }
        }
      }
      // Assign to the local list of elements or the halo
      if (ie >= nElements) {
        constraints_tmp.push_back(pot.constraints[ie-nElements]);
      } else if (el_proc[ie] == commRank) {
        elements_tmp.push_back(pot.elements[ie]);
      } else {
        pot.elements_halo.push_back(pot.elements[ie]);
      }
    }
    pot.constraints = constraints_tmp;
    pot.elements = elements_tmp;
    pot.distributed = true;
  }


  bool CommUnstructured::checkWellDistributed(int ndof) {
    if (commSize <= 1) return true;
    int nproc = irecv[commSize-1] + nrecv[commSize-1];
    bool wellDistributed = (nproc < ndof);
    #ifdef PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &wellDistributed, 1, MPI_C_BOOL, MPI_LAND, comm);
    #endif

    static bool warnBadDistr = true;
    if (!wellDistributed && warnBadDistr) {
      print("Warning: The state coordinates have not been effectively distributed. Reconsider if MPI is needed.");
      warnBadDistr = false;
    }

    return wellDistributed;
  }


  // Make the MPI datatypes to send and receive from each proc
  void CommUnstructured::makeMPITypes() {
    if (commSize <= 1) return;

    // Send and receive types
    // Send
    for (int i=0; i<commSize; i++) {
      if (send_lists[i].empty()) continue;
      vector<int> blocklens;
      vector<int> disps = {send_lists[i][0]};
      int previous = send_lists[i][0];
      for (int current : send_lists[i]) {
        if (current-previous > 1) {
          disps.push_back(current);
          blocklens.push_back(previous+1 - *std::prev(disps.end(), 2));
        }
        previous = current;
      }
      blocklens.push_back(send_lists[i].back() + 1 - disps.back());
      #ifdef PARALLEL
      MPI_Datatype* newSendType = new MPI_Datatype;
      //MPI_Type_create_indexed_block(send_lists[i].size(), 1, &send_lists[i][0], MPI_DOUBLE, &sendType);
      MPI_Type_indexed(blocklens.size(), &blocklens[0], &disps[0], MPI_DOUBLE, newSendType);
      MPI_Type_commit(newSendType);
      haloTypes.push_back({i, 0, std::shared_ptr<MPI_Datatype>(newSendType, mpiTypeDeleter)});
      #endif
    }
    // Receive
    for (int i=0; i<commSize; i++) {
      if (nrecv[i] == 0) continue;
      #ifdef PARALLEL
      MPI_Datatype* newRecvType = new MPI_Datatype;
      MPI_Type_indexed(1, &nrecv[i], &irecv[i], MPI_DOUBLE, newRecvType);
      MPI_Type_commit(newRecvType);
      edgeTypes.push_back({i, 0, std::shared_ptr<MPI_Datatype>(newRecvType, mpiTypeDeleter)});
      #endif
    }

    // Types for gather
    #ifdef PARALLEL
    // Local block type for sending
    MPI_Datatype* newBlockType = new MPI_Datatype;
    MPI_Type_contiguous(nblocks[commRank], MPI_DOUBLE, newBlockType);
    MPI_Type_commit(newBlockType);
    blockType = std::shared_ptr<MPI_Datatype>(newBlockType, mpiTypeDeleter);
    // Gather type for receiving
    this->nGather = nblocks;
    this->iGather = iblocks;
    MPI_Datatype* newGatherType = new MPI_Datatype;
    MPI_Type_contiguous(1, MPI_DOUBLE, newGatherType);
    MPI_Type_commit(newGatherType);
    gatherType = std::shared_ptr<MPI_Datatype>(newGatherType, mpiTypeDeleter);
    #endif

    mpiTypesCommitted = true;
  }

}
