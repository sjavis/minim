#include "Communicator.h"

#include <numeric>
#include <stdexcept>
#include "Potential.h"
#include "utils/vec.h"
#include "utils/mpi.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  int Communicator::rank() const {
    return commRank;
  }

  int Communicator::size() const {
    return commSize;
  }


  //===== Access data =====//
  double Communicator::get(const vector<double>& vector, int loc) const {
    #ifdef PARALLEL
    if (!usesThisProc) return 0;
    if (commSize > 1) {
      int root = getBlock(loc);
      double value;
      if (commRank == root) value = vector[getLocalIdx(loc, root)];
      MPI_Bcast(&value, 1, MPI_DOUBLE, root, comm);
      return value;
    }
    #endif
    return vector[loc];
  }


  //===== Communicate =====//
  void Communicator::communicate(vector<double>& vector) const {
    if (!usesThisProc || commSize==1) return;
    #ifdef PARALLEL
    int nRequest = haloTypes.size() + edgeTypes.size();
    MPI_Request requests[nRequest];
    int iRequest = 0;
    for (const auto& sendType : haloTypes) {
      MPI_Isend(&vector[0], 1, *sendType.type, sendType.rank, sendType.tag, comm, &requests[iRequest++]);
    }
    for (const auto& recvType : edgeTypes) {
      MPI_Irecv(&vector[0], 1, *recvType.type, recvType.rank, recvType.tag, comm, &requests[iRequest++]);
    }
    MPI_Waitall(nRequest, requests, MPI_STATUSES_IGNORE);
    #endif
  }

  void Communicator::communicateAccumulate(vector<double>& vector) const {
    if (!usesThisProc || commSize==1) return;
    #ifdef PARALLEL
    MPI_Win win;
    MPI_Win_create(vector.data(), vector.size()*sizeof(double), sizeof(double), MPI_INFO_NULL, comm, &win);

    // MPI_Win_fence(MPI_MODE_NOPRECEDE, win);
    for (int iDir=0; iDir<(int)edgeTypes.size(); iDir++) {
      const CommunicateObj& sendType = edgeTypes[iDir];
      const CommunicateObj& recvType = haloTypes[iDir];
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, sendType.rank, 0, win);
      MPI_Accumulate(vector.data(), 1, *sendType.type, sendType.rank, 0, 1, *recvType.type, MPI_SUM, win);
      MPI_Win_unlock(sendType.rank, win);
    }
    // MPI_Win_fence(MPI_MODE_NOSUCCEED, win);
    MPI_Barrier(comm);
    MPI_Win_free(&win);
    #endif
  }


  vector<double> Communicator::gather(const vector<double>& block, int root) const {
    if (!usesThisProc) return vector<double>();
    if (commSize == 1) return block;

  #ifdef PARALLEL
    vector<double> gathered;
    if (root == -1) {
      gathered = vector<double>(ndof);
      MPI_Allgatherv(&block[0], 1, *blockType,
                     &gathered[0], &nGather[0], &iGather[0], *gatherType,
                     comm);
    } else {
      if (commRank==root) gathered = vector<double>(ndof);
      MPI_Gatherv(&block[0], 1, *blockType,
                  &gathered[0], &nGather[0], &iGather[0], *gatherType, root,
                  comm);
    }
    return gathered;

  #else
    return block;
  #endif
  }


  vector<double> Communicator::scatter(const vector<double>& data, int root) const {
    if (!usesThisProc) return vector<double>();
  #ifdef PARALLEL
    if (commSize > 1) {
      if (root == -1) {
        // Note: Halo data may not be correct if 'data' is different on each processor
        return assignProc(data);
      } else {
        // Get copy of data on processor (potentially inefficient)
        vector<double> data_copy = (commRank==root) ? data : vector<double>(ndof);
        bcast(data_copy, root);
        return assignProc(data_copy);
      }
    }
  #endif
    return data;
  }


  void Communicator::bcast(int& value, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (commSize > 1) MPI_Bcast(&value, 1, MPI_INT, root, comm);
  #endif
  }


  void Communicator::bcast(double& value, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (commSize > 1) MPI_Bcast(&value, 1, MPI_DOUBLE, root, comm);
  #endif
  }


  void Communicator::bcast(vector<double>& vector, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (commSize > 1) MPI_Bcast(&vector[0], vector.size(), MPI_DOUBLE, root, comm);
  #endif
  }


  //===== MPI reduction functions =====//
  double Communicator::sum(double a) const {
    if (!usesThisProc) return 0;
    double result = a;
  #ifdef PARALLEL
    if (commSize > 1) MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
  #endif
    return result;
  }

  double Communicator::sum(const vector<double>& a) const {
    if (!usesThisProc) return 0;
    return sum(vec::sum(a)); // TODO: Ignore halo?
  }


  double Communicator::norm(const vector<double>& a) const {
    if (!usesThisProc) return 0;
    return sqrt(dotProduct(a, a));
  }


  double Communicator::dotProduct(const vector<double>& a, const vector<double>& b) const {
    if (!usesThisProc) return 0;
    return sum(std::inner_product(a.begin(), a.begin()+nblock, b.begin(), 0.0));
  }


  //===== Internal functions =====//

  #ifdef PARALLEL
  // Custom deleter to free MPI datatypes stored as shared_ptrs
  void Communicator::mpiTypeDeleter(MPI_Datatype* type) {
    if (*type!=MPI_DATATYPE_NULL) MPI_Type_free(type);
    delete type;
  }
  #endif


  void Communicator::defaultSetup(const Potential& pot, size_t ndof, vector<int> ranks) {
    // Default values in case of a serial run
    this->commSize = 1;
    this->commRank = 0;
    this->ndof = ndof;
    this->nproc = ndof;
    this->nblock = ndof;

    // Set the ranks
    if (ranks.empty()) {
      this->ranks = vector<int>(mpi.size);
      std::iota(this->ranks.begin(), this->ranks.end(), 0);
    } else {
      this->ranks = ranks;
    }
    this->usesThisProc = vec::isIn(this->ranks, mpi.rank);

    // Null values if not using this processor
    if (!usesThisProc) {
      this->commRank = -1;
      this->commSize = -1;
      this->nblock = -1;
      this->nproc = -1;
      return;
    }

    if (mpi.size==1 || ranks.size()==1 || pot.isSerial()) return;
    setComm(ranks);
  }


  // Define the MPI communicator, size, and processor rank
  void Communicator::setComm(vector<int> ranks) {
    #ifdef PARALLEL
    // Make a new communicator for processors with mpi.rank in ranks
    if (ranks.empty()) {
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    } else {
      MPI_Group world_group, group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      MPI_Group_incl(world_group, ranks.size(), &ranks[0], &group);
      int tag = rand()%1000000;
      MPI_Bcast(&tag, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Comm_create_group(MPI_COMM_WORLD, group, tag, &comm);
    }
    // Get the rank number and size of the communicator
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_rank(comm, &commRank);
      MPI_Comm_size(comm, &commSize);
    }
    #endif
  }

}
