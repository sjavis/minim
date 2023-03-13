#include "Communicator.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <set>
#include <numeric>
#include <limits>
#include <stdexcept>
#include "utils/vec.h"
#include "utils/mpi.h"
#include "utils/print.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;

  bool warnBadDistr = true;

  class Communicator::Priv {
    public:
      int commRank;
      int commSize;
      vector<int> nblocks; // Size of each block
      vector<int> iblocks; // Global index for the start of each block
      vector<int> nrecv;  // Number of halo coordinates to recieve from each proc
      vector<int> irecv;  // Starting indicies for each proc in halo
      vector<bool> send;  // States if data is to be sent to each proc
      vector2d<int> recv_lists; // List of indicies to recieve from each proc
      vector2d<int> send_lists; // List of block indicies to send to each other proc
  #ifdef PARALLEL
      MPI_Comm comm;
      vector<MPI_Datatype> sendtype; // MPI derived datatype to send to each proc
  #endif


      Priv() : commRank(mpi.rank), commSize(mpi.size) {};


      int getBlock(int loc) {
        for (int i=commSize-1; i>=0; i--) {
          if (loc >= iblocks[i]) return i;
        }
        throw std::invalid_argument("Invalid location");
      }


      // Define the MPI communicator, size, and processor rank
      void setComm(vector<int> ranks) {
#ifdef PARALLEL
        if (ranks.empty()) {
          MPI_Comm_dup(MPI_COMM_WORLD, &comm);
        } else {
          // Make a new communicator for processors with mpi.rank in ranks
          MPI_Group world_group, group;
          MPI_Comm_group(MPI_COMM_WORLD, &world_group);
          MPI_Group_incl(world_group, ranks.size(), &ranks[0], &group);
          int tag = vec::sum(vec::pow(2, ranks)); // Make tag unique for ranks to prevent clashes
          MPI_Comm_create_group(MPI_COMM_WORLD, group, tag, &comm);
        }
        if (comm != MPI_COMM_NULL) {
          MPI_Comm_rank(comm, &commRank);
          MPI_Comm_size(comm, &commSize);
        } else { // This processor is not used
          commRank = -1;
          commSize = -1;
        }
#endif
      }


      void setBlockSizes(int ndof) {
        iblocks = vector<int>(commSize);
        nblocks = vector<int>(commSize, ndof/commSize);
        for (int i=0; i<commSize-1; i++) {
          if (i < (int)ndof % commSize) nblocks[i]++;
          iblocks[i+1] = iblocks[i] + nblocks[i];
        }
      }


      // For each element's dofs get the block that contains it and if it is this block
      void getElementBlocks(const vector<Potential::Element>& elements, vector2d<int>& blocks, vector2d<bool>& in_block) {
        int nelements = elements.size();
        blocks = vector2d<int>(nelements);
        in_block = vector2d<bool>(nelements);
        for (int ie=0; ie<nelements; ie++) {
          auto e = elements[ie];
          int e_ndof = e.idof.size();
          blocks[ie] = vector<int>(e_ndof);
          in_block[ie] = vector<bool>(e_ndof);
          for (int i=0; i<e_ndof; i++) {
            blocks[ie][i] = getBlock(e.idof[i]);
            in_block[ie][i] = (blocks[ie][i] == commRank);
          }
        }
      }

      // Populate lists of indicies being sent to and received from each proc
      void setCommLists(const vector<Potential::Element>& elements, const vector2d<int>& blocks,
                        const vector2d<bool>& in_block)
      {
        recv_lists = vector2d<int>(commSize);
        send_lists = vector2d<int>(commSize);
        for (size_t ie=0; ie<elements.size(); ie++) {
          if (vec::all(in_block[ie]) || !vec::any(in_block[ie])) continue;
          auto e = elements[ie];
          int e_ndof = e.idof.size();
          for (int i=0; i<e_ndof; i++) {
            if (in_block[ie][i]) continue;
            vec::insert_unique(recv_lists[blocks[ie][i]], e.idof[i]);
            for (int j=0; j<e_ndof; j++) { // TODO: Improve this part, it will attempt to add the same values multiple times
              if (in_block[ie][j]) {
                vec::insert_unique(send_lists[blocks[ie][i]], e.idof[j] - iblocks[commRank]);
              }
            }
          }
        }
      }

      // Assign the number of coordinates to receive from each proc and the starting index
      void setRecvSizes() {
        nrecv = vector<int>(commSize);
        irecv = vector<int>(commSize);
        irecv[0] = nblocks[commRank];
        for (int i=0; i<commSize; i++) {
          nrecv[i] = recv_lists[i].size();
          if (i < commSize-1) irecv[i+1] = irecv[i] + nrecv[i];
        }
      }

      // Get the processor number for each element
      vector<int> assignElements(const vector<Potential::Element>& elements, const vector2d<int>& blocks) {
        vector<int> el_proc(elements.size());
        // Give to the proc with the most DoF contained in the element
        for (size_t ie=0; ie<elements.size(); ie++) {
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
        // for (size_t ie=0; ie<elements.size(); ie++) {
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

      // Equally distribute the elements among the processors
      void distributeElements(Potential& pot, const vector2d<int>& blocks, const vector2d<bool>& in_block) {
        vector<int> el_proc = assignElements(pot.elements, blocks);
        vector<Potential::Element> elements_tmp;
        for (size_t ie=0; ie<pot.elements.size(); ie++) {
          if (! vec::any(in_block[ie])) continue;
          // Update element.idof with local index
          auto e = pot.elements[ie];
          int idof_size = e.idof.size();
          for (int i=0; i<idof_size; i++) {
            if (in_block[ie][i]) {
              e.idof[i] = e.idof[i] - iblocks[commRank];
            } else {
              int block = blocks[ie][i];
              for (int j=0; j<nrecv[block]; j++) {
                if (e.idof[i] == recv_lists[block][j]) {
                  e.idof[i] = irecv[block] + j;
                  break;
                }
              }
            }
          }
          // Assign to the local list of elements or the halo
          if (el_proc[ie] == commRank) {
            elements_tmp.push_back(e);
          } else {
            pot.elements_halo.push_back(e);
          }
        }
        pot.elements = elements_tmp;
        pot.distributed = true;
      }

      // Make the MPI datatypes to send to each proc
      void setSendType() {
        if (commSize <= 1) return;
        send = vector<bool>(commSize);
#ifdef PARALLEL
        sendtype = vector<MPI_Datatype>(commSize);
#endif
        for (int i=0; i<commSize; i++) {
          if (send_lists[i].empty()) continue;
          send[i] = true;
          vector<int> blocklens;
          vector<int> disps = {send_lists[i][0]};
          int previous = send_lists[i][0];
          for (int current : send_lists[i]) {
            if (current-previous > 1) {
              disps.push_back(current);
              blocklens.push_back(disps.back() - *std::prev(disps.end(), 2));
            }
            previous = current;
          }
          blocklens.push_back(send_lists[i].back() + 1 - disps.back());
#ifdef PARALLEL
          //MPI_Type_create_indexed_block(send_lists[i].size(), 1, &send_lists[i][0], MPI_DOUBLE, &priv->sendtype[i]);
          MPI_Type_indexed(blocklens.size(), &blocklens[0], &disps[0], MPI_DOUBLE, &sendtype[i]);
          MPI_Type_commit(&sendtype[i]);
#endif
        }
      }


      void checkWellDistributed(int ndof, Potential& pot, vector2d<int>& blocks, vector2d<bool>& in_block) {
        if (commSize <= 1) return;
        int nproc = irecv[commSize-1] + nrecv[commSize-1];
        bool notDistributed = (nproc == ndof);
#ifdef PARALLEL
        MPI_Allreduce(MPI_IN_PLACE, &notDistributed, 1, MPI_C_BOOL, MPI_LOR, comm);
#endif
        if (notDistributed) {
          if (warnBadDistr) print("Warning: The state coordinates have not been effectively distributed. Reconsider if MPI is needed.");
          warnBadDistr = false;
          // Move all onto proc 0 to avoid issues arising from nproc == ndof
          nblocks = vector<int>(commSize, 0);
          nblocks[0] = ndof;
          iblocks = vector<int>(commSize, ndof);
          iblocks[0] = 0;
          getElementBlocks(pot.elements, blocks, in_block);
          setCommLists(pot.elements, blocks, in_block);
          setRecvSizes();
        }
      }


      void setup(Potential& pot, int ndof, vector<int> ranks) {
        if (pot.distributed) {
          print("Warning: Attempting to parallelise a potential that has already been distributed.");
          return;
        }

        setComm(ranks);
        if (commRank<0) return;

        setBlockSizes(ndof);

        vector2d<int> blocks;
        vector2d<bool> in_block;
        getElementBlocks(pot.elements, blocks, in_block);

        // Identify coordinates to send and receive based upon the list of energy elements
        setCommLists(pot.elements, blocks, in_block);
        setRecvSizes();

        checkWellDistributed(ndof, pot, blocks, in_block);

        distributeElements(pot, blocks, in_block);
        setSendType();
      }
  };


  Communicator::Communicator()
    : ndof(0), nproc(0), nblock(0), iblock(0), ranks({0}), p(std::unique_ptr<Priv>(new Priv))
  {}


  Communicator::Communicator(const Communicator& comm)
    : ndof(comm.ndof), nproc(comm.nproc), nblock(comm.nblock), iblock(comm.iblock),
      usesThisProc(comm.usesThisProc), ranks(comm.ranks), p(std::make_unique<Priv>(*comm.p))
  {
    if (usesThisProc) p->setSendType();
  }


  Communicator& Communicator::operator=(const Communicator& comm) {
    ndof = comm.ndof;
    nproc = comm.nproc;
    nblock = comm.nblock;
    iblock = comm.iblock;
    usesThisProc = comm.usesThisProc;
    ranks = comm.ranks;
    p = std::make_unique<Priv>(*comm.p);
    if (usesThisProc) p->setSendType();
    return *this;
  }


  Communicator::~Communicator() {
  #ifdef PARALLEL
    if (p->commSize <= 1) return;
    // Free any committed MPI datatypes
    for (int i=0; i<p->commSize; i++) {
      if (p->send[i]) {
        MPI_Type_free(&p->sendtype[i]);
      }
    }
  #endif
  }


  void Communicator::setup(Potential& pot, size_t ndof, vector<int> ranks) {
    // Serial parameters
    this->ndof = ndof;
    this->nproc = ndof;
    this->nblock = ndof;
    if (mpi.size == 1) return;

    // Parallel parameters
    p->setup(pot, ndof, ranks);
    usesThisProc = (p->commRank >= 0);
    this->ranks = ranks;
    if (ranks.empty()) {
      this->ranks = vector<int>(mpi.size);
      std::iota(this->ranks.begin(), this->ranks.end(), 0);
    }
    if (usesThisProc) {
      nblock = p->nblocks[p->commRank];
      nproc = p->irecv[p->commSize-1] + p->nrecv[p->commSize-1];
      iblock = p->iblocks[p->commRank];
    } else {
      nblock = -1;
      nproc = -1;
      iblock = -1;
    }
  }


  int Communicator::rank() const {
    return p->commRank;
  }

  int Communicator::size() const {
    return p->commSize;
  }


  template <typename T>
  vector<T> Communicator::assignBlock(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (p->commSize == 1) return in;

    vector<T> out(nblock);
    int i0 = (in.size() == ndof) ? iblock : 0;
    for (size_t i=0; i<nblock; i++) {
      out[i] = in[i0+i];
    }
    return out;
  }
  template vector<int> Communicator::assignBlock(const vector<int>&) const;
  template vector<bool> Communicator::assignBlock(const vector<bool>&) const;
  template vector<double> Communicator::assignBlock(const vector<double>&) const;


  template <typename T>
  vector<T> Communicator::assignProc(const vector<T>& in) const {
    if (!usesThisProc) return vector<T>();
    if (in.size() != ndof) throw std::invalid_argument("Input data has incorrect size. All degrees of freedom required.");
    if (p->commSize == 1) return in;

    vector<T> out(nproc);
    // Assign the main blocks
    for (size_t i=0; i<nblock; i++) {
      out[i] = in[iblock+i];
    }
    // Assign the halo regions
    for (int i=0; i<p->commSize; i++) {
      for (int j=0; j<p->nrecv[i]; j++) {
        out[p->irecv[i]+j] = in[p->recv_lists[i][j]];
      }
    }
    return out;
  }
  template vector<int> Communicator::assignProc(const vector<int>&) const;
  template vector<bool> Communicator::assignProc(const vector<bool>&) const;
  template vector<double> Communicator::assignProc(const vector<double>&) const;


  void Communicator::communicate(vector<double>& vector) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    for (int i=0; i<p->commSize; i++) {
      if (i == p->commRank) continue;
      // Send
      if (p->send[i]) {
        int tag = p->commRank*p->commSize + i;
        MPI_Send(&vector[0], 1, p->sendtype[i], i, tag, p->comm);
      }
      // Receive
      if (p->nrecv[i] > 0) {
        int tag = i*p->commSize + p->commRank;
        int irecv = std::accumulate(p->nrecv.begin(), p->nrecv.begin()+i, nblock);
        MPI_Recv(&vector[irecv], p->nrecv[i], MPI_DOUBLE, i, tag, p->comm, nullptr);
      }
    }
  #endif
  }


  double Communicator::get(const vector<double>& vector, int loc) const {
    if (!usesThisProc) return 0;
  #ifdef PARALLEL
    if (p->commSize > 1) {
      int i = p->getBlock(loc);
      double value = vector[loc - p->iblocks[i]];
      MPI_Bcast(&value, 1, MPI_DOUBLE, i, p->comm);
      return value;
    }
  #endif
    return vector[loc];
  }


  double Communicator::sum(double a) const {
    if (!usesThisProc) return 0;
    double result = a;
  #ifdef PARALLEL
    if (p->commSize > 1) MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, p->comm);
  #endif
    return result;
  }

  double Communicator::sum(const vector<double>& a) const {
    if (!usesThisProc) return 0;
    return sum(vec::sum(a));
  }


  double Communicator::dotProduct(const vector<double>& a, const vector<double>& b) const {
    if (!usesThisProc) return 0;
    return sum(std::inner_product(a.begin(), a.begin()+nblock, b.begin(), 0.0));
  }


  vector<double> Communicator::gather(const vector<double>& block, int root) const {
    if (!usesThisProc) return vector<double>();

  #ifdef PARALLEL
    if (p->commSize > 1) {
      vector<double> gathered;
      if (root == -1) {
        gathered = vector<double>(ndof);
        MPI_Allgatherv(&block[0], nblock, MPI_DOUBLE,
                       &gathered[0], &p->nblocks[0], &p->iblocks[0], MPI_DOUBLE,
                       p->comm);
      } else {
        if (p->commRank==root) gathered = vector<double>(ndof);
        MPI_Gatherv(&block[0], nblock, MPI_DOUBLE,
                    &gathered[0], &p->nblocks[0], &p->iblocks[0], MPI_DOUBLE, root,
                    p->comm);
      }
      return gathered;
    }
  #endif
    return block;
  }


  vector<double> Communicator::scatter(const vector<double>& data, int root) const {
    if (!usesThisProc) return vector<double>();
  #ifdef PARALLEL
    if (p->commSize > 1) {
      // Get copy of data on processor (potentially inefficient)
      vector<double> data_copy;
      if (root == -1) {
        data_copy = data;
        // Note: Halo data may not be correct if 'data' is different on each processor
      } else {
        data_copy = (p->commRank==root) ? data : vector<double>(ndof);
        bcast(data_copy, root);
      }
      return assignProc(data_copy);
    }
  #endif
    return data;
  }


  void Communicator::bcast(int& value, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (p->commSize > 1) MPI_Bcast(&value, 1, MPI_INT, root, p->comm);
  #endif
  }


  void Communicator::bcast(double& value, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (p->commSize > 1) MPI_Bcast(&value, 1, MPI_DOUBLE, root, p->comm);
  #endif
  }


  void Communicator::bcast(vector<double>& vector, int root) const {
    if (!usesThisProc) return;
  #ifdef PARALLEL
    if (p->commSize > 1) MPI_Bcast(&vector[0], vector.size(), MPI_DOUBLE, root, p->comm);
  #endif
  }

}
