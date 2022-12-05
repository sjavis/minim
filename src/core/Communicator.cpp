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

  typedef std::vector<double> Vector;

  bool warnBadDistr = true;

  class Communicator::Priv {
    public:
      Priv() : commRank(mpi.rank), commSize(mpi.size) {};

      int commRank;
      int commSize;
      std::vector<int> nblocks; // Size of each block
      std::vector<int> iblocks; // Global index for the start of each block
      std::vector<int> nrecv;  // Number of halo coordinates to recieve from each proc
      std::vector<int> irecv;  // Starting indicies for each proc in halo
      std::vector<bool> send;  // States if data is to be sent to each proc
      std::vector<std::vector<int>> recv_lists; // List of indicies to recieve from each proc
  #ifdef PARALLEL
      MPI_Comm comm;
      std::vector<MPI_Datatype> sendtype; // MPI derived datatype to send to each proc
  #endif


      int getBlock(int loc) {
        for (int i=commSize-1; i>=0; i--) {
          if (loc >= iblocks[i]) return i;
        }
        throw std::invalid_argument("Invalid location");
      }


      void setComm(bool useProc) {
        // Get the rank, size, and new MPI communicator
#ifdef PARALLEL
        MPI_Comm_split(MPI_COMM_WORLD, (useProc)?0:MPI_UNDEFINED, mpi.rank, &comm);
        if (useProc) {
          MPI_Comm_rank(comm, &commRank);
          MPI_Comm_size(comm, &commSize);
        } else {
          commRank = -1;
          commSize = -1;
        }
#endif
      }


      void setBlockSizes(int ndof) {
        iblocks = std::vector<int>(commSize);
        nblocks = std::vector<int>(commSize, ndof/commSize);
        for (int i=0; i<commSize-1; i++) {
          if (i < (int)ndof % commSize) nblocks[i]++;
          iblocks[i+1] = iblocks[i] + nblocks[i];
        }
      }


      // For each element's dofs get the block that contains it and if it is this block
      void getElementBlocks(const std::vector<Potential::Element>& elements,
                            std::vector<std::vector<int>>& blocks,
                            std::vector<std::vector<bool>>& in_block)
      {
        int nelements = elements.size();
        blocks = std::vector<std::vector<int>>(nelements);
        in_block = std::vector<std::vector<bool>>(nelements);
        for (int ie=0; ie<nelements; ie++) {
          auto e = elements[ie];
          int e_ndof = e.idof.size();
          blocks[ie] = std::vector<int>(e_ndof);
          in_block[ie] = std::vector<bool>(e_ndof);
          for (int i=0; i<e_ndof; i++) {
            blocks[ie][i] = getBlock(e.idof[i]);
            in_block[ie][i] = (blocks[ie][i] == commRank);
          }
        }
      }

      // Populate lists of indicies being sent to and received from each proc
      void setCommLists(const std::vector<Potential::Element>& elements,
                        const std::vector<std::vector<int>>& blocks,
                        const std::vector<std::vector<bool>>& in_block,
                        std::vector<std::vector<int>>& send_lists)
      {
        recv_lists = std::vector<std::vector<int>>(commSize);
        send_lists = std::vector<std::vector<int>>(commSize);
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
        nrecv = std::vector<int>(commSize);
        irecv = std::vector<int>(commSize);
        irecv[0] = nblocks[commRank];
        for (int i=0; i<commSize; i++) {
          nrecv[i] = recv_lists[i].size();
          if (i < commSize-1) irecv[i+1] = irecv[i] + nrecv[i];
        }
      }

      // Get the processor number for each element
      std::vector<int> assignElements(const std::vector<Potential::Element>& elements,
                                      const std::vector<std::vector<int>>& blocks)
      {
        std::vector<int> el_proc(elements.size());
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
        // std::vector<int> proc_ne(commSize);
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
      void distributeElements(Potential& pot,
                              const std::vector<std::vector<int>>& blocks,
                              const std::vector<std::vector<bool>>& in_block)
      {
        std::vector<int> el_proc = assignElements(pot.elements, blocks);
        std::vector<Potential::Element> elements_tmp;
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
      }

      // Make the MPI datatypes to send to each proc
      void setSendType(const std::vector<std::vector<int>>& send_lists) {
        send = std::vector<bool>(commSize);
        sendtype = std::vector<MPI_Datatype>(commSize);
        for (int i=0; i<commSize; i++) {
          if (send_lists[i].empty()) continue;
          send[i] = true;
          std::vector<int> blocklens;
          std::vector<int> disps = {send_lists[i][0]};
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


      void checkWellDistributed(int ndof, Potential& pot,
                                std::vector<std::vector<int>>& blocks,
                                std::vector<std::vector<bool>>& in_block,
                                std::vector<std::vector<int>>& send_lists) {
        if (commSize == 1) return;
        int nproc = irecv[commSize-1] + nrecv[commSize-1];
        bool notDistributed = (nproc == ndof);
#ifdef PARALLEL
        MPI_Allreduce(&notDistributed, &notDistributed, 1, MPI_C_BOOL, MPI_LAND, comm);
#endif
        if (notDistributed) {
          if (warnBadDistr) print("Warning: The state coordinates have not been effectively distributed. Reconsider if MPI is needed.");
          warnBadDistr = false;
          // Move all onto proc 0 to avoid issues arising from nproc == ndof
          nblocks = std::vector<int>(commSize, 0);
          nblocks[0] = ndof;
          iblocks = std::vector<int>(commSize, ndof);
          iblocks[0] = 0;
          getElementBlocks(pot.elements, blocks, in_block);
          setCommLists(pot.elements, blocks, in_block, send_lists);
          setRecvSizes();
        }
      }


      void setup(Potential& pot, int ndof, bool useProc) {
        setComm(useProc);
        if (!useProc) return;

        setBlockSizes(ndof);

        std::vector<std::vector<int>> blocks;
        std::vector<std::vector<bool>> in_block;
        getElementBlocks(pot.elements, blocks, in_block);

        // Identify coordinates to send and receive based upon the list of energy elements
        std::vector<std::vector<int>> send_lists; // Block indicies to send to each other block
        setCommLists(pot.elements, blocks, in_block, send_lists);
        setRecvSizes();

        checkWellDistributed(ndof, pot, blocks, in_block, send_lists);

        distributeElements(pot, blocks, in_block);
        setSendType(send_lists);
      }
  };


  Communicator::Communicator(Potential& pot, size_t ndof, std::vector<int> ranks)
    : ndof(ndof), nproc(ndof), nblock(ndof), iblock(0), p(std::unique_ptr<Priv>(new Priv))
  {
    if (mpi.size == 1) return;
    if (!ranks.empty()) usesThisProc = std::count(ranks.begin(), ranks.end(), mpi.rank);

    p->setup(pot, ndof, usesThisProc);
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


  Communicator::Communicator(const Communicator& comm)
    : ndof(comm.ndof), nproc(comm.nproc), nblock(comm.nblock),
      iblock(comm.iblock), p(std::make_unique<Priv>(*comm.p))
  {}


  Communicator& Communicator::operator=(const Communicator& comm) {
    ndof = comm.ndof;
    nproc = comm.nproc;
    nblock = comm.nblock;
    iblock = comm.iblock;
    p = std::make_unique<Priv>(*comm.p);
    return *this;
  }


  Communicator::~Communicator() {
  #ifdef PARALLEL
    // Free any committed MPI datatypes
    for (int i=0; i<p->commSize; i++) {
      if (p->send[i]) {
        MPI_Type_free(&p->sendtype[i]);
      }
    }
  #endif
  }


  int Communicator::rank() const {
    return p->commRank;
  }

  int Communicator::size() const {
    return p->commSize;
  }


  Vector Communicator::assignBlock(const Vector& in) const {
    Vector out = Vector(nblock);
    int i0 = (in.size() == ndof) ? iblock : 0;
    for (size_t i=0; i<nblock; i++) {
      out[i] = in[i0+i];
    }
    return out;
  }


  Vector Communicator::assignProc(const Vector& in) const {
    if (in.size() != ndof) throw std::invalid_argument("Input data has incorrect size. All degrees of freedom required.");
    Vector out = Vector(nproc);
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


  void Communicator::communicate(Vector& vector) const {
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


  double Communicator::get(const Vector& vector, int loc) const {
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
    double result = a;
  #ifdef PARALLEL
    MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, p->comm);
  #endif
    return result;
  }

  double Communicator::sum(const Vector& a) const {
    return sum(vec::sum(a));
  }


  double Communicator::dotProduct(const Vector& a, const Vector& b) const {
    return sum(std::inner_product(a.begin(), a.begin()+nblock, b.begin(), 0.0));
  }


  Vector Communicator::gather(const Vector& block, int root) const {

  #ifdef PARALLEL
    if (p->commSize > 1) {
      Vector gathered;
      if (root == -1) {
        gathered = Vector(ndof);
        MPI_Allgatherv(&block[0], nblock, MPI_DOUBLE,
                       &gathered[0], &p->nblocks[0], &p->iblocks[0], MPI_DOUBLE,
                       p->comm);
      } else {
        if (p->commRank==root) gathered = Vector(ndof);
        MPI_Gatherv(&block[0], nblock, MPI_DOUBLE,
                    &gathered[0], &p->nblocks[0], &p->iblocks[0], MPI_DOUBLE, root,
                    p->comm);
      }
      return gathered;
    }
  #endif
    return block;
  }


  Vector Communicator::scatter(const Vector& data, int root) const {
  #ifdef PARALLEL
    if (p->commSize > 1) {
      // Get copy of data on processor (potentially inefficient)
      Vector data_copy;
      if (root == -1) {
        data_copy = data;
        // Note: Halo data may not be correct if 'data' is different on each processor
      } else {
        data_copy = (p->commRank==root) ? data : Vector(ndof);
        bcast(data_copy, root);
      }
      return assignProc(data_copy);
    }
  #endif
    return data;
  }


  void Communicator::bcast(int& value, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&value, 1, MPI_INT, root, p->comm);
  #endif
  }


  void Communicator::bcast(double& value, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&value, 1, MPI_DOUBLE, root, p->comm);
  #endif
  }


  void Communicator::bcast(Vector& vector, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&vector[0], vector.size(), MPI_DOUBLE, root, p->comm);
  #endif
  }

}
