#include "Communicator.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <numeric>
#include <stdexcept>
#include "utils/vec.h"
#include "utils/mpi.h"

namespace minim {

  typedef std::vector<double> Vector;


  class Communicator::Priv {
    public:
      Priv() :
        nblocks(mpi.size),
        iblocks(mpi.size),
        nrecv(mpi.size),
        irecv(mpi.size),
        send(mpi.size),
        recv_lists(mpi.size)
      {};

      std::vector<int> nblocks; // Size of each block
      std::vector<int> iblocks; // Global index for the start of each block
      std::vector<int> nrecv;  // Number of halo coordinates to recieve from each proc
      std::vector<int> irecv;  // Starting indicies for each proc in halo
      std::vector<bool> send;  // States if data is to be sent to each proc
      std::vector<std::vector<int>> recv_lists; // List of indicies to recieve from each proc
  #ifdef PARALLEL
      // MPI derived datatype to send to each proc
      std::vector<MPI_Datatype> sendtype = std::vector<MPI_Datatype>(mpi.size);
  #endif

      int getBlock(int loc) {
        for (int i=mpi.size-1; i>=0; i--) {
          if (loc >= iblocks[i]) return i;
        }
        throw std::invalid_argument("Invalid location");
      }
  };


  Communicator::Communicator(int ndof, Potential::Args& args)
    : ndof(ndof), nproc(ndof), nblock(ndof), priv(std::unique_ptr<Priv>(new Priv))
  {
    priv->nblocks = std::vector<int>(mpi.size, ndof/mpi.size);
    if (mpi.size == 1) return;

    // Get size and index of each block
    for (int i=0; i<mpi.size-1; i++) {
      if (i < ndof % mpi.size) priv->nblocks[i]++;
      priv->iblocks[i+1] = priv->iblocks[i] + priv->nblocks[i];
    }
    nblock = priv->nblocks[mpi.rank];

    // Identify halo coordinates based upon the list of energy elements
    std::vector<std::vector<int>> send_lists(mpi.size);
    std::vector<std::vector<int>> blocks(args.elements.size());
    std::vector<std::vector<bool>> in_block(args.elements.size());
    int nelements = args.elements.size();
    for (int ie=0; ie<nelements; ie++) {
      auto e = args.elements[ie];
      int e_ndof = e.idof.size();
      blocks[ie] = std::vector<int>(e_ndof);
      in_block[ie] = std::vector<bool>(e_ndof);

      // Get block numbers for each element's dof
      for (int i=0; i<e_ndof; i++) {
        blocks[ie][i] = priv->getBlock(e.idof[i]);
        in_block[ie][i] = (blocks[ie][i] == mpi.rank);
      }

      if (vec::all(in_block[ie]) || !vec::any(in_block[ie])) continue;

      // Populate lists of indicies being sent to and received from each proc
      for (int i=0; i<e_ndof; i++) {
        if (in_block[ie][i]) continue;
        vec::insert_unique(priv->recv_lists[blocks[ie][i]], e.idof[i]);
        for (int j=0; j<e_ndof; j++) {
          if (in_block[ie][j]) {
            vec::insert_unique(send_lists[blocks[ie][i]], e.idof[j] - priv->iblocks[mpi.rank]);
          }
        }
      }
    } //end for elements

    // Store number and starting indicies of coords being recieved from each proc
    priv->irecv[0] = nblock;
    for (int i=0; i<mpi.size; i++) {
      priv->nrecv[i] = priv->recv_lists[i].size();
      if (i < mpi.size-1) priv->irecv[i+1] = priv->irecv[i] + priv->nrecv[i];
    }
    nproc = priv->irecv[mpi.size-1] + priv->nrecv[mpi.size-1];

    std::vector<int> nelements_blocks(mpi.size);
    std::vector<Potential::Args::Element> elements_tmp;
    for (int ie=0; ie<nelements; ie++) {
      // Assign each element to a proc
      int proc = blocks[ie][0];
      int fewest_elements = nelements_blocks[proc];
      for (int i : blocks[ie]) {
        if (nelements_blocks[i] < fewest_elements) {
          proc = i;
          fewest_elements = nelements_blocks[i];
        }
      }
      nelements_blocks[proc] ++;

      // Store the elements for this proc
      if (vec::any(in_block[ie])) {
        auto e = args.elements[ie];
        // Update element.idof with local index
        int idof_size = e.idof.size();
        for (int i=0; i<idof_size; i++) {
          if (in_block[ie][i]) {
            e.idof[i] = e.idof[i] - priv->iblocks[mpi.rank];
          } else {
            int block = blocks[ie][i];
            for (int j=0; j<priv->nrecv[block]; j++) {
              if (e.idof[i] == priv->recv_lists[block][j]) {
                e.idof[i] = priv->irecv[block] + j;
              }
            }
          }
        }
        if (proc == mpi.rank) {
          elements_tmp.push_back(e);
        } else {
          args.elements_halo.push_back(e);
        }
      }
    }
    args.elements = elements_tmp;

    // Make the MPI datatypes to send to each proc
    for (int i=0; i<mpi.size; i++) {
      if (send_lists[i].empty()) continue;
      priv->send[i] = true;
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
      MPI_Type_indexed(blocklens.size(), &blocklens[0], &disps[0], MPI_DOUBLE, &priv->sendtype[i]);
      MPI_Type_commit(&priv->sendtype[i]);
  #endif
    }
  }


  Communicator::Communicator(const Communicator& comm)
    : ndof(comm.ndof), nproc(comm.nproc), nblock(comm.nblock), priv(std::make_unique<Priv>(*comm.priv))
  {}


  Communicator& Communicator::operator=(const Communicator& comm) {
    ndof = comm.ndof;
    nproc = comm.nproc;
    nblock = comm.nblock;
    priv = std::make_unique<Priv>(*comm.priv);
    return *this;
  }


  Communicator::~Communicator() {
  #ifdef PARALLEL
    // Free any committed MPI datatypes
    for (int i=0; i<mpi.size; i++) {
      if (priv->send[i]) {
        MPI_Type_free(&priv->sendtype[i]);
      }
    }
  #endif
  }


  Vector Communicator::assignBlock(const Vector& in) const {
    Vector out = Vector(nproc);
    for (int i=0; i<nblock; i++) {
      out[i] = in[priv->iblocks[mpi.rank]+i];
    }
    return out;
  }


  void Communicator::communicate(Vector& vector) const {
  #ifdef PARALLEL
    for (int i=0; i<mpi.size; i++) {
      if (i == mpi.rank) continue;
      // Send
      if (priv->send[i]) {
        int tag = mpi.rank*mpi.size + i;
        MPI_Send(&vector[0], 1, priv->sendtype[i], i, tag, MPI_COMM_WORLD);
      }
      // Receive
      if (priv->nrecv[i] > 0) {
        int tag = i*mpi.size + mpi.rank;
        int irecv = std::accumulate(priv->nrecv.begin(), priv->nrecv.begin()+i, nblock);
        MPI_Recv(&vector[irecv], priv->nrecv[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD, nullptr);
      }
    }
  #endif
  }


  double Communicator::get(const Vector& vector, int loc) const {
    double value;

    if (mpi.size == 1) {
      value = vector[loc];
    } else {
  #ifdef PARALLEL
      int i = priv->getBlock(loc);
      value = vector[loc-priv->iblocks[i]];
      MPI_Bcast(&value, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
  #endif
    }

    return value;
  }


  double Communicator::dotProduct(const Vector& a, const Vector& b) const {
    double result = std::inner_product(a.begin(), a.begin()+nblock, b.begin(), 0.0);
  #ifdef PARALLEL
    MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
    return result;
  }


  Vector Communicator::gather(const Vector& block, int root) const {
  #ifdef PARALLEL
    Vector gathered;
    if (root == -1) {
      gathered = Vector(ndof);
      MPI_Allgatherv(&block[0], nblock, MPI_DOUBLE,
                     &gathered[0], &priv->nblocks[0], &priv->iblocks[0], MPI_DOUBLE,
                     MPI_COMM_WORLD);
    } else {
      if (mpi.rank==root) gathered = Vector(ndof);
      MPI_Gatherv(&block[0], nblock, MPI_DOUBLE,
                  &gathered[0], &priv->nblocks[0], &priv->iblocks[0], MPI_DOUBLE, root,
                  MPI_COMM_WORLD);
    }
    return gathered;
  #else
    return block;
  #endif
  }


  Vector Communicator::scatter(const Vector& data, int root) const {
  #ifdef PARALLEL
    // Get copy of data on processor (potentially inefficient)
    Vector data_copy;
    if (root == -1) {
      data_copy = data;
      // Note: Halo data may not be correct if 'data' is different on each processor
    } else {
      data_copy = (mpi.rank==root) ? data : Vector(ndof);
      bcast(data_copy, root);
    }
    // Assign the main blocks
    Vector scattered = assignBlock(data_copy);
    // Assign the halo regions
    for (int i=0; i<mpi.size; i++) {
      for (int j=0; j<priv->nrecv[i]; j++) {
        scattered[priv->irecv[i]+j] = data_copy[priv->recv_lists[i][j]];
      }
    }
    return scattered;

  #else
    return data;
  #endif
  }


  void Communicator::bcast(int& value, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&value, 1, MPI_INT, root, MPI_COMM_WORLD);
  #endif
  }


  void Communicator::bcast(double& value, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  #endif
  }


  void Communicator::bcast(Vector& vector, int root) const {
  #ifdef PARALLEL
    MPI_Bcast(&vector[0], vector.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
  #endif
  }

}
