#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>
#include <numeric>
#include "Communicator.h"
#include "Potential.h"
#include "minimMpi.h"
#include "vec.h"

using minim::mpi;


class Priv {
  public:
    Priv(int ndof) :
      nblock(mpi.size,ndof/mpi.size),
      iblock(mpi.size),
      nrecv(mpi.size),
      irecv(mpi.size),
      send(mpi.size),
      recv_lists(mpi.size)
    {
      ntot = nblock[mpi.rank];
    };

    int ntot; // Size of block + halo
    std::vector<int> nblock; // Size of each block
    std::vector<int> iblock; // Global index for the start of each block
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
        if (loc >= iblock[i]) return i;
      }
    }
};


Communicator::Communicator(int ndof, Args &args) : priv(new Priv(ndof)) {
  if (mpi.size == 1) return;

  // Get size and index of each block
  for (int i=0; i<mpi.size-1; i++) {
    if (i < ndof % mpi.size) priv->nblock[i]++;
    priv->iblock[i+1] = priv->iblock[i] + priv->nblock[i];
  }

  // Identify halo coordinates based upon the list of energy elements
  std::vector<std::vector<int>> send_lists(mpi.size);
  std::vector<std::vector<int>> blocks(args.elements.size());
  std::vector<std::vector<bool>> in_block(args.elements.size());
  for (int ie=0; ie<args.elements.size(); ie++) {
    Args::Element e = args.elements[ie];
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
          vec::insert_unique(send_lists[blocks[ie][i]], e.idof[j] - priv->iblock[mpi.rank]);
        }
      }
    }
  } //end for elements

  // Store number and starting indicies of coords being recieved from each proc
  priv->irecv[0] = priv->nblock[mpi.rank];
  for (int i=0; i<mpi.size; i++) {
    priv->nrecv[i] = priv->recv_lists[i].size();
    if (i < mpi.size-1) priv->irecv[i+1] = priv->irecv[i] + priv->nrecv[i];
  }
  priv->ntot = priv->irecv[mpi.size-1] + priv->nrecv[mpi.size-1];

  std::vector<int> nelements(mpi.size);
  std::vector<Args::Element> elements_tmp;
  for (int ie=0; ie<args.elements.size(); ie++) {
    // Assign each element to a proc
    int proc = blocks[ie][0];
    int fewest_elements = nelements[proc];
    for (int i : blocks[ie]) {
      if (nelements[i] < fewest_elements) {
        proc = i;
        fewest_elements = nelements[i];
      }
    }
    nelements[proc] ++;

    // Store the elements for this proc
    if (vec::any(in_block[ie])) {
      Args::Element e = args.elements[ie];
      // Update element.idof with local index
      for (int i=0; i<e.idof.size(); i++) {
        if (in_block[ie][i]) {
          e.idof[i] = e.idof[i] - priv->iblock[mpi.rank];
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
#ifdef PARALLEL
    MPI_Type_create_indexed_block(send_lists[i].size(), 1, &send_lists[i][0], MPI_DOUBLE, &priv->sendtype[i]);
    MPI_Type_commit(&priv->sendtype[i]);
#endif
  }
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
  delete priv;
}


std::vector<double> Communicator::assignBlock(std::vector<double> in) {
  std::vector<double> out = std::vector<double>(priv->ntot);
  for (int i=0; i<priv->nblock[mpi.rank]; i++) {
    out[i] = in[priv->iblock[mpi.rank]+i];
  }
  return out;
}


void Communicator::communicate(std::vector<double> &vector) {
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
      int irecv = std::accumulate(priv->nrecv.begin(), priv->nrecv.begin()+i, priv->nblock[mpi.rank]);
      MPI_Recv(&vector[irecv], priv->nrecv[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD, NULL);
    }
  }
#endif
}


double Communicator::get(std::vector<double> vector, int loc) {
  double value;

  if (mpi.size == 1) {
    value = vector[loc];
  } else {
#ifdef PARALLEL
    int i = priv->getBlock(loc);
    value = vector[loc-priv->iblock[i]];
    MPI_Bcast(&value, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif
  }

  return value;
}


double Communicator::dotProduct(std::vector<double> a, std::vector<double> b) {
  double result = std::inner_product(a.begin(), a.begin()+priv->nblock[mpi.rank], b.begin(), 0.0);
#ifdef PARALLEL
  MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return result;
}


std::vector<double> Communicator::gather(std::vector<double> block, int ndof) {
#ifdef PARALLEL
  std::vector<double> gathered(ndof);
  MPI_Allgatherv(&block[0], priv->nblock[mpi.rank], MPI_DOUBLE,
                 &gathered[0], &priv->nblock[0], &priv->iblock[0], MPI_DOUBLE,
                 MPI_COMM_WORLD);
  return gathered;
#else
  return block;
#endif
}


std::vector<double> Communicator::scatter(std::vector<double> &data) {
#ifdef PARALLEL
  std::vector<double> scattered(priv->ntot);
  // Scatter the main blocks
  MPI_Scatterv(&data[0], &priv->nblock[0], &priv->iblock[0], MPI_DOUBLE,
               &scattered[0], priv->nblock[mpi.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Assign the halo regions
  for (int i=0; i<mpi.size; i++) {
    for (int j=0; j<priv->nrecv[i]; j++) {
      scattered[priv->irecv[i]+j] = data[priv->recv_lists[i][j]];
    }
  }
  return scattered;

#else
  return data;
#endif
}
