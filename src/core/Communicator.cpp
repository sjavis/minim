#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>
#include <numeric>
#include "Communicator.h"
#include "minimMpi.h"

using minim::mpi;


Communicator::Communicator(int ndof)
  : _ndof(ndof), _nblock(mpi.size,ndof/mpi.size), _nrecv(mpi.size), _isend(mpi.size) {
  // Get size of each block
  for (int i=0; i<mpi.size-1; i++) {
    if (i < ndof % mpi.size) _nblock[i]++;
  }
}

//Communicator::Communicator(int ndof)
//  : _ndof(ndof), _nblock(mpi.size,ndof/mpi.size),
//    _nrecv(mpi.size), _isend(mpi.size), _sendtype(mpi.size) {
//#ifdef PARALLEL
//
//  // Get size of each block
//  for (int i=0; i<mpi.size-1; i++) {
//    if (i < ndof % mpi.size) _nblock[i]++;
//  }
//
//  // Get indices of block to send to each halo
//  for (int i=0; i<mpi.size; i++) {
//    if (i == mpi.rank) continue;
//    _isend[i] = {};
//
//    MPI_indexed(_isend[i], sendtype[i]);
//    MPI_commit(sendtype[i]);
//  }
//
//  // Get number being received from each block
//  for (int i=0; i<mpi.size; i++) {
//    if (i == mpi.rank) continue;
//    nrecv[i] = ...;
//  }
//#endif
//}


Communicator::~Communicator() {
#ifdef PARALLEL
  //for (i /= mpi.rank) {
  //  // Free any committed MPI datatypes
  //  MPI_free(sendtype[i]);
  //}
#endif
}


std::vector<double> Communicator::assignBlock(std::vector<double> in) {
  int istart = std::accumulate(_nblock.begin(), _nblock.begin()+mpi.rank, 0);
  return std::vector<double>(in.begin()+istart, in.begin()+istart+_nblock[mpi.rank]);
}


void Communicator::communicate(std::vector<double> &vector) {
#ifdef PARALLEL
  for (int i=0; i<mpi.size; i++) {
    if (i == mpi.rank) continue;
    // Send
    if (!_isend[i].empty()) {
      int tag = mpi.rank*mpi.size + i;
      MPI_Send(&vector[0], 1, _sendtype[i], i, tag, MPI_COMM_WORLD);
    }
    // Receive
    if (_nrecv[i] > 0) {
      int tag = i*mpi.size + mpi.rank;
      int irecv = std::accumulate(_nrecv.begin(), _nrecv.begin()+i, _nblock[mpi.rank]);
      MPI_Recv(&vector[irecv], _nrecv[i], MPI_DOUBLE, i, tag, MPI_COMM_WORLD, NULL);
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
    int i0 = _ndof;
    for (int i=mpi.size-1; i>=0; i--) {
      i0 -= _nblock[i];
      if (loc >= i0) {
        value = vector[loc-i0];
        MPI_Bcast(&value, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
        break;
      }
    }
#endif
  }

  return value;
}


double Communicator::dotProduct(std::vector<double> a, std::vector<double> b) {
  double result = std::inner_product(a.begin(), a.begin()+_nblock[mpi.rank], b.begin(), 0.0);
#ifdef PARALLEL
  MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return result;
}


std::vector<double> Communicator::gather(std::vector<double> block) {
#ifdef PARALLEL
  std::vector<double> gathered(_ndof);
  int dspls[mpi.size];
  std::partial_sum(_nblock.begin(), _nblock.end(), dspls);

  MPI_Allgatherv(&block[0], _nblock[mpi.rank], MPI_DOUBLE,
                 &gathered[0], &_nblock[0], dspls, MPI_DOUBLE,
                 MPI_COMM_WORLD);

  return gathered;

#else
  return block;
#endif
}
