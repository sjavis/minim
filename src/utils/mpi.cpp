#include "utils/mpi.h"

#include <iostream>
#include "utils/vec.h"


namespace minim {

  using std::vector;


  void mpiInit() {
    minim::mpi.init(0, 0);
  }

  void mpiInit(int* argc, char*** argv) {
    minim::mpi.init(argc, argv);
  }


  Mpi mpi;

  Mpi::Mpi() {
    size = 1;
    rank = 0;
  }


  void Mpi::init(int* argc, char*** argv) {
#ifdef PARALLEL
    MPI_Init(argc, argv);
    getSizeRank(MPI_COMM_WORLD);
    _init = true;
#else
    std::cerr << "Warning: MPI not initialised. -DPARALLEL flag must be passed to the compiler." << std::endl;
#endif
  }


#ifdef PARALLEL
  void Mpi::getSizeRank(MPI_Comm comm) {
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
  }
#endif


  Mpi::~Mpi() {
#ifdef PARALLEL
    if (_init == true) {
      MPI_Finalize();
    }
#endif
  }


  double Mpi::sum(double a) const {
    double result = a;
  #ifdef PARALLEL
    if (size > 1) MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
    return result;
  }

  double Mpi::sum(const vector<double>& a) const {
    return sum(vec::sum(a));
  }

  double Mpi::dotProduct(const vector<double>& a, const vector<double>& b) const {
    return sum(std::inner_product(a.begin(), a.end(), b.begin(), 0.0));
  }


  void Mpi::bcast(int& value, int root) const {
  #ifdef PARALLEL
    if (size > 1) MPI_Bcast(&value, 1, MPI_INT, root, MPI_COMM_WORLD);
  #endif
  }

  void Mpi::bcast(double& value, int root) const {
  #ifdef PARALLEL
    if (size > 1) MPI_Bcast(&value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  #endif
  }

  void Mpi::bcast(vector<double>& data, int root) const {
  #ifdef PARALLEL
    if (size > 1) {
      int nData = data.size();
      MPI_Bcast(&nData, 1, MPI_INT, root, MPI_COMM_WORLD);
      if (rank != root) data = vector<double>(nData);
      MPI_Bcast(&data[0], nData, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
  #endif
  }

}
