#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>


class Communicator {
  public:
    Communicator(int ndof);
    ~Communicator();

    std::vector<double> assignBlock(std::vector<double> in);

    void communicate(std::vector<double> &vector);
    double get(std::vector<double> vector, int i);
    double dotProduct(std::vector<double> a, std::vector<double> b);
    std::vector<double> gather(std::vector<double> block);
    
  private:
    int _ndof;
    std::vector<int> _nblock;             // Size of each block
    std::vector<int> _nrecv;              // Number of halo coordinates to recieve from each proc
    std::vector<std::vector<int>> _isend; // List of block indicies to send to each proc
#ifdef PARALLEL
    std::vector<MPI_Datatype> _sendtype;  // MPI derived datatype to send to each proc
#endif
};

#endif
