#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include <memory>

namespace minim {
  using std::vector;

  class Communicator {
    public:
      size_t ndof;   //!< Total number of degrees of freedom
      size_t nproc;  //!< Total number of degrees of freedom on processor (including halo)
      size_t nblock; //!< Number of degrees of freedom assigned to processor (excluding halo)
      int iblock;    //!< The starting index for this processor

      bool usesThisProc = true;
      vector<int> ranks = vector<int>();

      int rank() const;
      int size() const;

      virtual vector<int> assignBlock(const vector<int>& in) const = 0;
      virtual vector<bool> assignBlock(const vector<bool>& in) const = 0;
      virtual vector<double> assignBlock(const vector<double>& in) const = 0;

      virtual vector<int> assignProc(const vector<int>& in) const = 0;
      virtual vector<bool> assignProc(const vector<bool>& in) const = 0;
      virtual vector<double> assignProc(const vector<double>& in) const = 0;

      double get(const vector<double>& vector, int i) const;

      void communicate(vector<double>& vector) const;
      vector<double> gather(const vector<double>& block, int root=-1) const;
      vector<double> scatter(const vector<double>& data, int root=-1) const;

      void bcast(int& value, int root=0) const;
      void bcast(double& value, int root=0) const;
      void bcast(vector<double>& value, int root=0) const;

      double sum(double a) const;
      double sum(const vector<double>& a) const;
      double norm(const vector<double>& a) const;
      double dotProduct(const vector<double>& a, const vector<double>& b) const;

      ~Communicator();
      std::unique_ptr<Communicator> clone() const = 0;
      virtual void setup(Potential& pot, size_t ndof, vector<int> ranks) = 0;

    protected:
      void defaultSetup(size_t ndof, vector<int> ranks);

    private:
      int commSize = 1;
      int commRank = 0;
      vector<bool> send;  // States if data is to be sent to each proc
      vector<bool> recv;  // States if data is to be received from each proc

      #ifdef PARALLEL
      MPI_Comm comm;
      vector<MPI_Datatype> sendtype; // MPI derived datatype to send to each proc
      vector<MPI_Datatype> recvtype; // MPI derived datatype to recieve from each proc
      MPI_Datatype blockType;        // MPI derived datatype to send the local block
      MPI_Datatype gatherType;       // MPI derived datatype to receive the blocks for gathering
      #endif

      void setComm(vector<int> ranks);
      virtual void makeMPITypes() = 0;
  };

}

#endif
