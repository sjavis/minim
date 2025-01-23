#ifndef COMMUNSTRUCTURED_H
#define COMMUNSTRUCTURED_H

#include <vector>
#include <memory>
#include "Communicator.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;
  class Potential;

  class CommUnstructured : public Communicator {
    public:
      // Assign data
      vector<int> assignBlock(const vector<int>& in) const override;
      vector<char> assignBlock(const vector<char>& in) const override;
      vector<double> assignBlock(const vector<double>& in) const override;

      vector<int> assignProc(const vector<int>& in) const override;
      vector<char> assignProc(const vector<char>& in) const override;
      vector<double> assignProc(const vector<double>& in) const override;

      // Access data
      int getBlock(int loc) const override;
      int getLocalIdx(int loc, int block=-1) const override;

      // MPI reduction functions
      double dotProduct(const vector<double>& a, const vector<double>& b) const override;

      // Internal functions
      CommUnstructured();
      CommUnstructured(const CommUnstructured& other);
      CommUnstructured& operator=(const CommUnstructured& other);

      std::unique_ptr<Communicator> clone() const override;
      void setup(Potential& pot, size_t ndof, vector<int> ranks) override;

    private:
      int iblock;    // The starting index for this processor
      vector<int> nblocks; // Size of each block
      vector<int> iblocks; // Global index for the start of each block
      vector<int> nrecv;  // Number of halo coordinates to recieve from each proc
      vector<int> irecv;  // Starting indicies for each proc in halo
      vector2d<int> recv_lists; // List of indicies to recieve from each proc
      vector2d<int> send_lists; // List of block indicies to send to each other proc

      void getElementBlocks(const Potential& pot, vector2d<int>& blocks, vector2d<char>& in_block);
      void setCommLists(const Potential& pot, const vector2d<int>& blocks, const vector2d<char>& in_block);
      void setRecvSizes();
      vector<int> assignElements(int nElements, const vector2d<int>& blocks);
      void distributeElements(Potential& pot, const vector2d<int>& blocks, const vector2d<char>& in_block);
      bool checkWellDistributed(int ndof);

      void makeMPITypes() override;

      template<typename T> vector<T> assignBlockImpl(const vector<T>& in) const;
      template<typename T> vector<T> assignProcImpl(const vector<T>& in) const;
  };

}

#endif
