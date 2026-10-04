#ifndef COMMUNSTRUCTURED_H
#define COMMUNSTRUCTURED_H

#include <vector>
#include <memory>
#include "Communicator.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;
  class Potential;

  /// A communicator that distributes unstructured data over MPI
  /// processes.
  ///
  /// Used by element-based potentials. The degrees of freedom are split
  /// into equal-sized blocks in the order they are given, and the
  /// potential's energy elements are each assigned to the processor
  /// that holds the most of their degrees of freedom. Any remaining
  /// degrees of freedom of an element lie in the halo.
  class CommUnstructured : public Communicator {
    public:
      /// Assign the local block from global or block data.
      vector<int> assignBlock(const vector<int>& in) const override;
      /// Assign the local block from global or block data.
      vector<char> assignBlock(const vector<char>& in) const override;
      /// Assign the local block from global or block data.
      vector<double> assignBlock(const vector<double>& in) const override;

      /// Assign the local processor data (block plus halo) from global data.
      vector<int> assignProc(const vector<int>& in) const override;
      /// Assign the local processor data (block plus halo) from global data.
      vector<char> assignProc(const vector<char>& in) const override;
      /// Assign the local processor data (block plus halo) from global data.
      vector<double> assignProc(const vector<double>& in) const override;

      /// Get the block (processor) that owns a given global index.
      int getBlock(int loc) const override;
      /// Get the local index of a given global index, or -1 if not held
      /// by this processor.
      int getLocalIdx(int loc, int block=-1) const override;

      /// Compute the dot product of two vectors over all processors,
      /// using only the block entries.
      double dotProduct(const vector<double>& a, const vector<double>& b) const override;

      // Internal functions
      /// Construct an unstructured communicator.
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
