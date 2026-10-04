#ifndef COMMGRID_H
#define COMMGRID_H

#include <vector>
#include <memory>
#include "Communicator.h"

namespace minim {
  using std::vector;

  /// A communicator that distributes grid data over MPI processes.
  ///
  /// The grid is split into a commArray of sub-grids, one per
  /// processor, each surrounded by a halo of the given width. The final
  /// grid dimension is not parallelised.
  class CommGrid : public Communicator {
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
      /// excluding the halo entries.
      double dotProduct(const vector<double>& a, const vector<double>& b) const override;

      // Internal functions
      /// Construct a grid communicator with the given halo width.
      CommGrid(int haloWidth);
      CommGrid(const CommGrid& other);
      CommGrid& operator=(const CommGrid& other);

      std::unique_ptr<Communicator> clone() const override;
      void setup(Potential& pot, size_t ndof, vector<int> ranks) override;

      /// The number of grid dimensions.
      int nDim;
      /// The width of the halo region in grid nodes.
      int haloWidth;
      vector<int> commArray;   // The number of processors along each dimension (nDim)
      vector<int> commIndices; // The array indices of this MPI rank (nDim)
      vector<int> globalSizes; // The total global grid sizes along each dimension (nDim+1)
      vector<int> blockSizes;  // The local grid sizes along each dimension (nDim+1)
      vector<int> procSizes;   // The local grid sizes (including halo) along each dimension (nDim+1)
      vector<int> procStart;   // The index on the global grid where this processor starts (nDim+1)
      vector<int> haloWidths;  // The halo sizes along each dimension (nDim+1)

    private:
      vector<int> getCoords(int loc) const;
      void makeMPITypes() override;

      template<typename T> vector<T> assignBlockImpl(const vector<T>& in) const;
      template<typename T> vector<T> assignProcImpl(const vector<T>& in) const;
  };


  /// A CommGrid with dimension-specific optimisations for 2D grids.
  ///
  /// Used automatically for 2D grids; it specialises the dot product
  /// with an OpenMP-parallel loop collapsed over the two dimensions.
  class CommGrid2 : public CommGrid {
    public:
      /// Construct a 2D grid communicator with the given halo width.
      CommGrid2(int haloWidth);
      std::unique_ptr<Communicator> clone() const override;
      /// Compute the dot product of two vectors over all processors,
      /// excluding the halo entries.
      double dotProduct(const vector<double>& a, const vector<double>& b) const override;
  };

  /// A CommGrid with dimension-specific optimisations for 3D grids.
  ///
  /// Used automatically for 3D grids; it specialises the dot product
  /// with an OpenMP-parallel loop collapsed over the three dimensions.
  class CommGrid3 : public CommGrid {
    public:
      /// Construct a 3D grid communicator with the given halo width.
      CommGrid3(int haloWidth);
      std::unique_ptr<Communicator> clone() const override;
      /// Compute the dot product of two vectors over all processors,
      /// excluding the halo entries.
      double dotProduct(const vector<double>& a, const vector<double>& b) const override;
  };
}

#endif
