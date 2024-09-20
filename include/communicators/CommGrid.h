#ifndef COMMGRID_H
#define COMMGRID_H

#include <vector>
#include <memory>
#include "Communicator.h"

namespace minim {
  using std::vector;

  class CommGrid : public Communicator {
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

      // Internal functions
      CommGrid() = default;
      CommGrid(const CommGrid& other);
      CommGrid& operator=(const CommGrid& other);

      std::unique_ptr<Communicator> clone() const override;
      void setup(Potential& pot, size_t ndof, vector<int> ranks) override;

      int nDim;
      int haloWidth = 1;
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

}

#endif
