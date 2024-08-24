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
      vector<bool> assignBlock(const vector<bool>& in) const override;
      vector<double> assignBlock(const vector<double>& in) const override;

      vector<int> assignProc(const vector<int>& in) const override;
      vector<bool> assignProc(const vector<bool>& in) const override;
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
      vector<int> commArray;   // The number of processors along each dimension
      vector<int> commIndices; // The array indices of this MPI rank
      vector<int> globalSizes; // The total global grid sizes along each dimension
      vector<int> blockSizes;
      vector<int> blockStart;//
      vector<int> procSizes;
      vector<int> procStart;

    private:
      vector<int> getCoords(int loc) const;
      void makeMPITypes() override;

      template<typename T> vector<T> assignBlockImpl(const vector<T>& in) const;
      template<typename T> vector<T> assignProcImpl(const vector<T>& in) const;
  };

}

#endif
