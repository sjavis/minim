#ifndef MINIM_RANGE_H
#define MINIM_RANGE_H

#include <vector>

namespace minim {
  namespace range {

    /// An iterator over a (possibly haloed) multi-dimensional grid.
    template <typename TOut>
    struct Iterator {
      Iterator(std::vector<int> xSize, std::vector<int> xStart, std::vector<int> xEnd, std::vector<int> xHalo);
      Iterator<TOut> operator()(int i);
      /// Get the grid indices of the current position.
      TOut operator*() const;
      Iterator<TOut>& operator++();
      bool operator!=(const Iterator& other) const { return i != other.i; };

      private:
        int i, iEnd, nDim;
        std::vector<int> x;
        std::vector<int> xSize;
        std::vector<int> xStart;
        std::vector<int> xEnd;
        std::vector<int> xHalo;
        std::vector<int> stepSize;
    };

  }


  /// A range over a multi-dimensional grid, yielding the grid indices.
  class RangeX {
    using Iterator = range::Iterator<std::vector<int>>;

    public:
      /// Construct a range over a grid of the given size.
      RangeX(std::vector<int> xSize);
      /// Construct a range over a grid with a halo of the given width.
      RangeX(std::vector<int> xSize, int halo);
      /// Construct a range over a grid with the given halo widths.
      RangeX(std::vector<int> xSize, std::vector<int> xHalo);

      Iterator begin() { return iter(iStart); }
      Iterator end() { return iter(-1); }

    private:
      int iStart;
      Iterator iter;
  };


  /// A range over a multi-dimensional grid, yielding the flat index.
  class RangeI {
    using Iterator = range::Iterator<int>;

    public:
      /// Construct a range over a grid of the given size.
      RangeI(std::vector<int> xSize);
      /// Construct a range over a grid with a halo of the given width.
      RangeI(std::vector<int> xSize, int halo);
      /// Construct a range over a grid with the given halo widths.
      RangeI(std::vector<int> xSize, std::vector<int> xHalo);

      Iterator begin() { return iter(iStart); }
      Iterator end() { return iter(-1); }

    private:
      int iStart;
      Iterator iter;
  };

}

#include "range.hpp"

#endif
