#ifndef MINIM_RANGE_H
#define MINIM_RANGE_H

#include <vector>

namespace minim {
  namespace range {

    template <typename TOut>
    struct Iterator {
      Iterator(std::vector<int> xSize, std::vector<int> xStart, std::vector<int> xEnd, std::vector<int> xHalo);
      Iterator<TOut> operator()(int i);
      TOut operator*() const;
      Iterator<TOut> operator++();
      bool operator!=(const Iterator& other) const { return i != other.i; };

      private:
        int i, iEnd;
        std::vector<int> x;
        std::vector<int> xSize;
        std::vector<int> xStart;
        std::vector<int> xEnd;
        std::vector<int> xHalo;
    };

  }


  class RangeX {
    using Iterator = range::Iterator<std::vector<int>>;

    public:
      RangeX(std::vector<int> xSize);
      RangeX(std::vector<int> xSize, int halo);
      RangeX(std::vector<int> xSize, std::vector<int> xHalo);

      Iterator begin() { return iter(iStart); }
      Iterator end() { return iter(-1); }

    private:
      int iStart;
      Iterator iter;
  };


  class RangeI {
    using Iterator = range::Iterator<int>;

    public:
      RangeI(std::vector<int> xSize);
      RangeI(std::vector<int> xSize, int halo);
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
