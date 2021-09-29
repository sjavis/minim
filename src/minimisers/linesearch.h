#ifndef LINESEARCH_H
#define LINESEARCH_H

#include <vector>

namespace minim {
  class State;
  
  double backtrackingLinesearch(State &state, std::vector<double> &step, double de0);
}

#endif
