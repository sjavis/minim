#ifndef LINESEARCH_H
#define LINESEARCH_H

class State;

void backtrackingLinesearch(State &state, std::vector<double> &step, double slope);

#endif
