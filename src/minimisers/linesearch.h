#ifndef LINESEARCH_H
#define LINESEARCH_H

class State;

double backtrackingLinesearch(State &state, std::vector<double> &step, double de0);

#endif
