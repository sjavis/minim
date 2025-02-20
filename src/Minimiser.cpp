#include "Minimiser.h"

#include <stdexcept>
#include <sstream>
#include <iomanip>
#include "State.h"
#include "utils/vec.h"
#include "utils/mpi.h"
#include "utils/print.h"

namespace minim {

  Minimiser& Minimiser::setMaxIter(int maxIter) {
    this->maxIter = maxIter;
    return *this;
  }

  Minimiser& Minimiser::setLinesearch(std::string method) {
    if (!vec::isIn({"backtracking","none"}, method)) {
      throw std::invalid_argument("Invalid line search method.");
    }
    linesearch = method;
    return *this;
  }


  std::vector<double> Minimiser::minimise(State& state, std::function<void(int,State&)> adjustState) {
    if (!state.usesThisProc) return std::vector<double>();

    init(state);
    for (iter=0; iter<=maxIter; iter++) {
      if (adjustState) adjustState(iter, state);
      if (state.isFailed) break;
      iteration(state);
      if (checkConvergence(state)) break;
    }
    return state.coords();
  }


  std::vector<double> Minimiser::minimise(State& state, std::string logType) {
    std::function<void(int,State&)> logFn = nullptr;

    // Check if logType in the correct format, [fields]-[iteration]
    size_t divPos = logType.find("-");
    if (divPos < logType.length()) {
      std::string logFields = logType.substr(0, divPos);
      int logIter = std::stoi(logType.substr(divPos+1));
      if (logIter > 0) {
        int iterDigits = std::to_string(maxIter).length();
        // Define the log function
        logFn = [logIter, logFields, iterDigits](int i, State& s){
          if (i%logIter!=0) return;
          bool logE = (logFields.find('e') < logFields.length());
          bool logG = (logFields.find('g') < logFields.length());
          double e;
          vector<double> g;
          s.energyGradient(logE?(&e):nullptr, logG?(&g):nullptr);
          std::ostringstream log;
          log << "I: " << std::setw(iterDigits) << i;
          if (logE) log << "  E: " << e;
          if (logG) log << "  G: " << vec::rms(g);
          print(log.str());
        };
      }
    }

    return minimise(state, logFn);
  }

}
