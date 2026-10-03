#ifndef MINIMISER_H
#define MINIMISER_H

#include <vector>
#include <memory>
#include <string>
#include <functional>

namespace minim {
  class State;

  /// Abstract base class for minimisation procedures.
  ///
  /// Holds the settings common to different minimisation algorithms and
  /// provides the ``minimise`` entry point.
  class Minimiser {
    public:
      /// Maximum number of iterations before terminating.
      int maxIter = 100000;
      /// The line search method to use.
      std::string linesearch = "backtracking";

      typedef void (*AdjustFunc)(int, State&);
      int iter;

      virtual ~Minimiser() = default;
      /// Create a deep copy of the minimiser.
      virtual std::unique_ptr<Minimiser> clone() const = 0;

      /// Set the maximum number of iterations.
      virtual Minimiser& setMaxIter(int maxIter);
      /// Set the line search method.
      Minimiser& setLinesearch(std::string method);

      /// Minimise a state with an optional function to be run each iteration.
      ///
      /// @param state The state to be minimised.
      /// @param adjustState Optional function run each iteration.
      /// @return The final coordinates.
      std::vector<double> minimise(State& state, std::function<void(int,State&)> adjustState=nullptr);
      /// Minimise with a predefined log function.
      ///
      /// Format: [fields]-[iter], e.g. "e-100" logs the energy every 100 iterations.
      ///
      /// @param state The state to be minimised.
      /// @param logType The log format.
      /// @return The final coordinates.
      std::vector<double> minimise(State& state, std::string logType);

      /// Initialise the minimiser before the first iteration.
      virtual void init(State& state) {};
      /// Perform a single minimisation iteration.
      virtual void iteration(State& state) = 0;
      /// Check whether the minimisation is converged.
      virtual bool checkConvergence(const State& state) { return false; };
  };


  /// An intermediate class used to return the derived type for methods that return a Minimiser.
  template<typename Derived>
  class NewMinimiser : public Minimiser {
    public:
      std::unique_ptr<Minimiser> clone() const override {
        return std::make_unique<Derived>(static_cast<const Derived&>(*this));
      }
  };

}

#endif
