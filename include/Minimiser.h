#ifndef MINIMISER_H
#define MINIMISER_H

#include <vector>
#include <memory>
#include <functional>

namespace minim {
  class State;

  // Abstract class for minimisation proceedures
  class Minimiser {
    public:
      int maxIter = 100000;

      typedef void (*AdjustFunc)(int, State&);
      int iter;

      Minimiser() {};
      virtual ~Minimiser() {};
      virtual std::unique_ptr<Minimiser> clone() const = 0;

      virtual Minimiser& setMaxIter(int maxIter);

      std::vector<double> minimise(State& state, std::function<void(int,State&)> adjustState=nullptr);
      virtual void init(State& state) {};
      virtual void iteration(State& state) = 0;
      virtual bool checkConvergence(const State& state) { return false; };
  };


  // An intermediate class is used to return the derived type for methods that return a Minimiser
  template<typename Derived>
  class NewMinimiser : public Minimiser {
    public:
      std::unique_ptr<Minimiser> clone() const override {
        return std::make_unique<Derived>(static_cast<const Derived&>(*this));
      }
  };

}

#endif
